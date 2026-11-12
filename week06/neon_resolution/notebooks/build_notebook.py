"""Notebook authoring only; no scientific packages imported on Windows."""
import json
from pathlib import Path
import ast

HERE = Path(__file__).resolve().parent
cells = []
def md(s):
    cells.append(dict(cell_type='markdown', metadata={}, source=s.strip().splitlines(True)))
def code(s):
    cells.append(dict(cell_type='code', metadata={}, execution_count=None, outputs=[], source=s.strip().splitlines(True)))
def section(title, explanation, src):
    md('## '+title+'\n\n'+explanation)
    code(src)

md('''# NEON GRSM·ORNL 2021: 나무 탐지 → 수고·수관 → DBH → AGB·탄소

**실행 환경: Google Colab. Windows에서는 실행하지 않습니다.**

실행 순서: `런타임 → 런타임 유형 변경 → T4 GPU`를 선택하고 위에서 아래로 실행합니다.
첫 실행은 설치 5–15분, 실측·계수 다운로드 5–20분, RGB/CHM 다운로드 20분–수 시간,
탐지·비교 10–60분 정도를 예상합니다. 이는 계획용 추정이며 타일 수·Drive 속도·GPU에 따라 달라집니다.
10cm RGB 한 타일은 압축 전 약 300MB이므로 Drive 여유 공간 10–30GB 이상을 권장합니다.
다운로드 전에 예상 타일 수를 출력합니다. 작은 실행은 `max_plots_per_site=2`로 시작하세요.

NEON API 토큰이 필요할 수 있습니다. [NEON 계정](https://data.neonscience.org/)의 API Tokens에서 발급하고
Colab 보안 비밀 `NEON_TOKEN`에 등록하거나 숨김 입력창에 입력하세요. 토큰은 Drive에 저장하지 않습니다.

데이터·결과는 `/content/drive/MyDrive/neon_grsm_ornl/`에 저장합니다.
HF 모델 캐시는 Drive와 분리된 `/content/hf_cache`이며 새 런타임에서는 다시 받을 수 있습니다.
완료된 다운로드와 매개변수가 같은 중간 결과는 재사용합니다. 설정이 바뀌면 별도 실행 폴더에 결과를 만듭니다.

**비교 대상의 범위:** 기본은 2021년 실측, 살아 있는 DBH≥10cm의 나무입니다.
작은 나무의 중첩 소구획 조사와 전체 조사 면적을 혼합하지 않기 위한 공통 비교 범위입니다.
실측 좌표가 없는 개체도 조사구 참조 총량에 포함합니다. 좌표는 개체 매칭에만 필요합니다.
ORNL의 2021-06 RGB/CHM은 API에서 PROVISIONAL로 제공되므로 이를 명시적으로 포함합니다.
NEONForestAGB는 실측 DBH로 계산한 참조 추정값이며 직접 무게를 잰 정답이 아닙니다.
탐지 정확도는 수관의 정답 윤곽이 아닌 **줄기 위치 기반 대리지표**입니다.
모델은 NEON 사전학습 가중치를 그대로 사용하므로 완전히 독립적인 외부 지역 검증은 아닙니다.
''')
section('01. 패키지 설치', '첫 코드 셀입니다. 설치 후 Colab이 런타임 재시작을 요청하면 재시작하고 다음 셀부터 실행하세요.', r'''
%pip -q install "deepforest>=1.5,<3" "neonutilities>=2.0.1" rasterio geopandas shapely pyproj scikit-learn matplotlib scipy scikit-image requests tabulate
''')
section('02. Drive 마운트와 재사용 설정', '''한국 자료에는 아래 입력·공식 설정과 자료 로더를 바꾸고, 뒤의 탐지·추정·매칭·평가 함수를 재사용합니다.
`ALLOW_ADJACENT_YEARS=False`가 기본입니다. 켜면 2021년 조사구를 유지하고 2021년 없는 조사구만
2020→2022 순서로 추가합니다. 서로 다른 해의 개체를 한 조사구에 섞지 않습니다.
보완 자료는 2021 결과와 분리 집계합니다.''', r'''
from google.colab import drive
drive.mount('/content/drive')
from pathlib import Path
import os, json, hashlib, inspect, re, time, warnings, importlib.metadata, getpass, tempfile

ROOT = Path('/content/drive/MyDrive/neon_grsm_ornl')
CFG = dict(
    sites=['GRSM','ORNL'], year=2021, allow_adjacent_years=False,
    # 실측 연도 보완과 별개인 영상 대체 정책. 2021 공통 자료가 있으면 항상 우선합니다.
    imagery_allow_adjacent_years={'GRSM':False,'ORNL':True},
    imagery_adjacent_priority=[2022,2020],
    aop_include_provisional=True, api_refresh=True,
    adjacent_priority=[2020,2022], min_plots_per_site=5,
    max_plots_per_site=None, # 2로 바꾸면 소규모 점검
    release='current', buffer_m=50, patch_size=400, patch_overlap=0.25,
    score_threshold=0.15, nms_iou=0.15, min_height_m=3.0,
    expected_plot_count=23,
    comparison_plot_ids=['GRSM_001','GRSM_009','GRSM_012','GRSM_014','GRSM_016','GRSM_017',
        'GRSM_047','GRSM_049','GRSM_050','GRSM_051','GRSM_055','GRSM_057','GRSM_058','GRSM_059','GRSM_060',
        'ORNL_001','ORNL_004','ORNL_006','ORNL_009','ORNL_010','ORNL_027','ORNL_049','ORNL_057'],
    search_score_thresholds=[0.1,0.2,0.3,0.4,0.5],search_patch_sizes=[300,400,600],
    search_inner_folds=4,chm_sigma_m=1.0,chm_window_a=2.51,chm_window_b=0.00901,
    fusion_dedup_m=1.0,default_method_metric='carbon_RMSE_tC_ha',
    comparison_overlay_plots=['GRSM_016','GRSM_051','GRSM_012','ORNL_010'],
    min_field_dbh_cm=10.0, visible_height_ratio=0.8, visible_radius_m=3.0,
    match_size='stemDiameter', # 가장 큰 DBH 우선, 동률 때 중심거리
    jucker=dict(a=0.367, b=0.887, g=0.0, biome=11),
    # itcSegment/R/dbh.R biome 11: Nearctic-Temperate mixed forests-Angiosperm
    agb_default=dict(beta0=-2.2118, beta1=2.4133,
                     group='Chojnacky general hardwood surrogate (Cor/Eri/Lau/etc)'),
    carbon_fraction=0.5, carbon_alternative=0.47, cv_folds=5,
    model_name='weecology/deepforest-tree', model_revision='main',
    # 한국 자료용: site/plot/year/area_m2와 미터 투영좌표 polygon이 있는 GPKG/GeoJSON
    footprint_override=None,
    local_rgb={}, local_chm={}, # 예: {'KOREA_plot1':'/content/drive/.../rgb.tif'}
)
for d in ['data/raw','data/cache','data/crops','models','outputs']:
    (ROOT/d).mkdir(parents=True, exist_ok=True)
os.chdir(ROOT/'data/raw')
tempfile.tempdir = str(ROOT/'data/cache')
os.environ['HF_HOME'] = '/content/hf_cache'
os.environ['HF_HUB_CACHE'] = '/content/hf_cache/hub'
Path('/content/hf_cache').mkdir(parents=True,exist_ok=True)
os.environ['TORCH_HOME'] = str(ROOT/'models/torch')
os.environ['MPLCONFIGDIR'] = str(ROOT/'data/cache/matplotlib')
RUN_ID = hashlib.sha256(json.dumps(dict(config=CFG,versions={p:importlib.metadata.version(p)
    for p in ['deepforest','neonutilities','rasterio','geopandas']}),sort_keys=True).encode()).hexdigest()[:12]
OUT = ROOT/'outputs'/RUN_ID
OUT.mkdir(exist_ok=True)
(OUT/'config.json').write_text(json.dumps(CFG,ensure_ascii=False,indent=2))
try:
    from google.colab import userdata
    TOKEN = userdata.get('NEON_TOKEN')
except Exception:
    TOKEN = os.environ.get('NEON_TOKEN') or getpass.getpass('NEON API token (숨김 입력): ')
os.environ['NEON_TOKEN'] = TOKEN
print('결과 폴더:', OUT)
''')
section('03. 공통 라이브러리·캐시·버전 기록', '파일은 임시 파일로 받은 뒤 완료 시 이름을 바꿉니다. 위치 API의 성공 응답을 캐시하고 인증 오류를 숨기지 않습니다.', r'''
import math
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.windows import from_bounds, Window
from rasterio.plot import show
from shapely.geometry import Point, box, Polygon
from shapely.ops import unary_union
from pyproj import CRS, Transformer
from scipy.optimize import linear_sum_assignment
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import neonutilities as nu

session = requests.Session()
session.mount('https://', HTTPAdapter(max_retries=Retry(total=5,backoff_factor=1,
    status_forcelist=[429,500,502,503,504])))
VERSIONS = {p:importlib.metadata.version(p) for p in
    ['deepforest','neonutilities','numpy','pandas','rasterio','geopandas','shapely','pyproj','scikit-learn']}
(OUT/'versions.json').write_text(json.dumps(VERSIONS,indent=2))
AUDIT = []
def audit(stage, item, reason):
    AUDIT.append(dict(stage=stage,item=str(item),reason=str(reason)))

def cached_download(url, dest):
    dest=Path(dest)
    if dest.exists() and dest.stat().st_size>0: return dest
    tmp=dest.with_suffix(dest.suffix+'.part')
    with session.get(url,stream=True,timeout=180) as r:
        r.raise_for_status()
        with tmp.open('wb') as f:
            for block in r.iter_content(1024*1024): f.write(block)
        # requests가 gzip을 해제한 길이는 전송된 Content-Length와 다릅니다.
        if not r.headers.get('Content-Encoding') and r.headers.get('Content-Length') and tmp.stat().st_size != int(r.headers['Content-Length']):
            raise IOError('다운로드 크기 불일치: '+url)
    tmp.replace(dest)
    return dest

def cached_location(name):
    key=hashlib.sha256(name.encode()).hexdigest()
    path=ROOT/'data/cache'/f'location_{key}.json'
    if path.exists(): return json.loads(path.read_text())
    url='https://data.neonscience.org/api/v0/locations/'+requests.utils.quote(name,safe='._')
    r=session.get(url,headers={'X-API-Token':TOKEN},timeout=60)
    r.raise_for_status()
    data=r.json()['data']
    path.write_text(json.dumps(data,ensure_ascii=False))
    return data

def location_xy(d):
    zone=int(d['locationUtmZone'])
    south=d.get('locationUtmHemisphere','N')=='S'
    return float(d['locationUtmEasting']),float(d['locationUtmNorthing']),CRS.from_epsg((32700 if south else 32600)+zone)

def require_columns(df, cols, label):
    missing=set(cols)-set(df.columns)
    if missing: raise ValueError(f'{label} 필드 누락: {missing}; 실제 필드={list(df.columns)}')

def choose_table(tables, name):
    keys=[k for k in tables if k==name or k.endswith(name)]
    if len(keys)!=1: raise KeyError(f'{name}: {keys}; tables={list(tables)}')
    return tables[keys[0]].copy()

print(VERSIONS)
print('load_by_product',inspect.signature(nu.load_by_product))
print('by_tile_aop',inspect.signature(nu.by_tile_aop))
''')
section('04. 실측 자료 다운로드', '''공식 `load_by_product`로 연도별 자료를 받고 표를 Drive에 저장합니다.
2024년 이후 배포본은 과거 자료의 subplotID도 `31_100`, `21_400` 같은 형식으로 수정했습니다.
DBH≥10cm 복수 줄기는 NEON에서 개별 individualID가 부여되므로 서로 합치지 않습니다.
같은 ID의 중복 방문 기록은 마지막 측정만 사용합니다.''', r'''
TABLES=['vst_apparentindividual','vst_mappingandtagging','vst_perplotperyear']
def load_field_year(site, year):
    folder=ROOT/'data/raw'/f'vst_{site}_{year}_{CFG["release"]}'
    folder.mkdir(exist_ok=True)
    done=folder/'complete.json'
    if done.exists():
        return {t:pd.read_csv(folder/f'{t}.csv',low_memory=False) for t in TABLES}
    # load_by_product는 savepath가 없으며 사전 캐시 확인 후 실행합니다.
    tables=nu.load_by_product(dpid='DP1.10098.001',site=[site],
        startdate=f'{year}-01',enddate=f'{year}-12',release=CFG['release'],
        package='expanded',check_size=False,token=TOKEN)
    result={t:choose_table(tables,t) for t in TABLES}
    for t,df in result.items(): df.to_csv(folder/f'{t}.csv',index=False)
    # 변수·품질·인용 정보도 보존
    for k,v in tables.items():
        if isinstance(v,pd.DataFrame) and k not in TABLES:
            v.to_csv(folder/(re.sub(r'[^\w.-]','_',k)+'.csv'),index=False)
    done.write_text(json.dumps({'site':site,'year':year,'keys':list(tables)}))
    return result

def prepare_field(tables, site, year):
    a=tables[TABLES[0]].copy(); m=tables[TABLES[1]].copy(); p=tables[TABLES[2]].copy()
    require_columns(a,['individualID','plotID','date','stemDiameter','height','plantStatus','subplotID'],TABLES[0])
    require_columns(m,['individualID','pointID','stemDistance','stemAzimuth','taxonID'],TABLES[1])
    require_columns(p,['plotID','date','totalSampledAreaTrees'],TABLES[2])
    a=a[pd.to_datetime(a.date,errors='coerce').dt.year.eq(year)].copy()
    p=p[pd.to_datetime(p.date,errors='coerce').dt.year.eq(year)].copy()
    if 'dataCollected' in p:
        p=p[p.dataCollected.fillna('').str.lower().eq('allgrowthforms')]
    if 'samplingImpractical' in p:
        p=p[p.samplingImpractical.isna() | p.samplingImpractical.str.lower().eq('ok')]
    p=p.sort_values('date').drop_duplicates('plotID',keep='last')
    # 데이터 수집이 완료된 가장 최근 event만 사용
    if 'eventID' in a and 'eventID' in p:
        a=a.merge(p[['plotID','eventID']].drop_duplicates(),on=['plotID','eventID'],how='inner')
    else:
        a=a[a.plotID.isin(p.plotID)]
    for c in ['stemDiameter','height','maxCrownDiameter','ninetyCrownDiameter']:
        if c not in a: a[c]=np.nan
        a[c]=pd.to_numeric(a[c],errors='coerce')
    for c in ['stemDistance','stemAzimuth']:
        m[c]=pd.to_numeric(m[c],errors='coerce')
    # mapping에 여러 이력이 있으면 관측연도까지의 가장 최근 행
    if 'date' in m:
        m['_mapping_date']=pd.to_datetime(m.date,errors='coerce')
        m=m[m._mapping_date.isna() | m._mapping_date.dt.year.le(year)]
        m=m.sort_values('_mapping_date',na_position='first')
    key=['individualID']
    if 'stemID' in m and 'stemID' in a: key.append('stemID')
    m=m.drop_duplicates(key,keep='last')
    # apparent 표의 정체성 필드는 보존하고 위치·분류 필드만 보강
    mergecols=key+[c for c in ['pointID','stemDistance','stemAzimuth','taxonID','namedLocation'] if c in m and c not in key]
    for c in mergecols:
        if c in a and c not in key and c not in ['taxonID','namedLocation']: a=a.drop(columns=c)
    a=a.merge(m[mergecols],on=key,how='left',suffixes=('','_map'),validate='many_to_one')
    for c in ['taxonID','namedLocation']:
        if c+'_map' in a: a[c]=a[c].fillna(a[c+'_map']) if c in a else a[c+'_map']
    a['siteID']=site; a['field_year']=year
    a['plotNamedLocation']=a['plotID']+'.basePlot.vst'
    a['field_CD_m']=(a.maxCrownDiameter+a.ninetyCrownDiameter)/2
    # footprint 산정용: 살아 있는 나무로 제한하기 전의 모든 조사 나무
    footprint_evidence=a.copy()
    live=a.plantStatus.fillna('').str.startswith('Live')
    treeform=a.get('growthForm',pd.Series('tree',index=a.index)).fillna('').str.contains('tree',case=False)
    a=a[live & treeform & a.stemDiameter.ge(CFG['min_field_dbh_cm'])].copy()
    rows=[]
    for individual,g in a.groupby('individualID'):
        # 큰 나무의 multi-bole은 ID별로 구분된다. 같은 ID를 다시 합산하지 않는다.
        g=g.sort_values('date')
        r=g.iloc[-1].copy()
        r['duplicate_measurements']=len(g)
        r['field_DBH_cm']=r['stemDiameter']
        r['field_H_m']=r['height']
        rows.append(r)
    return pd.DataFrame(rows),p,footprint_evidence

FIELD_PARTS=[]; PLOT_PARTS=[]; EVIDENCE_PARTS=[]
for site in CFG['sites']:
    years=[CFG['year']]+(CFG['adjacent_priority'] if CFG['allow_adjacent_years'] else [])
    selected=set()
    for year in years:
        if year!=CFG['year'] and len(selected)>=CFG['min_plots_per_site']: break
        tables=load_field_year(site,year)
        f,p,e=prepare_field(tables,site,year)
        available=[x for x in p.plotID if x not in selected and x in set(f.plotID)]
        if year!=CFG['year']:
            available=available[:max(0,CFG['min_plots_per_site']-len(selected))]
        if CFG['max_plots_per_site'] is not None:
            available=available[:max(0,CFG['max_plots_per_site']-len(selected))]
        selected.update(available)
        FIELD_PARTS.append(f[f.plotID.isin(available)])
        PLOT_PARTS.append(p[p.plotID.isin(available)].assign(siteID=site,field_year=year))
        EVIDENCE_PARTS.append(e[e.plotID.isin(available)])
field=pd.concat(FIELD_PARTS,ignore_index=True)
plotmeta=pd.concat(PLOT_PARTS,ignore_index=True)
evidence=pd.concat(EVIDENCE_PARTS,ignore_index=True)
if field.empty: raise RuntimeError('해당 연도 실측 나무 없음. 인접연도 옵션/자료 가용성 확인')
display(field.groupby(['siteID','field_year']).agg(trees=('individualID','size'),plots=('plotID','nunique')))
''')
section('05. 기준점 좌표와 실제 조사 경계 복원', '''사용자 안내서 §3.1, §3.5, §6.1을 따릅니다.
동쪽 = 기준점 E + 거리×sin(방위각), 북쪽 = 기준점 N + 거리×cos(방위각).
큰 나무가 기록된 `100/400㎡ subplotID`의 API 중심점으로 정사각형 조사 범위를 합칩니다.
살아 있는 나무가 없어도 기록된 고사목을 경계 증거에 포함합니다. 면적은 perplot의
`totalSampledAreaTrees`와 검산합니다. **빈 소구획을 추측해서 확장하지 않습니다.**
모든 소구획을 복원하지 못하면 해당 조사구를 제외하고 진단 CSV에 남깁니다.
그 경우 정확한 조사 경계를 별도 파일로 제공하는 `footprint_override`를 사용할 수 있습니다.
경계의 소구획 모서리는 정북 기준이며 실제 경사면 면적 대신 수평 투영면적을 사용합니다.''', r'''
def offset_xy(e,n,distance,azimuth_deg):
    rad=np.deg2rad(azimuth_deg)
    return e+distance*np.sin(rad), n+distance*np.cos(rad)

def map_individuals(df):
    out=[]
    for _,r in df.iterrows():
        try:
            if not np.isfinite(r.stemDistance) or not np.isfinite(r.stemAzimuth) or pd.isna(r.pointID):
                raise ValueError('정확한 stem offset 없음')
            pid=str(r.pointID).removesuffix('.0')
            loc=cached_location(r.plotNamedLocation+'.'+pid)
            e,n,crs=location_xy(loc)
            x,y=offset_xy(e,n,r.stemDistance,r.stemAzimuth)
            rr=r.to_dict(); rr.update(x=x,y=y,source_epsg=crs.to_epsg())
            # 사이트마다 CRS가 달라도 공통 작업 CRS로 변환
            xx,yy=Transformer.from_crs(crs,WORK_CRS,always_xy=True).transform(x,y)
            rr['geometry']=Point(xx,yy)
            out.append(rr)
        except requests.HTTPError as exc:
            if exc.response.status_code in [401,403,429]: raise
            audit('coordinates',r.individualID,exc)
        except (ValueError,KeyError,TypeError) as exc: audit('coordinates',r.individualID,exc)
    if not out:
        pd.DataFrame(AUDIT).to_csv(OUT/'audit.csv',index=False)
        raise RuntimeError('매핑된 실측 나무가 없습니다. pointID·API 진단은 audit.csv 확인')
    return gpd.GeoDataFrame(out,geometry='geometry',crs=WORK_CRS)

WORK_CRS=CRS.from_epsg(32617) # GRSM·ORNL; 한국에는 현지 미터 투영 CRS 사용
field_geo=map_individuals(field)

def build_footprints(meta,evidence):
    plots=[]
    override=None
    if CFG['footprint_override']:
        override=gpd.read_file(CFG['footprint_override']).to_crs(WORK_CRS)
        require_columns(override,['plotID','field_year','geometry'],'footprint_override')
    for _,p in meta.iterrows():
        try:
            loc=cached_location(p.plotID+'.basePlot.vst')
            cx,cy,crs=location_xy(loc)
            ids=evidence.loc[(evidence.plotID==p.plotID)&(evidence.field_year==p.field_year)&
                (evidence.stemDiameter>=CFG['min_field_dbh_cm']),'subplotID'].dropna().unique()
            polygons=[]
            for sub in ids:
                match=re.fullmatch(r'(\d+)_(100|400)',str(sub))
                if not match: continue
                d=cached_location(p.plotID+'.basePlot.vst.'+str(sub))
                x,y,c=location_xy(d); half=np.sqrt(int(match.group(2)))/2
                geo=gpd.GeoSeries([box(x-half,y-half,x+half,y+half)],crs=c).to_crs(WORK_CRS).iloc[0]
                polygons.append(geo)
            expected=float(p.totalSampledAreaTrees)
            if override is not None:
                q=override[(override.plotID==p.plotID)&(override.field_year==p.field_year)]
                if not q.empty: polygons=list(q.geometry)
            if not polygons: raise ValueError('100/400㎡ 소구획 경계 증거 없음')
            geom=unary_union(polygons)
            if not np.isfinite(expected) or expected<=0 or abs(geom.area-expected)/expected>0.05:
                raise ValueError(f'조사 면적 불일치: polygon={geom.area:.1f}, totalSampledAreaTrees={expected}; 빈 소구획/조사 범위 확인')
            x,y=Transformer.from_crs(crs,WORK_CRS,always_xy=True).transform(cx,cy)
            plots.append(dict(plotID=p.plotID,siteID=p.siteID,field_year=int(p.field_year),
                area_m2=expected,reported_area_m2=expected,polygon_area_m2=geom.area,
                area_difference_m2=geom.area-expected,center_x=x,center_y=y,geometry=geom))
        except requests.HTTPError as exc:
            if exc.response.status_code in [401,403,429]: raise
            audit('footprints',p.plotID,exc)
        except (ValueError,KeyError,TypeError) as exc: audit('footprints',p.plotID,exc)
    pd.DataFrame(AUDIT).to_csv(OUT/'audit.csv',index=False)
    if not plots: raise RuntimeError('검증된 조사 경계 없음: audit.csv와 subplotID/API 확인')
    return gpd.GeoDataFrame(plots,geometry='geometry',crs=WORK_CRS)

plots=build_footprints(plotmeta,evidence)
requested_cohort=set(CFG['comparison_plot_ids'])
missing_cohort=requested_cohort-set(plots.plotID)
if missing_cohort: raise RuntimeError('기존 비교 조사구 경계 복원 실패: '+str(sorted(missing_cohort)))
plots=plots[plots.plotID.isin(requested_cohort)].copy()
if plots.empty: raise RuntimeError('검증된 조사 경계가 없습니다. AUDIT와 subplotID/API 확인 필요')
eligible=[]
for _,p in plots.iterrows():
    f=field_geo[(field_geo.plotID==p.plotID)&(field_geo.field_year==p.field_year)].copy()
    inside=f.geometry.apply(p.geometry.covers)
    for ind in f.loc[~inside,'individualID']: audit('outside_sampled_area',ind,p.plotID)
    eligible.append(f[inside])
field_geo=pd.concat(eligible,ignore_index=True)
field_geo=gpd.GeoDataFrame(field_geo,geometry='geometry',crs=WORK_CRS)
plots.to_file(OUT/'sampled_footprints.gpkg',driver='GPKG')
pd.DataFrame(AUDIT).to_csv(OUT/'audit.csv',index=False)
display(plots.drop(columns='geometry'))
''')
section('06. NEONForestAGB 참조값과 수종 계수', '''152MB 파일은 청크 단위로 읽습니다. Chojnacky 유효값을 우선하고 없으면 Jenkins를 선택합니다.
중복 개체는 날짜·공식 우선순위를 기록해 한 줄로 만듭니다. 참조값은 각 조사구의 선택된 실측 연도와 맞춥니다.
2020/2022 보완값은 2021년 값으로 바꾸어 부르지 않습니다.

Chojnacky에는 Jenkins의 “mixed hardwood”와 같은 단일 명칭의 계수가 없습니다.
따라서 기본은 **일반 활엽수 대체식** `exp(-2.2118 + 2.4133 ln(DBH_cm))`이며
Cor/Eri/Lau/Platanaceae/Rosaceae/Ulmaceae 그룹임을 명시합니다. 모든 활엽수에 검증된 혼합식이라는 뜻이 아닙니다.
실측 taxonID에 대응하는 β0·β1은 MasterTaxonList에서 읽고, 미대응 수종은 대체 여부를 표시합니다.''', r'''
AGB_URL='https://ndownloader.figshare.com/files/45710721'
TAXON_URL='https://raw.githubusercontent.com/atkinsjeff/NEONForestAGB/main/data/MasterTaxonList.csv'
agbp=cached_download(AGB_URL,ROOT/'data/raw/NEONForestAGBv2.csv')
taxp=cached_download(TAXON_URL,ROOT/'data/raw/MasterTaxonList.csv')
tax=pd.read_csv(taxp,encoding='cp1252')
require_columns(tax,['PLANTS.Code','beta0','beta1','chojnacky'],'MasterTaxonList')
for c in ['beta0','beta1']: tax[c]=pd.to_numeric(tax[c],errors='coerce')
valid=tax.dropna(subset=['beta0','beta1'])
conflicts=valid.groupby('PLANTS.Code')[['beta0','beta1']].nunique().max(axis=1)
if conflicts.gt(1).any():
    bad_taxa=conflicts[conflicts.gt(1)].index.tolist()
    warnings.warn('Chojnacky 계수 충돌(QUDO 등): 경고 후 해당 taxon 계수 행 제외 '+str(bad_taxa))
    valid[valid['PLANTS.Code'].isin(bad_taxa)].to_csv(OUT/'excluded_taxon_coefficients.csv',index=False)
    valid=valid[~valid['PLANTS.Code'].isin(bad_taxa)].copy()
taxlookup=valid.drop_duplicates('PLANTS.Code').set_index('PLANTS.Code')
default=CFG['agb_default']
assert ((valid.beta0-default['beta0']).abs()<1e-8).any(), '기본 계수가 출처 표에 없음'
parts=[]
for chunk in pd.read_csv(agbp,chunksize=100000,low_memory=False):
    require_columns(chunk,['individualID','siteID','date','plantStatus','allometry','AGB'],'NEONForestAGB')
    chunk['field_year']=pd.to_datetime(chunk.date,errors='coerce').dt.year
    z=chunk[chunk.siteID.isin(CFG['sites']) & chunk.field_year.isin(plots.field_year.unique()) &
        chunk.plantStatus.fillna('').str.startswith('Live') &
        chunk.allometry.isin(['AGBChojnacky','AGBJenkins'])].copy()
    z['AGB']=pd.to_numeric(z.AGB,errors='coerce')
    parts.append(z[z.AGB.gt(0)])
ref=pd.concat(parts,ignore_index=True)
ref['_priority']=ref.allometry.map({'AGBChojnacky':0,'AGBJenkins':1})
ref=ref.sort_values(['individualID','field_year','_priority','date'],ascending=[True,True,True,False])
ref=ref.drop_duplicates(['individualID','field_year'])
ref=ref[['individualID','field_year','AGB','allometry','date']].rename(columns={'AGB':'ref_AGB_kg','date':'ref_date'})
field_geo=field_geo.merge(ref,on=['individualID','field_year'],how='left',validate='one_to_one')
field_all=field.merge(ref,on=['individualID','field_year'],how='left',validate='one_to_one')
field_all=field_all[field_all.plotID.isin(plots.plotID)].copy()
field_all['ref_fallback_jenkins']=field_all.allometry.eq('AGBJenkins')
field_all['has_matching_coordinate']=field_all.individualID.isin(field_geo.individualID)
field_geo['ref_fallback_jenkins']=field_geo.allometry.eq('AGBJenkins')
for ind in field_all.loc[field_all.ref_AGB_kg.isna(),'individualID']: audit('missing_AGB',ind,'참조 AGB 없음; 총량 비교에서 완전성 확인')
field_geo['species_beta0']=field_geo.taxonID.map(taxlookup.beta0)
field_geo['species_beta1']=field_geo.taxonID.map(taxlookup.beta1)
field_geo['species_group']=field_geo.taxonID.map(taxlookup.chojnacky)
field_geo['species_coeff_missing']=field_geo.species_beta0.isna()|field_geo.species_beta1.isna()
field_geo['species_beta0']=field_geo.species_beta0.fillna(default['beta0'])
field_geo['species_beta1']=field_geo.species_beta1.fillna(default['beta1'])
print('참조값 coverage:',field_geo.ref_AGB_kg.notna().mean())
print('Jenkins 대체:',field_geo.ref_fallback_jenkins.sum())
''')
section('07. 공식·단위 검산', '''H·CD는 m, DBH는 cm, AGB는 kg, 면적은 ㎡입니다.
Jucker/itcSegment biome=11: DBH=0.367(H×CD)^0.887×exp(0²/2).
로그는 자연로그입니다. 수관 지름은 사용자가 요청한 상자 가로·세로 평균을 사용합니다.
이는 Jucker 원자료의 실제 수관 지름과 기하학적으로 다를 수 있습니다.
출처 원문과 SHA256을 보존해 추후 계수 변경을 추적합니다.''', r'''
JUCKER_SOURCE='https://raw.githubusercontent.com/cran/itcSegment/master/R/dbh.R'
jsp=cached_download(JUCKER_SOURCE,ROOT/'data/raw/itcSegment_dbh.R')
source=jsp.read_text()
assert re.search(r'11\s*,\s*0\.367\s*,\s*0\.887\s*,\s*0',source), 'Jucker 출처 계수 변경: 검토 필요'

def estimate_dbh(H_m, CD_m, coefficients=None):
    c=coefficients or CFG['jucker']
    H=np.asarray(H_m,dtype=float); CD=np.asarray(CD_m,dtype=float)
    return np.where((H>0)&(CD>0),c['a']*(H*CD)**c['b']*np.exp(c.get('g',0)**2/2),np.nan)

def estimate_agb(dbh_cm,beta0=None,beta1=None):
    b0=default['beta0'] if beta0 is None else beta0
    b1=default['beta1'] if beta1 is None else beta1
    D=np.asarray(dbh_cm,dtype=float)
    with np.errstate(invalid='ignore',divide='ignore'):
        return np.where(D>0,np.exp(b0+b1*np.log(D)),np.nan)

def kg_to_t_ha(kg,area_m2): return np.asarray(kg)/1000*10000/np.asarray(area_m2)
def carbon_kg(agb_kg,fraction=0.5): return np.asarray(agb_kg)*fraction

assert np.isclose(kg_to_t_ha(1000,10000),1)
assert np.isclose(kg_to_t_ha(100,400),2.5)
assert np.isclose(carbon_kg(100,0.5),50)
assert np.isclose(carbon_kg(100,0.47),47)
assert np.allclose(offset_xy(100,200,10,90),(110,200))
assert np.allclose(offset_xy(100,200,10,0),(100,210))
assert np.isclose(estimate_dbh([20],[5])[0],0.367*100**0.887)
assert np.isclose(estimate_agb([10])[0],np.exp(-2.2118+2.4133*np.log(10)))
provenance={str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in [taxp,jsp]}
(OUT/'formula_provenance.json').write_text(json.dumps(dict(urls=[JUCKER_SOURCE,TAXON_URL,AGB_URL],sha256=provenance),indent=2))
print('단위·방위각·공식 검산 통과')
''')
section('08. AOP 제공 연월 조회·영상 연도 선택·타일 다운로드', '''공식 products API의 사이트별 availableMonths와 availableReleases를 출력합니다.
ORNL의 2021-06 RGB/CHM은 PROVISIONAL이므로 `include_provisional=True`가 필요합니다.
2021년 공통 자료가 실제로 없을 때만 ORNL에 한해 2022→2020 순서로 영상 연도를 선택합니다.
실측 연도는 바꾸지 않습니다. ORNL은 UTM 16N, 작업 좌표는 UTM 17N이므로 다운로드 좌표를 변환합니다.
다운로드 후 TIFF 파일 수·확장자·실제 저장 폴더를 출력하고, 다음 셀에서 조사구별 커버 여부를 확인합니다.''', r'''
def tile_coordinates(centers,buffer_m=50):
    tiles=set()
    for x,y in centers:
        for e in range(int((x-buffer_m)//1000),int((x+buffer_m)//1000)+1):
            for n in range(int((y-buffer_m)//1000),int((y+buffer_m)//1000)+1): tiles.add((e*1000,n*1000))
    return sorted(tiles)

def product_site_entry(product,site):
    records=product.get('data',product).get('siteCodes',[])
    return next((r for r in records if r.get('siteCode')==site),dict(siteCode=site,availableMonths=[],availableReleases=[]))

def select_aop_pair(rgb_entry,chm_entry,target_year,allow_adjacent=False,priority=(2022,2020),include_provisional=True):
    common=sorted(set(rgb_entry.get('availableMonths',[])) & set(chm_entry.get('availableMonths',[])))
    years=[int(target_year)]+([int(y) for y in priority] if allow_adjacent else [])
    selected=next((y for y in years if any(m.startswith(str(y)+'-') for m in common)),None)
    if selected is None:
        raise ValueError(f'RGB/CHM 공통 연월 없음: 요청 연도 {years}; 공통 월 {common}')
    months=[m for m in common if m.startswith(str(selected)+'-')]
    provisional=False
    for entry in [rgb_entry,chm_entry]:
        stable={m for r in entry.get('availableReleases',[]) if r.get('release')!='PROVISIONAL' for m in r.get('availableMonths',[])}
        prov={m for r in entry.get('availableReleases',[]) if r.get('release')=='PROVISIONAL' for m in r.get('availableMonths',[])}
        provisional |= any(m in prov and m not in stable for m in months)
    if provisional and not include_provisional:
        raise ValueError('선택 연도에 PROVISIONAL 자료가 있습니다. aop_include_provisional=True 필요; 다른 해로 숨겨 대체하지 않습니다.')
    return dict(imagery_year=selected,common_months=months,include_provisional=bool(provisional),
                provisional=bool(provisional),adjacent_imagery=selected!=int(target_year))

def fetch_product_metadata(dpid):
    path=ROOT/'data/cache'/f'product_{dpid}.json'
    if path.exists() and not CFG.get('api_refresh',True): return json.loads(path.read_text())
    url='https://data.neonscience.org/api/v0/products/'+dpid
    r=session.get(url,headers={'X-API-Token':TOKEN},timeout=90); r.raise_for_status()
    payload=r.json(); path.write_text(json.dumps(payload,ensure_ascii=False))
    return payload

def aop_tiff_files(folder):
    return sorted(p for p in Path(folder).rglob('*') if p.is_file() and p.suffix.lower() in ['.tif','.tiff'])

PRODUCTS=['DP3.30010.001','DP3.30015.001']
product_metadata={product:fetch_product_metadata(product) for product in PRODUCTS}
AOP_SELECTION={}; AOP_DIRS={}; availability_rows=[]; inventory_rows=[]
for site,g in plots.groupby('siteID'):
    entries={product:product_site_entry(product_metadata[product],site) for product in PRODUCTS}
    for product,entry in entries.items():
        print(site,product,'실제 제공 연월:',entry.get('availableMonths',[]))
        print('  release별:',entry.get('availableReleases',[]))
        for month in entry.get('availableMonths',[]):
            releases=[r['release'] for r in entry.get('availableReleases',[]) if month in r.get('availableMonths',[])]
            availability_rows.append(dict(siteID=site,product=product,month=month,releases=';'.join(releases)))
    selection=select_aop_pair(entries[PRODUCTS[0]],entries[PRODUCTS[1]],CFG['year'],
        CFG['imagery_allow_adjacent_years'].get(site,False),CFG['imagery_adjacent_priority'],CFG['aop_include_provisional'])
    AOP_SELECTION[site]=selection
    print(site,'영상 선택:',selection)
    if all(p in CFG['local_rgb'] and p in CFG['local_chm'] for p in g.plotID):
        print(site,'사용자 입력 RGB/CHM 사용: 다운로드 생략')
        continue
    loc=cached_location(g.iloc[0].plotID+'.basePlot.vst')
    _,_,site_crs=location_xy(loc)
    transform=Transformer.from_crs(WORK_CRS,site_crs,always_xy=True)
    centers=[transform.transform(r.center_x,r.center_y) for _,r in g.iterrows()]
    nt=len(tile_coordinates(centers,CFG['buffer_m']))
    print(site,'다운로드 CRS:',site_crs.to_string(),'첫 중심 좌표:',centers[0])
    print(site,'요청 예상 1km 타일:',nt,'RGB 압축 전 대략',round(nt*0.30,2),'GB')
    for product in PRODUCTS:
        legacy=ROOT/'data/raw/aop'/site
        # 기존 2021 캐시를 재사용. 인접연도는 별도 폴더에 격리.
        dest=legacy if selection['imagery_year']==CFG['year'] and (legacy/product).exists() else legacy/str(selection['imagery_year'])
        dest.mkdir(parents=True,exist_ok=True)
        args=dict(dpid=product,site=site,year=selection['imagery_year'],easting=[c[0] for c in centers],
            northing=[c[1] for c in centers],buffer=CFG['buffer_m'],check_size=False,
            savepath=str(dest),token=TOKEN,skip_if_exists=True,overwrite='yes',
            include_provisional=selection['include_provisional'],verbose=True)
        if 'skip_if_exists' not in inspect.signature(nu.by_tile_aop).parameters:
            raise RuntimeError('neonutilities를 업데이트하세요: 안전한 다운로드 재사용 API 없음')
        nu.by_tile_aop(**args)
        product_folder=dest/product
        # 확장자는 대소문자 구분 없이 .tif와 .tiff 모두 검색. 제품 디렉터리는 재귀적으로 검사.
        paths=aop_tiff_files(product_folder)
        if not paths:
            candidates=[q for q in dest.rglob('*') if q.is_dir() and q.name==product]
            if len(candidates)==1: product_folder=candidates[0]; paths=aop_tiff_files(product_folder)
        AOP_DIRS[(site,product)]=product_folder
        extensions=sorted({p.suffix for p in paths})
        print(site,product,'다운로드 후 TIFF 타일 수:',len(paths),'확장자:',extensions,'폴더:',product_folder)
        inventory_rows.append(dict(siteID=site,product=product,imagery_year=selection['imagery_year'],
            site_crs=site_crs.to_string(),expected_tiles=nt,tiff_tiles=len(paths),extensions=';'.join(extensions),folder=str(product_folder)))
        if not paths:
            print('폴더 진단:',[str(p.relative_to(dest)) for p in list(dest.rglob('*'))[:20]])
            raise RuntimeError(f'{site} {product}: 다운로드 후 TIFF 0개. 제공 월·PROVISIONAL·좌표·by_tile_aop 로그를 확인하세요.')
pd.DataFrame(availability_rows).to_csv(OUT/'aop_availability.csv',index=False)
pd.DataFrame(inventory_rows).to_csv(OUT/'aop_download_inventory.csv',index=False)
(OUT/'aop_selection.json').write_text(json.dumps(AOP_SELECTION,ensure_ascii=False,indent=2))
plots['imagery_year']=plots.siteID.map(lambda s:AOP_SELECTION[s]['imagery_year'])
plots['imagery_provisional']=plots.siteID.map(lambda s:AOP_SELECTION[s]['provisional'])
''')
section('09. 타일 색인·조사구 주변 잘라내기', '''1km RGB 전체를 한 번에 메모리에 읽지 않습니다.
조사구 주변 창만 mosaic하며, CRS·해상도·nodata를 검사합니다. 조사구 경계를 완전히 덮지 못하는 자료는
평가하지 않고 사유를 저장합니다. 여러 촬영일 자료가 겹치면 중복 합성을 피하기 위해 촬영 폴더를 먼저 확인하세요.
RGB·CHM 파일 경로와 CRS는 manifest에 남깁니다.''', r'''
def raster_index(folder,expected_resolution):
    rows=[]
    for path in sorted(Path(folder).rglob('*')):
        if path.suffix.lower() not in ['.tif','.tiff']: continue
        with rasterio.open(path) as ds:
            if not ds.crs or not ds.crs.is_projected: raise ValueError(f'미터 투영 CRS 필요: {path}')
            if not np.allclose(ds.res,expected_resolution,atol=expected_resolution*0.05):
                audit('resolution',path,f'예상 {expected_resolution}; 실제 {ds.res}'); continue
            rows.append(dict(path=str(path),crs=ds.crs,geometry=box(*ds.bounds)))
    return rows

def crop_mosaic(index,region,region_crs,dest,bands):
    if dest.exists(): return dest
    datasets=[]; coverage=[]
    try:
        for r in index:
            geom=gpd.GeoSeries([region],crs=region_crs).to_crs(r['crs']).iloc[0]
            if r['geometry'].intersects(geom):
                datasets.append(rasterio.open(r['path'])); coverage.append(r['geometry'])
        if not datasets: raise ValueError('해당 위치를 덮는 타일 없음')
        crs=datasets[0].crs
        if any(d.crs!=crs for d in datasets): raise ValueError('mosaic 입력 CRS 불일치')
        geom=gpd.GeoSeries([region],crs=region_crs).to_crs(crs).iloc[0]
        if not unary_union(coverage).buffer(0.2).covers(geom): raise ValueError('조사구 주변 영상 coverage 불완전')
        data,trans=merge(datasets,bounds=geom.bounds,indexes=bands)
        profile=datasets[0].profile.copy()
        profile.update(driver='GTiff',height=data.shape[1],width=data.shape[2],count=len(bands),
            transform=trans,compress='deflate',tiled=True,blockxsize=256,blockysize=256)
        with rasterio.open(dest,'w',**profile) as out: out.write(data)
    finally:
        for d in datasets: d.close()
    return dest

def raster_covers_geometry(index,geometry,geometry_crs):
    if not index: return False
    # 타일 CRS가 여러 개일 수 있으므로 모든 경계를 작업 CRS로 옮겨 합칩니다.
    bounds=[gpd.GeoSeries([r['geometry']],crs=r['crs']).to_crs(geometry_crs).iloc[0] for r in index]
    return unary_union(bounds).buffer(0.2).covers(geometry)

RASTERS={}; MANIFEST=[]; coverage_rows=[]
for site,g in plots.groupby('siteID'):
    rgb_idx=raster_index(AOP_DIRS.get((site,'DP3.30010.001'),ROOT/'data/raw/aop'/site/'DP3.30010.001'),0.1)
    chm_idx=raster_index(AOP_DIRS.get((site,'DP3.30015.001'),ROOT/'data/raw/aop'/site/'DP3.30015.001'),1.0)
    print(site,'해상도 검사 후 RGB/CHM 타일 수:',len(rgb_idx),len(chm_idx))
    for _,p in g.iterrows():
        key=p.plotID
        status=dict(siteID=site,plotID=key,imagery_year=p.imagery_year,field_year=p.field_year,
            rgb_sampled_coverage=False,chm_sampled_coverage=False,usable=False,reason='')
        try:
            if key in CFG['local_rgb'] and key in CFG['local_chm']:
                rgb=Path(CFG['local_rgb'][key]); chm=Path(CFG['local_chm'][key])
                with rasterio.open(rgb) as ds:
                    if ds.dtypes[0]!='uint8' or ds.count<3: raise ValueError('사용자 RGB는 uint8 3밴드 필요')
                RASTERS[key]=(rgb,chm)
                MANIFEST.append(dict(plotID=key,rgb=str(rgb),chm=str(chm),field_year=p.field_year,
                    imagery_year=p.imagery_year,imagery_provisional=p.imagery_provisional))
                status.update(usable=True,reason='사용자 입력 파일')
                coverage_rows.append(status)
                continue
            status['rgb_sampled_coverage']=raster_covers_geometry(rgb_idx,p.geometry,WORK_CRS)
            status['chm_sampled_coverage']=raster_covers_geometry(chm_idx,p.geometry,WORK_CRS)
            if not rgb_idx or not chm_idx: raise ValueError(f'{site}: {p.imagery_year} RGB/CHM TIFF 없음. AOP_DIRS 확인')
            region=box(p.center_x-CFG['buffer_m'],p.center_y-CFG['buffer_m'],
                       p.center_x+CFG['buffer_m'],p.center_y+CFG['buffer_m'])
            rgb=crop_mosaic(rgb_idx,region,WORK_CRS,ROOT/'data/crops'/f'{key}_{p.imagery_year}_rgb_{RUN_ID}.tif',[1,2,3])
            chm=crop_mosaic(chm_idx,region,WORK_CRS,ROOT/'data/crops'/f'{key}_{p.imagery_year}_chm_{RUN_ID}.tif',[1])
            with rasterio.open(rgb) as ds:
                if ds.dtypes[0]!='uint8': raise ValueError('RGB는 uint8 입력 필요; 스케일링을 명시하세요')
            RASTERS[key]=(rgb,chm)
            MANIFEST.append(dict(plotID=key,rgb=str(rgb),chm=str(chm),field_year=p.field_year,
                imagery_year=p.imagery_year,imagery_provisional=p.imagery_provisional))
            status['usable']=True
        except ValueError as exc:
            audit('raster_coverage',key,exc); status['reason']=str(exc)
        coverage_rows.append(status)
pd.DataFrame(MANIFEST).to_csv(OUT/'raster_manifest.csv',index=False)
coverage=pd.DataFrame(coverage_rows)
coverage.to_csv(OUT/'aop_plot_coverage.csv',index=False)
display(coverage)
print('사이트별 분석 가능 조사구:'); display(coverage.groupby('siteID').usable.agg(['sum','count']))
pd.DataFrame(AUDIT).to_csv(OUT/'audit.csv',index=False)
for site in CFG['sites']:
    if not coverage.loc[coverage.siteID==site,'usable'].any():
        raise RuntimeError(f'{site} 분석 가능 조사구 0개. aop_plot_coverage.csv 확인 후 다운로드/좌표 문제를 해결하세요.')
plots=plots[plots.plotID.isin(RASTERS)].copy()
field_geo=field_geo[field_geo.plotID.isin(RASTERS)].copy()
field_all=field_all[field_all.plotID.isin(RASTERS)].copy()
if plots.empty: raise RuntimeError('영상 coverage를 만족하는 조사구가 없습니다')
''')
section('10. DeepForest 사전학습 모델과 버전 호환', '''설치된 버전과 실제 메서드 존재 여부를 확인합니다.
새 버전은 `load_model`, 이전 버전은 `use_release`를 사용합니다.
네트워크·가중치 다운로드 오류를 구버전 오류로 오해해 숨기지 않습니다. 추가 학습은 하지 않습니다.''', r'''
from deepforest import main
def load_detector():
    version=importlib.metadata.version('deepforest')
    m=main.deepforest()
    if callable(getattr(m,'load_model',None)):
        m.load_model(model_name=CFG['model_name'],revision=CFG['model_revision'])
        method='load_model'
    elif callable(getattr(m,'use_release',None)):
        m.use_release(); method='use_release'
    else: raise RuntimeError(f'DeepForest {version}: 지원 모델 로더 없음')
    print('DeepForest',version,'로더',method,'predict_tile',inspect.signature(m.predict_tile))
    (OUT/'model.json').write_text(json.dumps(dict(version=version,loader=method,
        model=CFG['model_name'],revision=CFG['model_revision']),indent=2))
    return m
model=load_detector()
''')
section('11. 탐지·CHM 수고·수관 지름 추출', '''RGB 픽셀 좌표를 raster transform으로 지도 좌표로 바꿉니다.
수고는 상자 내부 CHM 최대값이고, 수관 지름은 상자 가로·세로 평균입니다.
CHM nodata/음수는 제외하며 유효 CHM이 없는 탐지는 따로 기록합니다.
조사구 안/밖 판정은 탐지 상자의 중심점 기준입니다. 실측은 줄기 위치 기준이므로 경계 편향이 남습니다.''', r'''
def zonal_max(path, geom, geom_crs):
    with rasterio.open(path) as ds:
        g=gpd.GeoSeries([geom],crs=geom_crs).to_crs(ds.crs).iloc[0]
        try: arr,_=mask(ds,[g],crop=True,filled=False,indexes=1)
        except ValueError: return np.nan
        v=arr.compressed(); v=v[np.isfinite(v)&(v>=0)]
        return float(v.max()) if len(v) else np.nan

def predict_crowns(m,rgb,chm,plot):
    cache=OUT/f'{plot.plotID}_raw_predictions.csv'
    if cache.exists(): raw=pd.read_csv(cache)
    else:
        params=inspect.signature(m.predict_tile).parameters
        path_key='path' if 'path' in params else 'raster_path'
        args={path_key:str(rgb),'patch_size':CFG['patch_size'],'patch_overlap':CFG['patch_overlap']}
        if 'iou_threshold' in params: args['iou_threshold']=CFG['nms_iou']
        raw=m.predict_tile(**args)
        if isinstance(raw,tuple): raw=raw[0]
        if raw is None or len(raw)==0:
            raw=pd.DataFrame(columns=['xmin','ymin','xmax','ymax','score'])
        raw.to_csv(cache,index=False)
    require_columns(raw,['xmin','ymin','xmax','ymax','score'],'DeepForest predictions')
    rows=[]
    with rasterio.open(rgb) as ds:
        trans=ds.transform; crs=ds.crs
        for idx,r in raw[raw.score>=CFG['score_threshold']].iterrows():
            corners=[trans*(x,y) for x,y in [(r.xmin,r.ymin),(r.xmax,r.ymin),(r.xmax,r.ymax),(r.xmin,r.ymax)]]
            geom=Polygon(corners)
            geom=gpd.GeoSeries([geom],crs=crs).to_crs(WORK_CRS).iloc[0]
            if not plot.geometry.covers(geom.centroid): continue
            h=zonal_max(chm,geom,WORK_CRS)
            if not np.isfinite(h): audit('CHM_nodata',f'{plot.plotID}_{idx}','CHM 없음'); continue
            if h<CFG['min_height_m']: continue
            width=(r.xmax-r.xmin)*abs(trans.a)
            length=(r.ymax-r.ymin)*abs(trans.e)
            rows.append(dict(det_id=f'{plot.plotID}_{idx}',plotID=plot.plotID,siteID=plot.siteID,
                field_year=plot.field_year,score=r.score,xmin_px=r.xmin,ymin_px=r.ymin,
                xmax_px=r.xmax,ymax_px=r.ymax,H_m=h,CD_m=(width+length)/2,
                width_m=width,length_m=length,geometry=geom))
    return rows

records=[]
for _,p in plots.iterrows():
    print('탐지:',p.plotID)
    records.extend(predict_crowns(model,*RASTERS[p.plotID],p))
DET_COLUMNS=['det_id','plotID','siteID','field_year','score','xmin_px','ymin_px','xmax_px','ymax_px','H_m','CD_m','width_m','length_m','geometry']
detections=gpd.GeoDataFrame(records,columns=DET_COLUMNS,geometry='geometry',crs=WORK_CRS)
detections['imagery_year']=detections.siteID.map(lambda s:AOP_SELECTION[s]['imagery_year'])
detections['DBH_jucker_cm']=estimate_dbh(detections.H_m,detections.CD_m)
detections['AGB_default_kg']=estimate_agb(detections.DBH_jucker_cm)
detections['C_default_kg']=carbon_kg(detections.AGB_default_kg,CFG['carbon_fraction'])
print('평가 영역 내 탐지:',len(detections))
''')
section('12. 보이는 나무 정의와 1:1 매칭', '''실측 수관위치가 Full sun / Open grown / Partially shaded 이거나,
실측 수고가 주변 반경 3m CHM 최대값의 80% 이상이면 “보이는 나무”입니다.
이 정의는 관측 가능성의 근사이며 높이 검증 표본 선택 편향을 일으킬 수 있어 수관위치만 쓴 플래그도 보존합니다.

헝가리안 매칭의 후보는 **줄기 점이 상자 안에 있는 조합**입니다.
먼저 가능한 매칭 수를 최대화하고, 그 안에서 큰 실측 DBH를 우선하며 중심거리는 동률 해소에만 씁니다.
보이는 나무/전체 나무 두 모집단을 각각 매칭합니다. 전체 나무 기준은 선택된 공통 비교 범위(DBH≥10cm)의 전체입니다.
점 정답은 실제 수관 IoU 평가와 같지 않습니다.''', r'''
def label_visibility(f,rasters):
    f=f.copy(); local=[]
    for _,r in f.iterrows():
        local.append(zonal_max(rasters[r.plotID][1],r.geometry.buffer(CFG['visible_radius_m']),WORK_CRS))
    f['local_CHM_m']=local
    canopy=f.get('canopyPosition',pd.Series('',index=f.index)).fillna('').str.strip().str.lower()
    f['visible_by_position']=canopy.isin(['full sun','open grown','partially shaded'])
    f['visible_by_height']=f.field_H_m.ge(CFG['visible_height_ratio']*f.local_CHM_m)&f.local_CHM_m.gt(0)
    f['visible']=f.visible_by_position|f.visible_by_height
    return f

def match_one_to_one(det,truth):
    if det.empty or truth.empty: return pd.DataFrame(columns=['det_id','individualID','plotID','distance_m'])
    d=det.reset_index(drop=True); t=truth.reset_index(drop=True)
    n,m=len(d),len(t)
    weights=np.zeros((n,m+n)) # dummy columns = unmatched
    valid=np.zeros((n,m),dtype=bool)
    sizes=pd.to_numeric(t[CFG['match_size']],errors='coerce').fillna(0).to_numpy()
    rank=pd.Series(sizes).rank(method='dense').to_numpy()/(m+1)
    for i,r in d.iterrows():
        for j,s in t.iterrows():
            if r.geometry.covers(s.geometry):
                valid[i,j]=True
                distance=r.geometry.centroid.distance(s.geometry)
                weights[i,j]=(n+1)*10 + rank[j] + 1e-6/(1+distance)
    ri,ci=linear_sum_assignment(-weights)
    rows=[]
    for i,j in zip(ri,ci):
        if j<m and valid[i,j]:
            rows.append(dict(det_id=d.iloc[i].det_id,individualID=t.iloc[j].individualID,
                plotID=d.iloc[i].plotID,distance_m=d.iloc[i].geometry.centroid.distance(t.iloc[j].geometry)))
    return pd.DataFrame(rows,columns=['det_id','individualID','plotID','distance_m'])

field_geo=label_visibility(field_geo,RASTERS)
MATCHES={}; detmetrics=[]
for population in ['visible','all']:
    items=[]
    for _,p in plots.iterrows():
        d=detections[detections.plotID==p.plotID]
        t=field_geo[field_geo.plotID==p.plotID]
        if population=='visible': t=t[t.visible]
        pairs=match_one_to_one(d,t); items.append(pairs)
        tp=len(pairs); precision=tp/len(d) if len(d) else np.nan
        recall=tp/len(t) if len(t) else np.nan
        f1=2*tp/(len(d)+len(t)) if len(d)+len(t) else np.nan
        detmetrics.append(dict(plotID=p.plotID,siteID=p.siteID,field_year=p.field_year,
            population=population,TP=tp,FP=len(d)-tp,FN=len(t)-tp,
            n_detections=len(d),n_truth=len(t),precision=precision,recall=recall,F1=f1))
    MATCHES[population]=pd.concat(items,ignore_index=True)
detmetrics=pd.DataFrame(detmetrics)
pairs=MATCHES['visible']
matched=pairs.merge(detections.drop(columns='geometry'),on=['det_id','plotID'],validate='one_to_one')
matched=matched.merge(field_geo.drop(columns=['geometry','plotID','siteID','field_year']),on='individualID',validate='one_to_one')
display(detmetrics)
''')
section('13. 조사구 분리 DBH 회귀와 수종 민감도', '''log(DBH) ~ log(H) + log(CD)를 적합합니다. GroupKFold의 그룹은 조사구입니다.
동일 조사구 나무가 학습과 평가에 동시에 들어가지 않습니다. 로그 역변환의 Duan smearing도 학습 fold에서만 계산합니다.
공식 비교용 회귀값은 out-of-fold 값만 사용합니다. 전체 표본 적합 모델은 재사용 목적으로 별도 저장하며 평가값으로 쓰지 않습니다.
유효 조사구가 2개 미만이면 회귀 지표는 계산하지 않습니다.''', r'''
def grouped_dbh_regression(df,n_splits=5):
    result=pd.Series(np.nan,index=df.index); fold_id=pd.Series(np.nan,index=df.index)
    good=df[['H_m','CD_m','field_DBH_cm']].gt(0).all(axis=1)
    good &= np.isfinite(df[['H_m','CD_m','field_DBH_cm']]).all(axis=1)
    z=df.loc[good]; groups=z.plotID
    if groups.nunique()<2:
        audit('regression','GroupKFold','독립 조사구 2개 미만'); return result,fold_id,None
    X=np.log(z[['H_m','CD_m']].to_numpy()); y=np.log(z.field_DBH_cm.to_numpy())
    folds=GroupKFold(n_splits=min(n_splits,groups.nunique()))
    for fold,(tr,te) in enumerate(folds.split(X,y,groups)):
        if len(tr)<4:
            audit('regression',fold,'학습 나무 4개 미만'); continue
        m=LinearRegression().fit(X[tr],y[tr])
        smear=float(np.mean(np.exp(y[tr]-m.predict(X[tr]))))
        result.loc[z.index[te]]=np.exp(m.predict(X[te]))*smear
        fold_id.loc[z.index[te]]=fold
    m=LinearRegression().fit(X,y)
    fit=dict(intercept=float(m.intercept_),log_H=float(m.coef_[0]),log_CD=float(m.coef_[1]),
             smearing=float(np.mean(np.exp(y-m.predict(X)))),n=len(z),n_plots=groups.nunique())
    return result,fold_id,fit

if not matched.empty:
    matched['DBH_regression_oof_cm'],matched['cv_fold'],final_fit=grouped_dbh_regression(matched,CFG['cv_folds'])
    matched['AGB_regression_oof_kg']=estimate_agb(matched.DBH_regression_oof_cm)
    matched['AGB_known_species_kg']=estimate_agb(matched.DBH_jucker_cm,matched.species_beta0,matched.species_beta1)
    matched['species_effect_kg']=matched.AGB_known_species_kg-matched.AGB_default_kg
    matched['species_effect_pct']=100*(matched.AGB_known_species_kg/matched.AGB_default_kg-1)
    if final_fit: (OUT/'dbh_regression_full_fit.json').write_text(json.dumps(final_fit,indent=2))
else:
    for c in ['DBH_regression_oof_cm','AGB_regression_oof_kg','AGB_known_species_kg','species_effect_kg','species_effect_pct']:
        matched[c]=pd.Series(dtype=float)
    audit('matching','visible','매칭 0건: 회귀·개체 오차표 계산 불가')
''')
section('14. 단계별 오차와 조사구 총량·누락 분해', '''Bias=추정−참조입니다. 상관계수는 Pearson r이며 표본 수와 유효 조사구 수를 함께 봅니다.
수관 지름은 실측 두 축 평균과 비교합니다. 수관 측정이 없는 행은 그 단계에서만 제외합니다.

조사구 참조 총량은 좌표 필터 전 전체 살아 있는 DBH≥10cm 나무를 합칩니다.
좌표 누락은 총량을 NaN으로 만들지 않습니다. 참조 AGB 자체가 없는 나무가 있을 때만 NaN 처리합니다.
분모는 totalSampledAreaTrees를 우선하며 폴리곤 면적도 함께 기록합니다.
차이는 `매칭 나무 추정오차 + 미매칭 탐지량 − 매칭되지 않은 전체 실측량`으로 분해합니다.
매칭 불가 좌표 누락과 실제 탐지 누락을 혼동하지 않도록 별도 열로 나눕니다.
영상 연도와 실측 연도가 다르면 사이트·연도 조합별로 분리합니다.
수종 민감도는 같은 매칭 나무 표본끼리만 비교합니다. 미매칭 탐지의 수종을 아는 척하지 않습니다.''', r'''
def errors(truth,pred):
    a=np.asarray(truth,dtype=float); b=np.asarray(pred,dtype=float)
    ok=np.isfinite(a)&np.isfinite(b); a=a[ok]; b=b[ok]
    if not len(a): return dict(n=0,MAE=np.nan,RMSE=np.nan,Bias=np.nan,r=np.nan)
    e=b-a
    r=np.corrcoef(a,b)[0,1] if len(a)>1 and np.std(a)>0 and np.std(b)>0 else np.nan
    return dict(n=len(a),MAE=float(np.mean(abs(e))),RMSE=float(np.sqrt(np.mean(e**2))),Bias=float(np.mean(e)),r=r)

def plot_biomass_totals(detections_all,truth_all,truth_coord,pairs,plot,carbon_fraction=0.5,carbon_alternative=0.47):
    # 일반 list[dict] 입력으로 동작: 좌표/지리 패키지 없이도 총량 계산을 재사용 가능.
    def number(value):
        try: return float(value)
        except (ValueError,TypeError): return float('nan')
    def sum_valid(rows,column):
        values=[number(r.get(column)) for r in rows]
        return sum(values) if all(math.isfinite(v) for v in values) else float('nan')
    area=number(plot.get('reported_area_m2'))
    if not math.isfinite(area) or area<=0: area=number(plot.get('area_m2'))
    if not math.isfinite(area) or area<=0: raise ValueError('양수 조사면적 필요')
    scale=lambda value:value/1000*10000/area
    by_id={r['individualID']:r for r in truth_all}
    if len(by_id)!=len(truth_all): raise ValueError('조사구 참조 개체 중복')
    det_by_id={r['det_id']:r for r in detections_all}
    paired_truth={r['individualID'] for r in pairs}
    paired_det={r['det_id'] for r in pairs}
    coord_ids={r['individualID'] for r in truth_coord}
    if not paired_truth<=by_id.keys() or not paired_det<=det_by_id.keys(): raise ValueError('매칭 ID가 조사구 모집단에 없음')
    reference=sum_valid(truth_all,'ref_AGB_kg')
    pred=sum_valid(detections_all,'AGB_default_kg')
    coord_reference=sum_valid(truth_coord,'ref_AGB_kg')
    missing_agb=sum(not math.isfinite(number(r.get('ref_AGB_kg'))) for r in truth_all)
    matched_ref=sum_valid([by_id[k] for k in paired_truth],'ref_AGB_kg')
    matched_pred=sum_valid([det_by_id[k] for k in paired_det],'AGB_default_kg')
    unmatched_pred=sum_valid([r for r in detections_all if r['det_id'] not in paired_det],'AGB_default_kg')
    missed_ref=sum_valid([r for r in truth_all if r['individualID'] not in paired_truth],'ref_AGB_kg')
    no_coord=sum_valid([r for r in truth_all if r['individualID'] not in coord_ids],'ref_AGB_kg')
    missed_coord=sum_valid([r for r in truth_coord if r['individualID'] not in paired_truth],'ref_AGB_kg')
    field_year=int(plot['field_year']); imagery_year=int(plot.get('imagery_year',field_year))
    row=dict(plotID=plot['plotID'],siteID=plot['siteID'],field_year=field_year,imagery_year=imagery_year,
        same_year=field_year==imagery_year,imagery_provisional=plot.get('imagery_provisional',False),
        area_m2=area,reported_area_m2=number(plot.get('reported_area_m2',area)),
        polygon_area_m2=number(plot.get('polygon_area_m2',area)),
        area_difference_m2=number(plot.get('polygon_area_m2',area))-area,
        n_truth=len(truth_all),n_truth_coord=len(truth_coord),n_detected=len(detections_all),n_matched=len(pairs),
        n_field_before_coordinate_filter=len(truth_all),n_without_matching_coordinate=len(truth_all)-len(truth_coord),
        n_missing_ref_AGB=missing_agb,AGB_reference_coverage=1-missing_agb/len(truth_all) if truth_all else 1.0,
        reference_complete=missing_agb==0,coordinates_complete=len(truth_all)==len(truth_coord),
        AGB_est_t_ha=scale(pred),AGB_ref_t_ha=scale(reference),AGB_ref_coordonly_t_ha=scale(coord_reference),
        matched_error_t_ha=scale(matched_pred-matched_ref),unmatched_detection_t_ha=scale(unmatched_pred),
        missed_truth_t_ha=scale(missed_ref),unmapped_truth_t_ha=scale(no_coord),missed_mapped_truth_t_ha=scale(missed_coord))
    row['C_est_tC_ha']=row['AGB_est_t_ha']*carbon_fraction
    row['C_ref_tC_ha']=row['AGB_ref_t_ha']*carbon_fraction
    row['C_ref_coordonly_tC_ha']=row['AGB_ref_coordonly_t_ha']*carbon_fraction
    row['C_est_047_tC_ha']=row['AGB_est_t_ha']*carbon_alternative
    row['C_ref_047_tC_ha']=row['AGB_ref_t_ha']*carbon_alternative
    row['C_bias_tC_ha']=row['C_est_tC_ha']-row['C_ref_tC_ha']
    row['C_missed_tC_ha']=row['missed_truth_t_ha']*carbon_fraction
    row['estimate_reference_ratio']=pred/reference if reference>0 else float('nan')
    if math.isfinite(reference) and math.isfinite(pred):
        assert math.isclose(pred-reference,matched_pred-matched_ref+unmatched_pred-missed_ref,rel_tol=1e-9,abs_tol=1e-6)
    return row

def comparison_statistics(truth,pred):
    result=errors(truth,pred)
    a=np.asarray(truth,dtype=float); b=np.asarray(pred,dtype=float)
    ok=np.isfinite(a)&np.isfinite(b); a=a[ok]; b=b[ok]
    mean_ref=float(a.mean()) if len(a) else np.nan
    mean_est=float(b.mean()) if len(b) else np.nan
    result.update(mean_reference=mean_ref,mean_estimate=mean_est,
        estimate_reference_ratio=mean_est/mean_ref if mean_ref>0 else np.nan)
    return result

stage_rows=[]
comparisons=[('height','m','field_H_m','H_m'),('crown_diameter','m','field_CD_m','CD_m'),
 ('DBH_Jucker','cm','field_DBH_cm','DBH_jucker_cm'),('DBH_regression_OOF','cm','field_DBH_cm','DBH_regression_oof_cm'),
 ('AGB_default','kg','ref_AGB_kg','AGB_default_kg'),('AGB_known_species','kg','ref_AGB_kg','AGB_known_species_kg'),
 ('AGB_regression_OOF','kg','ref_AGB_kg','AGB_regression_oof_kg')]
for year in sorted(plots.field_year.unique()):
    z=matched[matched.field_year==year]
    for name,unit,refcol,predcol in comparisons:
        stage_rows.append(dict(stage=name,unit=unit,field_year=year,population='visible_matched',
            n_plots=z.plotID.nunique(),**errors(z[refcol],z[predcol])))
stage_errors=pd.DataFrame(stage_rows)

plot_rows=[]
for _,p in plots.iterrows():
    d=detections[detections.plotID==p.plotID]
    fa=field_all[(field_all.plotID==p.plotID)&(field_all.field_year==p.field_year)]
    fc=field_geo[(field_geo.plotID==p.plotID)&(field_geo.field_year==p.field_year)]
    ma=MATCHES['all']; ma=ma[ma.plotID==p.plotID]
    plot_rows.append(plot_biomass_totals(d.to_dict('records'),fa.to_dict('records'),fc.to_dict('records'),
        ma.to_dict('records'),p.to_dict(),CFG['carbon_fraction'],CFG['carbon_alternative']))
plot_comparison=pd.DataFrame(plot_rows)
for (year,image_year),z in plot_comparison.groupby(['field_year','imagery_year']):
    for name,unit,refcol,predcol in [('plot_AGB','t/ha','AGB_ref_t_ha','AGB_est_t_ha'),
      ('plot_carbon','tC/ha','C_ref_tC_ha','C_est_tC_ha')]:
        stage_rows.append(dict(stage=name,unit=unit,siteID='ALL',field_year=year,imagery_year=image_year,
            population='plot',n_plots=len(z),**comparison_statistics(z[refcol],z[predcol])))
        for site,q in z.groupby('siteID'):
            stage_rows.append(dict(stage=name,unit=unit,siteID=site,field_year=year,imagery_year=image_year,
                population='plot',n_plots=len(q),**comparison_statistics(q[refcol],q[predcol])))
stage_errors=pd.DataFrame(stage_rows)
for (year,pop),z in detmetrics.groupby(['field_year','population']):
    tp,fp,fn=z[['TP','FP','FN']].sum()
    stage_rows.append(dict(stage='detection',unit='ratio',field_year=year,population=pop,
        n=int(tp+fn),n_plots=z.plotID.nunique(),
        precision=tp/(tp+fp) if tp+fp else np.nan,recall=tp/(tp+fn) if tp+fn else np.nan,
        F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan))
stage_errors=pd.DataFrame(stage_rows)
display(stage_errors); display(plot_comparison)
''')
improvement_source=(HERE/'detection_improvements.py').read_text(encoding='utf-8')
improvement_helpers,improvement_run=improvement_source.split('# Same plot cohort for all methods;',1)
section('14A. 나무 찾기 개선 함수: DF 탐색·CHM 꼭대기·watershed·결합', '''이전 실행의 유지사항(QUDO 충돌 경고/제외, gzip 검사 예외, /content/hf_cache, ORNL provisional,
좌표 없는 나무 포함 참조 총량)을 유지합니다. 새로운 함수들은 이 셀에서 정의만 합니다.
A는 고정 사전학습 모델의 score 5개×patch 3개 조합을 탐색합니다. 학습 가중치는 바꾸지 않습니다.
B는 nodata 정규화 Gaussian(sigma=1m) → 높이 가변 원형 최대점 → 시드 watershed입니다.
창 지름=2.51+0.00901H²(m): Popescu & Wynne (2004) eq.3의 2.51503을 요청대로 반올림합니다.
B의 H는 평활 후 탐지한 픽셀의 원본 CHM 높이, CD는 실제 watershed 면적의 등가원 지름입니다.
매칭은 모든 방법에서 기존의 줄기 점이 상자 안에 있는 헝가리안 규칙을 유지하므로 B의 매칭 상자는 수관 polygon의 envelope입니다.
C는 A의 상자 안에서 최고 CHM 꼭대기로 H를 보완(CD는 상자 기준), 밖의 꼭대기는 B 수관으로 추가하고 1m 중복을 제거합니다.
상자에 꼭대기가 없으면 원래 DF 탐지를 보존합니다. 이 선택은 보고서에도 명시합니다.''',improvement_helpers)
section('14B. 동일 23개 조사구에서 네 방법 비교·기본값 선택', '''A: 바깥 GroupKFold 평가 조사구의 정답은 조합 선택에 쓰지 않습니다.
각 바깥 학습 조사구 안의 안쪽 GroupKFold **보이는 나무 기준 F1**로 조합을 고르고, 바깥 평가 탐지만 합쳐 OOF 결과를 만듭니다.
전 표본으로 고른 단일 조합은 향후 새 자료 적용용으로 별도 저장하며 현재 성능표에는 쓰지 않습니다.
PR 그림의 곡선은 조합 선택 점수, 별표는 바깥 OOF 성능입니다. B는 고정 설정, C는 A의 바깥 OOF 상자를 사용합니다.
4개 방법×전체/사이트 비교표를 저장하고, 23개 공통 조사구 탄소 RMSE가 최소인 방법을 이후 민감도·보고서의 기본값으로 씁니다.
이는 같은 23개에서 선택한 기본값이며 새 지역의 독립 검증 결과라는 뜻은 아닙니다.
추가 GPU 추론은 patch별 1회씩(총 69 조사구 추론)이며 threshold는 동일 예측을 필터링해 재사용합니다.
정확히 23개가 아니거나 요청한 4개 그림 조사구가 없으면 멈추므로 비교 모집단이 조용히 바뀌지 않습니다.''',
    '# Same plot cohort for all methods;'+improvement_run)

section('15. 수고 오차 ±10%·±20%의 탄소 영향', '''탐지·매칭·수관 지름을 고정하고 H만 바꾸는 통제된 민감도입니다.
수고 오차가 탐지 필터나 매칭까지 바꾸는 효과는 포함하지 않습니다.
Jucker DBH와 기본 AGB 식을 연쇄 적용하므로 탄소 배율은 이론적으로 `(1+오차)^(0.887×2.4133)`입니다.
이는 한국의 낮은 정확도 높이 자료를 검토하는 시나리오이며 실제 한국 자료의 정확도 검증을 대신하지 않습니다.''', r'''
sensitivity=[]
for _,p in plots.iterrows():
    z=detections[detections.plotID==p.plotID]
    baseline=float(kg_to_t_ha(z.AGB_default_kg.sum(),p.area_m2)*CFG['carbon_fraction'])
    for delta in [-0.2,-0.1,0.0,0.1,0.2]:
        D=estimate_dbh(z.H_m*(1+delta),z.CD_m)
        C=float(kg_to_t_ha(estimate_agb(D).sum(),p.area_m2)*CFG['carbon_fraction'])
        factor=(1+delta)**(CFG['jucker']['b']*default['beta1'])
        if baseline>0: assert np.isclose(C/baseline,factor)
        sensitivity.append(dict(plotID=p.plotID,field_year=p.field_year,height_error_pct=100*delta,
            C_tC_ha=C,baseline_tC_ha=baseline,change_tC_ha=C-baseline,
            change_pct=100*(C/baseline-1) if baseline>0 else np.nan,theoretical_factor=factor))
sensitivity=pd.DataFrame(sensitivity)
display(sensitivity)
''')
section('16. 조사구 지도와 단계별 산점도', '''RGB 위의 상자는 조사 범위 내 탐지, 점은 실측 줄기 위치입니다.
녹색은 보이는 나무, 자홍색은 가려진 나무이며 조사 경계를 노란색으로 표시합니다.
그림의 라벨은 Colab 기본 폰트에서도 깨지지 않도록 영어로 씁니다. 보고서는 한국어입니다.''', r'''
FIGDIR=OUT/'figures'; FIGDIR.mkdir(exist_ok=True)
for _,p in plots.iterrows():
    fig,ax=plt.subplots(figsize=(9,9))
    rgb,_=RASTERS[p.plotID]
    with rasterio.open(rgb) as ds:
        show(ds,ax=ax)
        crs=ds.crs
    gpd.GeoSeries([p.geometry],crs=WORK_CRS).to_crs(crs).boundary.plot(ax=ax,color='yellow',linewidth=2)
    d=detections[detections.plotID==p.plotID].to_crs(crs)
    if not d.empty: d.boundary.plot(ax=ax,color='cyan',linewidth=0.7)
    f=field_geo[field_geo.plotID==p.plotID].to_crs(crs)
    for v,color,label in [(True,'lime','Visible field stem'),(False,'magenta','Occluded field stem')]:
        q=f[f.visible==v]
        if not q.empty: ax.scatter(q.geometry.x,q.geometry.y,c=color,s=18,label=label,edgecolors='black',linewidths=0.3)
    ax.set_title(f'{p.plotID}: imagery {p.imagery_year} / field {p.field_year}')
    ax.legend(); fig.savefig(FIGDIR/f'{p.plotID}_overlay.png',dpi=160,bbox_inches='tight'); plt.close(fig)

def scatter_identity(df,x,y,label,filename):
    if x not in df or y not in df: return
    z=df[[x,y]].replace([np.inf,-np.inf],np.nan).dropna()
    fig,ax=plt.subplots(figsize=(5,5))
    if len(z):
        ax.scatter(z[x],z[y],s=12,alpha=0.6)
        lo=min(z.min()); hi=max(z.max())
        ax.plot([lo,hi],[lo,hi],'k--',label='1:1')
    ax.set(xlabel='Reference '+label,ylabel='Estimated '+label,title=f'{label} (n={len(z)})')
    fig.savefig(FIGDIR/filename,dpi=170,bbox_inches='tight'); plt.close(fig)

# 같은 해 검증 산점도와 인접연도 보완 그림을 섞지 않음
for year in sorted(plots.field_year.unique()):
    z=matched[matched.field_year==year]
    for x,y,label,fn in [('field_H_m','H_m','Height (m)','height'),
      ('field_CD_m','CD_m','Crown diameter (m)','crown'),
      ('field_DBH_cm','DBH_jucker_cm','DBH Jucker (cm)','dbh_jucker'),
      ('field_DBH_cm','DBH_regression_oof_cm','DBH OOF (cm)','dbh_oof')]:
        scatter_identity(z,x,y,label,f'{fn}_{year}.png')

def plot_carbon_scatter(plot_table,output_folder):
    for (year,image_year),cohort in plot_table.groupby(['field_year','imagery_year']):
        z=cohort.replace([np.inf,-np.inf],np.nan).dropna(subset=['C_ref_tC_ha','C_est_tC_ha'])
        fig,ax=plt.subplots(figsize=(6,6))
        for site,q in z.groupby('siteID'):
            ax.scatter(q.C_ref_tC_ha,q.C_est_tC_ha,s=35,alpha=0.8,label=site,
                color={'GRSM':'tab:blue','ORNL':'tab:orange'}.get(site))
        if len(z):
            hi=max(z[['C_ref_tC_ha','C_est_tC_ha']].max())*1.05
            hi=max(hi,1.0); ax.plot([0,hi],[0,hi],'k--',label='1:1')
            ax.set_xlim(0,hi); ax.set_ylim(0,hi)
        ax.set(xlabel='Field-based reference carbon (tC/ha)',ylabel='Estimated carbon (tC/ha)',
            title=f'Plot carbon: field {year} / imagery {image_year} (n={len(z)})')
        ax.legend(); ax.set_aspect('equal',adjustable='box')
        filename=f'plot_carbon_{year}.png' if year==image_year else f'plot_carbon_field{year}_imagery{image_year}.png'
        fig.savefig(Path(output_folder)/filename,dpi=170,bbox_inches='tight'); plt.close(fig)
plot_carbon_scatter(plot_comparison,FIGDIR)
''')
section('17. CSV·한국어 보고서·재사용 함수 저장', '''CSV의 기하정보는 WKT로 저장하며 좌표계는 별도 기록합니다.
현재 실행 결과는 `outputs/실행ID/`, 최근 결과 안내는 `outputs/report.md`에 기록합니다.
노트북에 정의한 함수는 독립 Python 모듈로도 내보내며, 전역 CFG·ROOT 등을 쓰는 NEON 전용 로더와
입력 배열·계수만 받는 일반 추정 함수를 구분합니다.''', r'''
def write_geo_csv(df,path):
    df=df.copy()
    if 'geometry' in df: df['geometry']=df.geometry.apply(lambda g:g.wkt if g is not None else None)
    pd.DataFrame(df).to_csv(path,index=False)

def build_site_summary(plot_table,detection_metrics,matched_table,requested_sites):
    rows=[]
    for site in requested_sites:
        site_plots=plot_table[plot_table.siteID==site]
        if site_plots.empty:
            rows.append(dict(siteID=site,n_plots=0,comparison_note='분석 결과 없음')); continue
        for (field_year,image_year),q in site_plots.groupby(['field_year','imagery_year']):
            ids=set(q.plotID)
            dm=detection_metrics[detection_metrics.plotID.isin(ids)]
            mt=matched_table[matched_table.plotID.isin(ids)]
            row=dict(siteID=site,field_year=field_year,imagery_year=image_year,
                comparison_note='영상 연도 = 실측 연도' if field_year==image_year else f'{site}은 영상 연도 ≠ 실측 연도',
                n_plots=len(q),n_field_trees=int(q.n_truth.sum()),n_coord_trees=int(q.n_truth_coord.sum()),
                n_detections=int(q.n_detected.sum()))
            for population,label in [('visible','F1_visible'),('all','F1_all')]:
                z=dm[dm.population==population]; tp=z.TP.sum(); fp=z.FP.sum(); fn=z.FN.sum()
                row[label]=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan
            for label,refcol,predcol in [('height_Bias_m','field_H_m','H_m'),
                ('DBH_Bias_cm','field_DBH_cm','DBH_jucker_cm'),('individual_AGB_Bias_kg','ref_AGB_kg','AGB_default_kg')]:
                row[label]=errors(mt[refcol],mt[predcol])['Bias']
            stats=comparison_statistics(q.C_ref_tC_ha,q.C_est_tC_ha)
            row.update(n_carbon_plots=stats['n'],C_est_mean_tC_ha=stats['mean_estimate'],
                C_ref_mean_tC_ha=stats['mean_reference'],estimate_reference_ratio=stats['estimate_reference_ratio'])
            rows.append(row)
    return pd.DataFrame(rows)

SUMMARY_LABELS={'siteID':'사이트','field_year':'실측 연도','imagery_year':'영상 연도','comparison_note':'연도 구분',
    'n_plots':'조사구 수','n_field_trees':'실측 나무 수(전체)','n_coord_trees':'매칭 좌표 나무 수',
    'n_detections':'탐지 수','F1_visible':'F1(보이는 나무)','F1_all':'F1(전체 매핑 나무)',
    'height_Bias_m':'수고 Bias(m)','DBH_Bias_cm':'DBH Bias(cm)','individual_AGB_Bias_kg':'개체 AGB Bias(kg)',
    'n_carbon_plots':'탄소 비교 조사구 수','C_est_mean_tC_ha':'평균 추정(tC/ha)',
    'C_ref_mean_tC_ha':'평균 참조(tC/ha)','estimate_reference_ratio':'추정/참조'}
site_summary=build_site_summary(plot_comparison,detmetrics,matched,CFG['sites'])
summary_display=site_summary.rename(columns=SUMMARY_LABELS)
site_summary.to_csv(OUT/'site_summary.csv',index=False)

write_geo_csv(detections,OUT/'detected_trees.csv')
write_geo_csv(field_geo,OUT/'field_trees.csv')
field_all.to_csv(OUT/'field_trees_all.csv',index=False)
matched.to_csv(OUT/'matched_trees_visible.csv',index=False)
for key,df in MATCHES.items(): df.to_csv(OUT/f'matches_{key}.csv',index=False)
detmetrics.to_csv(OUT/'detection_metrics.csv',index=False)
stage_errors.to_csv(OUT/'stage_errors.csv',index=False)
plot_comparison.to_csv(OUT/'plot_comparison.csv',index=False)
sensitivity.to_csv(OUT/'height_sensitivity.csv',index=False)
pd.DataFrame(AUDIT,columns=['stage','item','reason']).to_csv(OUT/'audit.csv',index=False)
species_valid=matched[~matched.get('species_coeff_missing',pd.Series(True,index=matched.index)).astype(bool)]
species_plot=species_valid.groupby(['plotID','field_year']).agg(n=('individualID','size'),
    matched_default_AGB_kg=('AGB_default_kg','sum'),matched_known_species_AGB_kg=('AGB_known_species_kg','sum')).reset_index()
species_plot['change_pct']=100*(species_plot.matched_known_species_AGB_kg/species_plot.matched_default_AGB_kg-1)
species_plot['change_C_kg']=(species_plot.matched_known_species_AGB_kg-species_plot.matched_default_AGB_kg)*CFG['carbon_fraction']
species_plot.to_csv(OUT/'species_sensitivity_matched.csv',index=False)

def micro_detection(df):
    rows=[]
    for (year,pop),z in df.groupby(['field_year','population']):
        tp=z.TP.sum(); fp=z.FP.sum(); fn=z.FN.sum()
        rows.append(dict(field_year=year,population=pop,TP=tp,FP=fp,FN=fn,
            precision=tp/(tp+fp) if tp+fp else np.nan,recall=tp/(tp+fn) if tp+fn else np.nan,
            F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan))
    return pd.DataFrame(rows)

species=matched[~matched.get('species_coeff_missing',pd.Series(True,index=matched.index)).astype(bool)]
species_summary=species.groupby('field_year')[['AGB_default_kg','AGB_known_species_kg','species_effect_kg','species_effect_pct']].mean()
report=f"""# GRSM·ORNL 2021 나무·탄소 비교 보고서

## 나무 찾기 방법 비교 및 이후 기본값
기본 방법: **{DEFAULT_DETECTION_METHOD}**, 선택 기준: 같은 23개 조사구의 탄소 RMSE 최소.
{method_comparison.to_markdown(index=False,floatfmt='.3f')}
A/C는 바깥 GroupKFold 평가 탐지만 합친 OOF 결과다. 각 fold의 조합은 그 fold의 학습 조사구에서만 보이는 나무 F1로 선택했다.
전 표본 선택 조합은 A_deployment_grid.csv와 selected_detection_method.json에 새 자료 적용용으로 별도 저장했다.
기본 방법의 선택 오차는 동일 비교 표본의 값이며 독립적인 지역 일반화 성능으로 해석하지 않는다.
B: sigma={CFG['chm_sigma_m']}m Gaussian, 창 지름=2.51+0.00901H², H=꼭대기 픽셀 원본 CHM 높이,
CD=2√(watershed 면적/π). Popescu & Wynne (2004), [원 논문](https://doi.org/10.14358/pers.70.5.589).
매칭은 B도 수관 envelope 상자를 사용해 기존 point-in-box 1:1 규칙을 유지했다.
C: A OOF 상자 안 최고 CHM 꼭대기로 수고 보완, 상자 CD 유지; 꼭대기 없는 DF 상자는 보존;
상자 밖 CHM 꼭대기는 B의 수고·수관으로 추가, {CFG['fusion_dedup_m']}m 중복 제거.
시드 watershed는 [scikit-image 공식 API](https://scikit-image.org/docs/stable/api/skimage.segmentation.html#skimage.segmentation.watershed)를 사용한다.

## 사이트별 결과 요약
{summary_display.to_markdown(index=False,floatfmt='.3f')}
F1은 보이는 나무 기준을 기본으로 제시하고 전체 매핑 나무 F1도 함께 표시한다.
수고·DBH·개체 AGB Bias는 보이는 나무 매칭 표본 기준이다.
조사구 탄소 평균·비율은 추정과 참조가 모두 유효한 같은 조사구 집합의 평균이다.

실행 ID: {RUN_ID}; 작업 CRS: {WORK_CRS.to_string()}
분석 범위: 실제 조사 소구획의 살아 있는 DBH≥{CFG['min_field_dbh_cm']}cm. 조사구 참조 총량에는 좌표 없는 나무도 포함.
2021년 실측과 인접연도 보완을 분리하며, 보완 연도={sorted(set(plots.field_year)-{CFG['year']})}.
실행한 조사구={len(plots)}, 전체 실측 개체={len(field_all)}, 매칭 좌표 개체={len(field_geo)}, 탐지={len(detections)}, 보이는 나무 매칭={len(matched)}.

## 영상 연월·잠정 자료 여부
{pd.DataFrame([dict(siteID=s,**v) for s,v in AOP_SELECTION.items()]).to_markdown(index=False)}
1차 수정 시 API 조회에서 ORNL 2021-06은 PROVISIONAL로 확인했다. 현재 선택 상태는 위 표를 따른다.
영상 대체가 필요할 때만 ORNL 2022→2020을 검토하며 실측 연도는 유지한다.
실제 제공 연월은 aop_availability.csv, 타일 수는 aop_download_inventory.csv,
조사구별 커버 여부는 aop_plot_coverage.csv에 기록했다.

## 탐지 정확도 (마이크로 집계)
{micro_detection(detmetrics).to_markdown(index=False)}

## 단계별 오차
Bias=추정−참조; Pearson r. DBH 회귀는 조사구 GroupKFold의 OOF 결과.
{stage_errors.to_markdown(index=False)}

## 조사구 AGB·탄소 비교
{plot_comparison.to_markdown(index=False)}
총량 차이 = 매칭 개체 추정오차 + 미매칭 탐지 추정량 − 놓친 실측 개체 참조량.
참조 AGB 자체가 누락된 경우에만 총량은 NaN 처리. 면적은 totalSampledAreaTrees를 우선한다.
area_m2·reported_area_m2는 분모, polygon_area_m2·area_difference_m2는 기하 면적 진단이다.
AGB_ref_coordonly_t_ha는 좌표 있는 나무만의 보조 참조값이다.
unmapped_truth_t_ha는 좌표/경계 문제로 매칭할 수 없던 실측량이며 실제 탐지 누락과 구분해야 한다.

## 수종 미지정 영향 (계수가 확인된 동일 매칭 표본의 평균)
{species_summary.to_markdown()}
수종 미대응 개체는 이 표에서 제외하고 개별 CSV에 표시.
조사구별 동일 매칭 표본 합계 비교는 species_sensitivity_matched.csv 참고.

## 높이 민감도
{sensitivity.groupby('height_error_pct')[['change_pct','change_tC_ha']].mean().to_markdown()}
탐지와 CD를 고정한 가상 H 변화. CHM 오차가 탐지/매칭에 미치는 영향은 포함하지 않음.
탄소계수 0.47 사용 시 0.5 대비 모든 탄소량은 6% 감소.

## 공식·계수·출처
- DBH(cm)=0.367×(H(m)×CD(m))^0.887. itcSegment biome 11; g=0.
  [Jucker 구현 원문]({JUCKER_SOURCE}); Jucker et al. (2017), GCB 23:177–190.
- AGB(kg)=exp(β0+β1×ln(DBH(cm))). 기본 β0={default['beta0']}, β1={default['beta1']}.
  이는 Chojnacky의 일반 활엽수 Cor/Eri/Lau/etc 대체 그룹. ‘mixed hardwood’라는 독립 공식 아님.
  [MasterTaxonList]({TAXON_URL}), [Chojnacky 원 논문](https://doi.org/10.1093/forestry/cpt053).
- 수종 기지 시 동일 식에 taxonID별 beta0/beta1 적용; 계수 결측은 개별 표시.
- C=AGB×0.5 (민감도 0.47). kg/1000×10000/㎡ = t/ha.
- [NEONForestAGBv2]({AGB_URL}), [저장소](https://github.com/atkinsjeff/NEONForestAGB),
  [데이터 DOI](https://doi.org/10.6084/m9.figshare.25625352.v1).
- [NEON 조사 사용자 안내서](https://data.neonscience.org/api/v0/documents/NEON_vegStructure_userGuide_vE),
  [neonutilities 공식 API](https://neon-utilities-python.readthedocs.io/en/latest/api.html),
  [DeepForest 공식 예측 안내](https://deepforest.readthedocs.io/en/stable/getting_started/intro_tutorials/04_predict_large_tile.html).
  설치 버전은 versions.json, 원문 해시는 formula_provenance.json 참고.

## 해석상 한계
- NEONForestAGB는 실측 DBH 기반 모델 추정값이다. 개체 바이오매스를 파괴 조사한 실측 정답이 아니다.
- 나무 아래층은 영상에 가려진다. 개체 매칭의 ‘전체’는 DBH≥10cm 좌표 확인 나무이며 조사구 총량은 좌표 없는 나무도 포함한다.
- 줄기 위치 기반 매칭은 수관 IoU 검증을 대체하지 않으며 수평 위치 오차·기울어진 줄기·겹친 수관 영향을 받는다.
- 탐지는 중심점, 실측은 줄기점 기준으로 경계를 적용한다. 경계 근처의 외부 나무 수관 유입/유출 편향이 남는다.
- CHM 최대값은 이웃 큰 나무 높이를 담을 수 있다. 상자 기반 CD는 실제 수관 직경과 다르다.
- 현재 사전학습 가중치는 NEON 자료를 사용한 모델이다. 모델 로더/버전에 따라 2021 논문의 모델과 다를 수 있다.
- subplotID와 총 조사 면적이 일치하지 않는 조사구는 추측하지 않고 제외했다. audit.csv 확인.
- 2021 현장 조사일과 항공 촬영일은 일치하지 않을 수 있다. field date/ref_date와 영상 파일 경로를 확인한다.
- 품질 플래그는 원본 표와 field_trees.csv에 보존했다. 최신 실측 배포본과 2024년 고정 AGBv2의 개정 차이가 있을 수 있다.
- 다줄기 나무는 NEON의 줄기별 individualID를 유지한다. 중복 방문은 마지막 측정으로 정리한다.
- AGB 참조값이 없는 개체의 총량을 0으로 간주하지 않는다. 좌표 미확인 개체는 audit에 기록한다.
- coordinates_complete=False는 개체 매칭의 자료 완전성 진단이다. 조사구 참조 총량을 제외하는 조건이 아니다.
- 한국에서는 CRS·센서 해상도·수관 구조·수종식·CHM 정확도가 달라 지역 검증과 필요 시 모델 재학습이 필요하다.
"""
(OUT/'report.md').write_text(report,encoding='utf-8')
(ROOT/'outputs/report.md').write_text(f'최근 완료 결과: [{RUN_ID}/report.md]({RUN_ID}/report.md)\n',encoding='utf-8')
print('완료:',OUT/'report.md')
display(micro_detection(detmetrics))
''')
export_names={'require_columns','offset_xy','estimate_dbh','estimate_agb','kg_to_t_ha','carbon_kg','errors',
    'grouped_dbh_regression','match_one_to_one','zonal_max','predict_crowns','label_visibility',
    'crop_mosaic','raster_index','tile_coordinates','product_site_entry','select_aop_pair','aop_tiff_files',
    'raster_covers_geometry','plot_biomass_totals','comparison_statistics','build_site_summary','plot_carbon_scatter',
    'finalize_detection_frame','deepforest_patch_candidates','variable_window_peaks','chm_watershed_crowns',
    'combine_boxes_and_peaks','detection_counts','nested_group_parameter_selection','evaluate_detection_method','method_comparison_rows',
    'detection_algorithm_smoke_checks'}
module_header='''"""Reusable inference/evaluation functions. Requires Colab scientific dependencies.
Call configure(config, output_directory, projected_crs) before inference.
DBH inputs: H, CD in metres. AGB input DBH in centimetres; output kilograms.
"""
from pathlib import Path
import inspect
import math
import json, importlib.metadata
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.mask import mask
from shapely.geometry import Point, box, Polygon
from shapely.ops import unary_union
from scipy.optimize import linear_sum_assignment
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, maximum_filter, label as ndi_label
from skimage.segmentation import watershed
from rasterio.features import shapes as raster_shapes
from shapely.geometry import shape
AUDIT=[]
def configure(config, output_directory, projected_crs):
    global CFG, OUT, WORK_CRS, default
    CFG=config
    OUT=Path(output_directory)
    OUT.mkdir(parents=True,exist_ok=True)
    WORK_CRS=projected_crs
    default=config['agb_default']
def audit(stage,item,reason):
    AUDIT.append(dict(stage=stage,item=str(item),reason=str(reason)))
'''
exports=[]
for c in cells:
    if c['cell_type']!='code' or ''.join(c['source']).startswith('%'): continue
    for node in ast.parse(''.join(c['source'])).body:
        if isinstance(node,ast.FunctionDef) and node.name in export_names:
            exports.append(ast.unparse(node))
assert len(exports)==len(export_names)
module_source=module_header+'\n\n'+'\n\n'.join(exports)+'\n'
ast.parse(module_source)
(HERE/'neon_pipeline.py').write_text(module_source,encoding='utf-8')
section('18. 재사용 함수 모듈 내보내기', '''탐지·공식·매칭·회귀·평가 함수를 `neon_pipeline.py`로 저장합니다.
다른 Colab에서 모듈을 불러온 후 `configure(설정dict, 결과폴더, 작업CRS)`를 호출하세요.
`estimate_dbh`와 `estimate_agb`에는 공식 계수를 직접 전달할 수도 있습니다.''',
    "(OUT/'neon_pipeline.py').write_text("+repr(module_source)+",encoding='utf-8')\nprint(OUT/'neon_pipeline.py')")

md('''## 19. 한국 자료로 재사용하기

1. 설정의 `WORK_CRS`를 한국 입력의 미터 투영 CRS로 바꾸세요.
2. `RASTERS={plotID:(RGB경로, CHM경로)}`, `plots`(실제 조사 geometry·area_m2·siteID·field_year·중심 좌표),
   `field_all`(좌표 유무와 무관한 전체 실측·참조 AGB)과 `field_geo`(매칭 가능한 실측 점)를 같은 필드로 준비하면
   NEON 다운로드 단계를 건너뛸 수 있습니다. `plots.imagery_year`로 영상 연도를 별도로 지정하세요.
3. Jucker 계수와 AGB 계수를 한국 수종·지역 검증식으로 교체하세요. 북미 계수를 그대로 적용한 값은 탐색적 추정입니다.
4. `predict_crowns`, `estimate_dbh`, `estimate_agb`, `label_visibility`, `match_one_to_one`,
   `grouped_dbh_regression`, `errors`, `kg_to_t_ha`를 재사용하세요.
5. 실측 DBH≥10cm라는 조사 범위와 탐지 가능한 모든 수관의 모집단 차이를 고려하세요.
   정확한 전체 산림 탄소 추정에는 작은 나무·관목·아래층을 포함하는 별도 설계가 필요합니다.

**결과 확인 순서:** audit.csv → sampled_footprints.gpkg/overlay 그림 → detection_metrics.csv →
stage_errors.csv → plot_comparison.csv → height_sensitivity.csv → report.md.

이 노트북은 로컬에서 문법·구조를 검사했으며 Colab에서의 다운로드/GPU 실행 결과는 포함하지 않습니다.
''')

section('20. 사람이 읽는 최종 요약', '''사이트별 조사구·전체 실측·매칭 좌표 나무·탐지 수와 정확도를 한 표로 봅니다.
보고서 맨 위와 동일한 요약이며, 영상 연도와 실측 연도가 다르면 행을 분리합니다.
Bias는 추정−참조, F1 기본은 보이는 나무, 탄소 비율은 동일한 유효 조사구들의 평균 추정/평균 참조입니다.''', r'''
display(summary_display.round(3))
display(method_comparison.round(3))
print('이후 기본 나무 찾기 방법:',DEFAULT_DETECTION_METHOD)
print('전체 참조값에는 좌표 없는 나무가 포함됩니다. 자세한 자료는',OUT/'field_trees_all.csv')
print('보고서:',OUT/'report.md')
''')
section('21. 결과 ZIP 저장 및 다운로드', '''현재 실행 결과 폴더만 ZIP으로 묶습니다.
`/content/neon_results.zip`을 만들고 Drive의 outputs/에도 사본을 남긴 뒤 브라우저 다운로드를 시작합니다.''', r'''
import shutil
from google.colab import files
archive=shutil.make_archive('/content/neon_results','zip',root_dir=str(OUT))
shutil.copy2(archive,ROOT/'outputs/neon_results.zip')
print('ZIP:',archive,'/ Drive 사본:',ROOT/'outputs/neon_results.zip')
files.download('/content/neon_results.zip')
''')

from resolution_extension import apply_resolution_extension
cells = apply_resolution_extension(cells)

# Every code cell has Korean explanatory Markdown directly above it.
for i,c in enumerate(cells):
    if c['cell_type']=='code':
        src=''.join(c['source'])
        if not src.startswith('%'): ast.parse(src,filename=f'cell_{i}')
        assert i>0 and cells[i-1]['cell_type']=='markdown'
nb=dict(nbformat=4,nbformat_minor=5,metadata=dict(
    kernelspec=dict(display_name='Python 3',language='python',name='python3'),
    language_info=dict(name='python',version='3.11'),
    colab=dict(name='NEON_DeepForest_AGB_2021.ipynb',provenance=[]),accelerator='GPU'),cells=cells)
for i,c in enumerate(cells): c['id']=f'neon-{i:03}'
path=HERE/'NEON_DeepForest_AGB_2021.ipynb'
path.write_text(json.dumps(nb,ensure_ascii=False,indent=1),encoding='utf-8')
print(path)
print(f'{len(cells)} cells; {sum(c["cell_type"]=="code" for c in cells)} code cells; syntax OK')
