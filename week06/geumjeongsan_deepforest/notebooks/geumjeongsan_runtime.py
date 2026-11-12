"""Colab scientific runtime, embedded in GEUMJEONGSAN_DeepForest_2021.ipynb.
Never imported/executed by the local Windows validation script.
"""
from contextlib import ExitStack
from types import SimpleNamespace
from xml.sax.saxutils import escape
import importlib.util, inspect, json, os, re, shutil, tempfile, time, warnings, gzip
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.enums import Resampling
from rasterio.plot import show
from rasterio.warp import reproject, transform_bounds
from rasterio.transform import from_origin
from rasterio.errors import NotGeoreferencedWarning
from shapely.geometry import Point, box, Polygon, shape
from shapely.ops import unary_union
from shapely import wkt
from pyproj import Transformer, CRS
from sklearn.model_selection import GroupKFold
import matplotlib.pyplot as plt
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

CRS_WORK=CRS.from_epsg(5187)
TO_XY=Transformer.from_crs(4326,5187,always_xy=True)
TO_LL=Transformer.from_crs(5187,4326,always_xy=True)
JUCKER_URL='https://raw.githubusercontent.com/cran/itcSegment/master/R/dbh.R'
META_PREFIX='forests/v1/alsgedi_global_v6_float/'
META_BUCKET='https://dataforgood-fb-data.s3.amazonaws.com'
SCHEMA='geumjeongsan2021-v1'
AUDIT=[]

def audit(stage,item,message):
    AUDIT.append(dict(stage=stage,item=str(item),message=str(message)))

def json_save(path,obj):
    path=Path(path);tmp=path.with_name(path.name+'.part')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str),encoding='utf-8');tmp.replace(path)

def geo_csv(frame,path):
    frame=frame.copy()
    if 'geometry' in frame:frame['geometry']=frame.geometry.map(lambda g:g.wkt if g is not None else '')
    pd.DataFrame(frame).to_csv(path,index=False)

def geo_read(path):
    df=pd.read_csv(path);df['geometry']=df.geometry.map(wkt.loads)
    return gpd.GeoDataFrame(df,geometry='geometry',crs=CRS_WORK)

def stat_signature(path):
    path=Path(path);s=path.stat()
    return dict(path=str(path),bytes=s.st_size,mtime_ns=s.st_mtime_ns)

def cache_valid(stamp,fingerprint,files):
    return all(Path(p).is_file() for p in files) and Path(stamp).is_file() and json.loads(Path(stamp).read_text())==fingerprint

def download(url,dest,optional=False):
    dest=Path(dest)
    if dest.exists() and dest.stat().st_size:return dest
    dest.parent.mkdir(parents=True,exist_ok=True);tmp=dest.with_name(dest.name+'.part')
    for attempt in range(5):
        try:
            with SESSION.get(url,stream=True,timeout=(30,180)) as response:
                if response.status_code==404 and optional:return None
                response.raise_for_status()
                with tmp.open('wb') as f:
                    for block in response.iter_content(1024*1024):
                        if block:f.write(block)
                if not response.headers.get('Content-Encoding') and response.headers.get('Content-Length'):
                    if tmp.stat().st_size!=int(response.headers['Content-Length']):raise IOError('전송 크기 불일치')
            tmp.replace(dest);return dest
        except Exception:
            if attempt==4:
                if optional:return None
                raise
            time.sleep(min(2**attempt,8))

def initialize():
    global DATA,OUT,CACHE,SESSION,NEON,JUCKER,ALLO,PLOTS,LAYERS,INPUT_SIGNATURE,MODEL
    MODEL=None
    DATA=ROOT/'gedi_geumjeongsan/data'
    if not DATA.exists() and (ROOT/'data').exists():DATA=ROOT/'data'
    needed=[DATA/'field'/n for n in ('seo2011_plots.csv','seo2011_trees.csv','seo2011_layers.csv')]
    needed += [DATA/'reference'/n for n in ('README.md','nifos2014_allometry.csv','nifos2014_carbon_factors.csv','nifos2010_carbon_fraction.csv')]
    missing=[str(p) for p in needed if not p.is_file()]
    if missing:raise FileNotFoundError('Drive 업로드 파일 누락:\n'+'\n'.join(missing))
    print((DATA/'reference/README.md').read_text(encoding='utf-8'))
    for module in (np,pd,gpd,rasterio):print(module.__name__,module.__version__)
    import importlib.metadata
    assert importlib.metadata.version('deepforest')=='2.1.0'
    SESSION=requests.Session();SESSION.mount('https://',HTTPAdapter(max_retries=Retry(total=5,backoff_factor=1,status_forcelist=(429,500,502,503,504))))
    CACHE=ROOT/'data/cache/geumjeongsan';CACHE.mkdir(parents=True,exist_ok=True)
    source=download(JUCKER_URL,CACHE/'reference/itcSegment_dbh.R')
    JUCKER=parse_jucker_r(source.read_text(encoding='utf-8'))
    for biome,expected in ((20,dict(a=.694,b=.73,g=0)),(21,dict(a=1.004,b=.73,g=0))):
        if JUCKER[biome]!=expected:warnings.warn(f'Jucker biome{biome} 원문값 변경: {JUCKER[biome]}; 원문값을 사용합니다.')
    ALLO=KoreanAllometry(read_csv_rows(DATA/'reference/nifos2014_allometry.csv'),read_csv_rows(DATA/'reference/nifos2014_carbon_factors.csv'),read_csv_rows(DATA/'reference/nifos2010_carbon_fraction.csv'))
    plot_table=pd.read_csv(DATA/'field/seo2011_plots.csv');assert len(plot_table)==50 and plot_table.plot_id.is_unique
    assert set(plot_table.zone)=={'중앙','동측','서측','남측','북측'} and plot_table.groupby('zone').size().eq(10).all()
    x,y=TO_XY.transform(plot_table['경도'].to_numpy(),plot_table['위도'].to_numpy())
    PLOTS=gpd.GeoDataFrame(plot_table,geometry=gpd.points_from_xy(x,y),crs=CRS_WORK)
    PLOTS['plotID']=PLOTS.plot_id;PLOTS['center_x']=x;PLOTS['center_y']=y
    PLOTS['area_m2']=400.;PLOTS['field_period']='2009-05..2010-09'
    LAYERS=pd.read_csv(DATA/'field/seo2011_layers.csv')
    t1=LAYERS[LAYERS['기호']=='T1'];assert len(t1)==50 and t1.plot_id.is_unique
    PLOTS['T1_min_m']=PLOTS.plot_id.map(t1.set_index('plot_id')['최소수고_m'])
    assert PLOTS.T1_min_m.notna().all()
    ortho=DATA/'ortho/2021'
    for sheet in SHEETS:
        path=ortho/f'(B060)정사영상_2021_{sheet}.tif';xmls=list(ortho.glob('*'+sheet+'.xml'))
        if not path.is_file() or len(xmls)!=1:raise FileNotFoundError(f'도엽{sheet} TIF/XML 업로드 필요')
        needed += [xmls[0]]
    # Pin the verified NEON resolution model; an optional uploaded provenance must agree.
    candidates=[ROOT/'neon_grsm_ornl/outputs_v4_resolution/resolution_model_provenance.json',ROOT/'resolution_model_provenance.json']
    for path in candidates:
        if path.exists():
            actual=json.loads(path.read_text())
            if actual['resolved_revision']!=HF_MODEL_PROVENANCE['resolved_revision']:raise ValueError('업로드한 HF 커밋과 내장 NEON 커밋이 다릅니다. 설정을 검토하세요.')
    reference_hashes={str(p.relative_to(DATA)):sha_file(p) for p in needed}
    spec=dict(schema=SCHEMA,code_sha=GEUM_CODE_SHA,files=reference_hashes,jucker_sha=sha_file(source),
              hf=HF_MODEL_PROVENANCE,neon_25cm=NEON_25CM_BASELINE,
              vectors=[stat_signature(p) for category in ('2025_임상도','등산로')
                       for p in sorted((DATA/'forest_map'/category).rglob('*')) if p.suffix.lower() in ('.shp','.dbf','.shx','.prj','.cpg')],
              ortho={s:stat_signature(ortho/f'(B060)정사영상_2021_{s}.tif') for s in SHEETS})
    INPUT_SIGNATURE=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    OUT=ROOT/'outputs'/('geumjeongsan_'+INPUT_SIGNATURE[:12]);OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'figures').mkdir(exist_ok=True)
    json_save(OUT/'provenance.json',dict(spec,run_id=OUT.name,Jucker={k:JUCKER[k] for k in (20,21)},
               neon_pipeline_sha256=hashlib.sha256(NEON_PIPELINE_SOURCE.encode()).hexdigest(),
               units='D cm; H/CD m; AGB dry kg; aboveground C only; tC/ha',
               no_roots=True,field_period='2009-05..2010-09',image_year=2021))
    pd.DataFrame([dict(path=k,sha256=v) for k,v in reference_hashes.items()]).to_csv(OUT/'input_manifest.csv',index=False)
    # Self-contained notebook: no local module upload is required.
    module_dir=Path('/content/geumjeongsan_runtime');module_dir.mkdir(exist_ok=True)
    module_path=module_dir/'neon_pipeline.py';module_path.write_text(NEON_PIPELINE_SOURCE,encoding='utf-8')
    spec_module=importlib.util.spec_from_file_location('geum_neon_pipeline',module_path)
    NEON=importlib.util.module_from_spec(spec_module);spec_module.loader.exec_module(NEON)
    dimensional_checks()
    print('지상부만 계산; 뿌리 제외. 결과:',OUT,'\nHF 고정 커밋:',HF_MODEL_PROVENANCE['resolved_revision'])

def georeference():
    global ORTHO,AOI
    AOI=box(PLOTS.center_x.min()-80,PLOTS.center_y.min()-80,PLOTS.center_x.max()+80,PLOTS.center_y.max()+80)
    ortho=DATA/'ortho/2021';directory=CACHE/'vrt';directory.mkdir(exist_ok=True)
    ORTHO=[];rows=[]
    to5179=Transformer.from_crs(5187,5179,always_xy=True)
    for sheet in SHEETS:
        path=ortho/f'(B060)정사영상_2021_{sheet}.tif';xml_path=next(ortho.glob('*'+sheet+'.xml'))
        xml=ET.parse(xml_path).getroot()
        assert float(xml.findtext('.//지상표본거리'))==.25
        assert xml.findtext('.//투영원점')=='동부' and xml.findtext('.//기준타원체')=='GRS80'
        assert float(xml.findtext('.//Easting'))==200000 and float(xml.findtext('.//Northing'))==600000
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',NotGeoreferencedWarning)
            with rasterio.open(path) as ds:
                width,height=ds.width,ds.height
                if ds.count<3 or ds.dtypes[0]!='uint8':raise ValueError('8bit RGB 필요')
        w,s,e,n=sheet_bounds(sheet)
        coordinates=[TO_XY.transform(lon,lat) for lon in (w,e) for lat in (s,n)]
        minx=min(x for x,y in coordinates);maxx=max(x for x,y in coordinates)
        miny=min(y for x,y in coordinates);maxy=max(y for x,y in coordinates)
        left,top=minx-50,maxy+50;right=left+width*.25;bottom=top-height*.25
        transform=from_origin(left,top,.25,.25)
        row=dict(sheet=sheet,width_px=width,height_px=height,gsd_m=.25,west=w,south=s,east=e,north=n,
                 left_5187=left,bottom_5187=bottom,right_5187=right,top_5187=top,
                 width_difference_m=width*.25-(maxx-minx+100),height_difference_m=height*.25-(maxy-miny+100),
                 source_crs='EPSG:5187 manual assignment',georef_method='sheet bbox +50m border, north-up',xml_sha=sha_file(xml_path))
        assert abs(row['width_difference_m'])<5 and abs(row['height_difference_m'])<5
        points=[to5179.transform(x,y) for x in (left,right) for y in (bottom,top)]
        mapped=(min(x for x,y in points),min(y for x,y in points),max(x for x,y in points),max(y for x,y in points))
        for name,value in zip(('left','bottom','right','top'),mapped):row[name+'_5179']=value
        if sheet=='35913003':
            errors=[a-b for a,b in zip(mapped,CONTROL_5179)]
            for name,err in zip(('left','bottom','right','top'),errors):row['control_error_'+name+'_m']=err
            row['control_max_abs_error_m']=max(map(abs,errors))
            assert row['control_max_abs_error_m']<=3,f'35913003 좌표 검증 실패: {errors}'
        # A VRT assigns metadata in the cache; the original TIFF is never modified.
        vrt=directory/(sheet+'.vrt')
        bands=''.join(f'''<VRTRasterBand dataType="Byte" band="{b}"><ColorInterp>{color}</ColorInterp><SimpleSource>
<SourceFilename relativeToVRT="0">{escape(str(path))}</SourceFilename><SourceBand>{b}</SourceBand>
<SrcRect xOff="0" yOff="0" xSize="{width}" ySize="{height}"/><DstRect xOff="0" yOff="0" xSize="{width}" ySize="{height}"/>
</SimpleSource></VRTRasterBand>''' for b,color in enumerate(('Red','Green','Blue'),1))
        vrt.write_text(f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}"><SRS>EPSG:5187</SRS><GeoTransform>{left},0.25,0,{top},0,-0.25</GeoTransform>{bands}</VRTDataset>',encoding='utf-8')
        ORTHO.append(dict(sheet=sheet,path=path,vrt=vrt,geometry=box(left,bottom,right,top),transform=transform,width=width,height=height))
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT/'georef_check.csv',index=False);display(pd.DataFrame(rows))
    coverage=unary_union([o['geometry'] for o in ORTHO])
    for _,p in PLOTS.iterrows():
        if not coverage.buffer(.001).covers(box(p.center_x-60,p.center_y-60,p.center_x+60,p.center_y+60)):
            raise ValueError(f'{p.plot_id} 60m crop 전체를 덮는 영상이 없습니다. 조사구를 임의 제외하지 않습니다.')
    return pd.DataFrame(rows)

def read_vectors():
    global FOREST,TRAILS
    bounds5179=transform_bounds(5187,5179,*AOI.bounds,densify_pts=41)
    forest=[];trails=[]
    for category,target in (('2025_임상도',forest),('등산로',trails)):
        for district in ('금정구','북구','양산시'):
            paths=sorted((DATA/'forest_map'/category/district).glob('*.shp'))
            if not paths:raise FileNotFoundError(f'필수 shapefile 없음: {category}/{district}')
            for path in paths:
                df=gpd.read_file(path,bbox=bounds5179,encoding='cp949')
                print(category,district,'속성:',list(df.columns));display(df.drop(columns='geometry').head(3))
                if df.crs is None:
                    warnings.warn(f'{path.name} PRJ 없음: 입력 지시 EPSG5179 사용');df=df.set_crs(5179)
                if category=='2025_임상도':
                    if not ({'FRTP_NM','FRTP_CD'}&set(df)) or not ({'KOFTR_NM','KOFTR_GROU'}&set(df)):
                        raise ValueError('임상/수종 필드 확인 필요; 임의 필드로 진행하지 않습니다.')
                df=df.to_crs(CRS_WORK);df=df[df.geometry.notna()&~df.geometry.is_empty].copy()
                invalid=~df.geometry.is_valid
                if invalid.any():
                    audit('geometry_repair',path,f'{int(invalid.sum())}개 make_valid')
                    df.loc[invalid,'geometry']=df.loc[invalid,'geometry'].make_valid()
                df['source_district']=district;df['source_file']=path.name
                target.append(df[df.geometry.intersects(AOI)].copy())
    FOREST=gpd.GeoDataFrame(pd.concat(forest,ignore_index=True),geometry='geometry',crs=CRS_WORK)
    TRAILS=gpd.GeoDataFrame(pd.concat(trails,ignore_index=True),geometry='geometry',crs=CRS_WORK)
    assert len(FOREST)>0
    fields=pd.DataFrame([dict(layer=label,field=c) for label,df in [('forest',FOREST),('trails',TRAILS)] for c in df.columns])
    fields.to_csv(OUT/'vector_fields.csv',index=False)

def text(value):
    return '' if value is None or pd.isna(value) else str(value).strip()

SPECIES_CODES={'11':'소나무','12':'잣나무','13':'낙엽송','14':'리기다소나무','15':'곰솔','17':'편백나무','18':'삼나무',
 '10':'기타침엽수','31':'상수리나무','32':'신갈나무','33':'굴참나무','34':'기타참나무류','35':'오리나무',
 '36':'고로쇠나무','37':'자작나무','38':'박달나무','39':'밤나무','40':'물푸레나무','41':'서어나무',
 '42':'때죽나무','44':'백합나무','45':'포플러','46':'벚나무','47':'느티나무','49':'아까시나무','30':'기타활엽수','77':'침활혼효림','78':'죽림'}
SPECIES_CODES.update({'16':'잔나무','19':'가문비나무','20':'비자나무','21':'은행나무','43':'호두나무','48':'층층나무',
 '61':'가시나무','62':'구실잣밤나무','63':'녹나무','64':'굴거리나무','65':'황칠나무','66':'사스레피나무','67':'후박나무',
 '68':'새덕이','60':'기타상록활엽수','81':'미립목지','82':'제지','83':'관목덤불','91':'주거지','92':'초지','93':'경작지','94':'수체','95':'과수원','99':'기타'})

def forest_at(point):
    hits=list(FOREST.sindex.query(point,predicate='intersects'))
    def priority(i):
        r=FOREST.iloc[i];return (text(r.get('FRTP_CD')).removesuffix('.0') not in ('1','2','3'),('금정구','북구','양산시').index(r.source_district),int(i))
    if not hits:return dict(forest_code='',forest_type='임상도 밖/미피복',map_species='',forest_hit_count=0,forest_fallback=True)
    r=FOREST.iloc[min(hits,key=priority)]
    code=text(r.get('FRTP_CD')).removesuffix('.0');name=text(r.get('FRTP_NM'))
    if not code:
        code={'침엽수림':'1','활엽수림':'2','혼효림':'3','침활혼효림':'3','죽림':'4'}.get(name,'')
    code_species=text(r.get('KOFTR_GROU')).removesuffix('.0')
    species=text(r.get('KOFTR_NM')) or SPECIES_CODES.get(code_species,'')
    return dict(forest_code=code,forest_type=name or {'1':'침엽수림','2':'활엽수림','3':'혼효림','4':'죽림','0':'무립목지/비산림'}.get(code,'임상미상'),
                map_species=species,forest_hit_count=len(hits),forest_fallback=code not in ('1','2','3'),
                forest_source_district=r.source_district,forest_source_file=r.source_file)

def field_reference():
    rows=[];lookup=PLOTS.set_index('plot_id')
    for i,t in enumerate(read_csv_rows(DATA/'field/seo2011_trees.csv'),1):
        p=lookup.loc[t['plot_id']];d=float(t['흉고직경_cm']);h=float(t['수고_m'])
        included=t['수종']!='등나무' and d>=6
        estimate=ALLO.predict(t['수종'],d) if included else dict(AGB_kg=np.nan,carbon_kg=np.nan,carbon_cf_kg=np.nan,formula='excluded',mapping='등나무/DBH<6 제외',range_flag='',CF=.5)
        r=dict(tree_id=i,**t,zone=p.zone,DBH_cm=d,H_m=h,included_ge6=included,included_ge10=included and d>=10,
               exclusion_reason='등나무' if t['수종']=='등나무' else 'DBH<6cm' if d<6 else '',
               T1_min_m=p.T1_min_m,visible_by_T1=h>=p.T1_min_m,field_area_m2=400,**estimate)
        rows.append(r)
    result=pd.DataFrame(rows);assert len(result)==1750
    result.to_csv(OUT/'field_reference_trees.csv',index=False)
    return result

def rgb_crop(plot):
    directory=CACHE/'crops'/INPUT_SIGNATURE[:16]/plot.plot_id;directory.mkdir(parents=True,exist_ok=True)
    original=directory/'rgb_25cm.tif';expanded=directory/'rgb_10cm_bilinear.tif';stamp=directory/'rgb.json'
    region=box(plot.center_x-60,plot.center_y-60,plot.center_x+60,plot.center_y+60)
    selected=[r for r in ORTHO if r['geometry'].intersects(region)]
    fingerprint=dict(schema=SCHEMA,bounds=list(region.bounds),sources=[stat_signature(o['path']) for o in selected],
                     transforms=[list(o['transform']) for o in selected],resampling='bilinear to exact original extent at10cm')
    if not cache_valid(stamp,fingerprint,[original,expanded]):
        with ExitStack() as stack:
            sources=[stack.enter_context(rasterio.open(o['vrt'])) for o in selected]
            data,trans=merge(sources,bounds=region.bounds,res=.25,indexes=[1,2,3])
        assert data.shape==(3,480,480)
        profile=dict(driver='GTiff',height=480,width=480,count=3,dtype='uint8',crs=CRS_WORK,transform=trans,compress='deflate',tiled=True,blockxsize=256,blockysize=256)
        tmp=directory/'rgb_25cm.part.tif'
        with rasterio.open(tmp,'w',**profile) as ds:ds.write(data)
        tmp.replace(original)
        with rasterio.open(original) as ds:up=ds.read(out_shape=(3,1200,1200),resampling=Resampling.bilinear)
        profile.update(height=1200,width=1200,transform=trans*rasterio.Affine.scale(.4,.4))
        tmp=directory/'rgb_10cm.part.tif'
        with rasterio.open(tmp,'w',**profile) as ds:ds.write(up)
        tmp.replace(expanded);json_save(stamp,fingerprint)
    with rasterio.open(original) as lo,rasterio.open(expanded) as hi:
        assert lo.crs==hi.crs==CRS_WORK and lo.width==lo.height==480 and hi.width==hi.height==1200
        assert np.allclose(lo.res,.25) and np.allclose(hi.res,.1)
        assert np.allclose(lo.bounds,hi.bounds,atol=1e-8,rtol=0)
    return original,expanded,directory

def list_s3(prefix,delimiter=None,max_keys=100):
    params={'list-type':'2','prefix':prefix,'max-keys':max_keys}
    if delimiter:params['delimiter']=delimiter
    response=SESSION.get(META_BUCKET,params=params,timeout=60);response.raise_for_status()
    root=ET.fromstring(response.content)
    return [e.text for e in root.findall('.//{*}Contents/{*}Key')],[e.text for e in root.findall('.//{*}CommonPrefixes/{*}Prefix')]

def meta_sources():
    global META_SOURCES,META_METADATA,META_SIGNATURES
    bounds=transform_bounds(5187,4326,*AOI.bounds,densify_pts=41)
    quadkeys=quadkeys_for_bounds(bounds,9);print('금정산 AOI quadkey(zoom9):',quadkeys)
    sources={};metadata={};signatures={}
    search_dirs=[DATA/'chm',DATA/'raw/meta_chm',ROOT/'data/raw/meta_chm',ROOT/'gedi_hongcheon/data/raw/chm']+[Path(p) for p in EXTRA_CHM_DIRS]
    try:
        root_keys,prefixes=list_s3(META_PREFIX,delimiter='/')
        date_prefixes=[p for p in prefixes if any(k in p.lower().removeprefix(META_PREFIX) for k in ('date','geojson','metadata','footprint'))]
    except Exception as exc:
        root_keys=[];date_prefixes=[];audit('CHM_date_listing','S3',exc)
    for quad in quadkeys:
        cached=[p for d in search_dirs if d.exists() for p in d.rglob('*'+quad+'*.tif')]
        selected=None;target=box(*bounds).intersection(box(*quadkey_bounds(quad)))
        for candidate in sorted(cached):
            try:
                with rasterio.open(candidate) as ds:
                    footprint=box(*transform_bounds(ds.crs,4326,*ds.bounds,densify_pts=21))
                    if ds.crs.is_projected and max(ds.res)<=2 and footprint.buffer(1e-5).covers(target):
                        selected=candidate;break
            except Exception as exc:audit('CHM_cached_header',candidate,exc)
        sources[quad]=str(selected) if selected else META_BUCKET+'/'+META_PREFIX+'chm/'+quad+'.tif'
        if selected:signatures[quad]=stat_signature(selected)
        else:
            response=SESSION.head(sources[quad],timeout=60);response.raise_for_status()
            signatures[quad]=dict(url=sources[quad],etag=response.headers.get('ETag'),last_modified=response.headers.get('Last-Modified'),
                                  content_length=response.headers.get('Content-Length'))
        info=dict(quadkey=quad,source=sources[quad],date_status='unknown; no observation metadata found',date_files=[])
        candidates=[k for k in root_keys if quad in k and k.endswith(('.geojson','.geojson.gz','.json','.json.gz'))]
        for prefix in date_prefixes:
            try:
                keys,_=list_s3(prefix+quad,max_keys=20)
                candidates += [k for k in keys if k.endswith(('.geojson','.json','.geojson.gz','.json.gz'))]
            except Exception as exc:audit('CHM_date',quad,exc)
        for key in sorted(set(candidates)):
            path=download(META_BUCKET+'/'+key,CACHE/'meta_dates'/Path(key).name,optional=True)
            if path:
                info['date_files'].append(str(path));info['date_status']='associated metadata downloaded; acquisition fields checked per plot'
        metadata[quad]=info
    META_SOURCES=sources;META_METADATA=metadata;META_SIGNATURES=signatures
    for q in metadata:metadata[q]['source_signature']=signatures[q]
    json_save(OUT/'meta_chm_sources.json',metadata)
    return sources

def chm_crop(plot,directory):
    native=directory/'meta_chm_native.tif';local=directory/'meta_chm_5187_1m.tif';stamp=directory/'chm.json'
    region=box(plot.center_x-60,plot.center_y-60,plot.center_x+60,plot.center_y+60)
    fingerprint=dict(schema=SCHEMA,bounds=list(region.bounds),sources=META_SOURCES,source_signatures=META_SIGNATURES)
    if not cache_valid(stamp,fingerprint,[native,local]):
        # COG HTTP range windows only; never download a global quad in full.
        for attempt in range(5):
            try:
                with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_MAX_RETRY=5,GDAL_HTTP_RETRY_DELAY=2,VSI_CACHE=True):
                    with ExitStack() as stack:
                        opened=[(q,stack.enter_context(rasterio.open(path))) for q,path in META_SOURCES.items()]
                        crs=opened[0][1].crs
                        assert all(ds.crs==crs for q,ds in opened)
                        bounds=transform_bounds(CRS_WORK,crs,*region.bounds,densify_pts=41)
                        selected=[ds for q,ds in opened if box(*ds.bounds).intersects(box(*bounds))]
                        if not selected:raise ValueError('금정산 범위를 덮는 Meta CHM 없음')
                        coverage=unary_union([box(*ds.bounds) for ds in selected])
                        if not coverage.buffer(.01).covers(box(*bounds)):raise ValueError('CHM source coverage 불완전')
                        arr,trans=merge(selected,bounds=bounds,indexes=[1],nodata=np.nan,target_aligned_pixels=True)
                        tags={q:dict(tags=ds.tags(),band_tags=ds.tags(1),crs=ds.crs.to_string(),res=list(ds.res)) for q,ds in opened}
                        dtype=str(arr.dtype);assert np.issubdtype(arr.dtype,np.floating),'float형 Meta CHM 필요'
                break
            except Exception:
                if attempt==4:raise
                time.sleep(min(2**attempt,8))
        arr=arr.astype('float32');arr[(arr<0)|~np.isfinite(arr)]=np.nan
        profile=dict(driver='GTiff',count=1,height=arr.shape[1],width=arr.shape[2],crs=crs,transform=trans,dtype='float32',nodata=np.nan,compress='deflate')
        tmp=directory/'chm_native.part.tif'
        with rasterio.open(tmp,'w',**profile) as ds:ds.write(arr)
        tmp.replace(native)
        # Also deliver a local-metric 1m crop for inspection; inference uses native
        # pixels so zonal maximum H is not changed by CHM reprojection.
        out=np.full((120,120),np.nan,'float32');local_transform=from_origin(plot.center_x-60,plot.center_y+60,1,1)
        reproject(arr[0],out,src_transform=trans,src_crs=crs,dst_transform=local_transform,dst_crs=CRS_WORK,
                  src_nodata=np.nan,dst_nodata=np.nan,resampling=Resampling.max)
        profile.update(height=120,width=120,crs=CRS_WORK,transform=local_transform)
        tmp=directory/'chm_5187.part.tif'
        with rasterio.open(tmp,'w',**profile) as ds:ds.write(out,1)
        tmp.replace(local)
        json_save(directory/'chm_source_tags.json',tags);json_save(stamp,fingerprint)
    with rasterio.open(native) as ds:
        # Every comparison pixel must be valid (including valid zero height).
        geom=gpd.GeoSeries([plot.geometry.buffer(31)],crs=CRS_WORK).to_crs(ds.crs).iloc[0]
        from rasterio.mask import mask
        data,masked_transform=mask(ds,[geom],crop=True,filled=False,indexes=1)
        from rasterio.features import geometry_mask
        inside=geometry_mask([geom],out_shape=data.shape,transform=masked_transform,invert=True)
        valid=(~np.ma.getmaskarray(data))&np.isfinite(data.data)
        valid_fraction=float(np.count_nonzero(inside&valid)/np.count_nonzero(inside))
        if valid_fraction<.999:raise ValueError(f'{plot.plot_id} 30m 비교범위 CHM미피복 {100*(1-valid_fraction):.3f}%; 탄소0으로처리하지않습니다.')
    dates=[];file_dates=[]
    tag_file=directory/'chm_source_tags.json'
    if tag_file.exists():
        for quad,payload in json.loads(tag_file.read_text()).items():
            META_METADATA[quad]['native_raster_metadata']=payload
            for namespace in ('tags','band_tags'):
                for k,v in payload.get(namespace,{}).items():
                    if any(t in k.lower() for t in ('date','year','time')):
                        if any(t in k.lower() for t in ('acq','observ','image','capture','sensing')):dates.append(f'tag:{k}={v}')
                        else:file_dates.append(f'tag:{k}={v}')
    region_ll=gpd.GeoSeries([region],crs=CRS_WORK).to_crs(4326).iloc[0]
    for quad,info in META_METADATA.items():
        for path in info['date_files']:
            try:
                with Path(path).open('rb') as stream:is_gzip=stream.read(2)==b'\x1f\x8b'
                obj=json.loads(gzip.open(path,'rt',encoding='utf-8').read() if is_gzip else Path(path).read_text(encoding='utf-8'))
                for f in obj.get('features',[]):
                    if f.get('geometry') and not shape(f['geometry']).intersects(region_ll):continue
                    for k,v in f.get('properties',{}).items():
                        if any(t in k.lower() for t in ('date','year','time','acq')):
                            dates.append(f'{k}={v}')
            except Exception as exc:audit('CHM_date_parse',path,exc)
    info=dict(plot_id=plot.plot_id,native_path=str(native),local_5187_path=str(local),
              observation_dates=';'.join(sorted(set(dates))) or 'unknown',valid_fraction=valid_fraction,
              file_timestamps_not_acquisition=';'.join(sorted(set(file_dates))),
              height_extraction='maximum of native Meta pixels; 5187 1m crop is display-only')
    return native,local,info

def build_crops():
    global RASTERS,CHM_INFO
    RASTERS={};rows=[]
    meta_sources()
    for _,p in PLOTS.iterrows():
        low,high,directory=rgb_crop(p);chm,display_chm,info=chm_crop(p,directory)
        RASTERS[p.plot_id]=(low,high,chm);rows.append(dict(info,rgb_25cm=str(low),rgb_10cm=str(high)))
        print(p.plot_id,'crop 준비:',low.name,high.name,'CHM 날짜:',info['observation_dates'])
    CHM_INFO=pd.DataFrame(rows);CHM_INFO.to_csv(OUT/'raster_manifest.csv',index=False)
    json_save(OUT/'meta_chm_sources.json',META_METADATA)

def load_model():
    global MODEL
    if globals().get('MODEL') is not None:return MODEL
    from deepforest import main
    MODEL=main.deepforest()
    MODEL.load_model(model_name=HF_MODEL_PROVENANCE['model_name'],revision=HF_MODEL_PROVENANCE['resolved_revision'])
    changed=False
    for obj in (MODEL.model,getattr(MODEL.model,'model',None)):
        if obj is not None and hasattr(obj,'score_thresh'):obj.score_thresh=.05;changed=True
    for cfg in (MODEL.config,getattr(MODEL.config,'retinanet',None)):
        if cfg is None:continue
        try:previous=cfg['score_thresh'];cfg['score_thresh']=.05;changed=True
        except (TypeError,KeyError,AttributeError):pass
    if not changed:raise RuntimeError('NEON3차와같이내부score floor.05를설정하지못했습니다.임계값누락을숨기지않습니다.')
    json_save(OUT/'model_provenance.json',HF_MODEL_PROVENANCE)
    return MODEL

def estimate_boxes():
    frames=[]
    for _,p in PLOTS.iterrows():
        low,rgb,chm=RASTERS[p.plot_id]
        fp=dict(schema=SCHEMA,run_signature=INPUT_SIGNATURE,geometry=p.geometry.wkt,
                rgb_sha=sha_file(rgb),chm_sha=sha_file(chm),hf=HF_MODEL_PROVENANCE,score=.2,patch=600,min_height=3.)
        ident=hashlib.sha256(json.dumps(fp,sort_keys=True).encode()).hexdigest()[:16]
        directory=CACHE/'detections'/ident;directory.mkdir(parents=True,exist_ok=True)
        processed=directory/'detected.csv';stamp=directory/'complete.json'
        if cache_valid(stamp,fp,[processed]):
            frame=geo_read(processed);print(p.plot_id,'탐지 캐시 재사용');frames.append(frame);continue
        raw=directory/f'{p.plot_id}_raw_predictions.csv'
        raw_stamp=directory/'raw_complete.json'
        if raw.exists() and (not raw_stamp.exists() or json.loads(raw_stamp.read_text()).get('sha256')!=sha_file(raw)):
            raw.replace(raw.with_suffix('.incomplete.csv'))
        # The original NEON inference, center filter and native CHM maximum are reused.
        cfg=dict(score_threshold=.2,patch_size=600,patch_overlap=.25,nms_iou=.15,min_height_m=3.,
                 agb_default=dict(beta0=-2.2118,beta1=2.4133))
        NEON.configure(cfg,directory,CRS_WORK)
        plot=SimpleNamespace(plotID=p.plot_id,siteID='GEUMJEONGSAN',field_year='2009-2010',
                            geometry=box(p.center_x-60,p.center_y-60,p.center_x+60,p.center_y+60))
        detector=None if raw.is_file() else load_model()
        records=NEON.predict_crowns(detector,rgb,chm,plot)
        json_save(raw_stamp,dict(sha256=sha_file(raw)))
        rows=[]
        for row in records:
            attrs=forest_at(row['geometry'].centroid)
            d,biome=jucker_diameter(row['H_m'],row['CD_m'],attrs['forest_code'],JUCKER)
            value=ALLO.predict(attrs['map_species'],d,attrs['forest_code'],'map')
            r=dict(row,zone=p.zone,plot_id=p.plot_id,**attrs,DBH_cm=d,Jucker_biome=biome,**value)
            for delta,label in ((-.1,'height_minus10'),(.1,'height_plus10')):
                changed,_=jucker_diameter(row['H_m']*(1+delta),row['CD_m'],attrs['forest_code'],JUCKER)
                v=ALLO.predict(attrs['map_species'],changed,attrs['forest_code'],'map')
                r['carbon_'+label+'_kg']=v['carbon_kg'];r['DBH_'+label+'_cm']=changed
            rows.append(r)
        if rows:frame=gpd.GeoDataFrame(rows,geometry='geometry',crs=CRS_WORK)
        else:
            columns=['det_id','plotID','plot_id','zone','score','H_m','CD_m','DBH_cm','AGB_kg','carbon_kg','carbon_cf_kg',
                     'carbon_height_minus10_kg','carbon_height_plus10_kg','forest_type','forest_code','map_species','formula','mapping','range_flag','geometry']
            frame=gpd.GeoDataFrame(columns=columns,geometry='geometry',crs=CRS_WORK)
        geo_csv(frame,processed.with_suffix('.part.csv'));processed.with_suffix('.part.csv').replace(processed)
        json_save(stamp,fp);frames.append(frame)
        print(p.plot_id,'탐지:',len(frame))
    result=gpd.GeoDataFrame(pd.concat(frames,ignore_index=True),geometry='geometry',crs=CRS_WORK)
    geo_csv(result,OUT/'detected_trees.csv')
    return result

def metric(reference,estimate):
    a=np.asarray(reference,float);b=np.asarray(estimate,float);good=np.isfinite(a)&np.isfinite(b);a=a[good];b=b[good]
    if not len(a):return dict(n=0,reference_mean=np.nan,estimate_mean=np.nan,estimate_reference_ratio=np.nan,r=np.nan,RMSE=np.nan,MAE=np.nan,Bias=np.nan)
    delta=b-a
    return dict(n=len(a),reference_mean=float(a.mean()),estimate_mean=float(b.mean()),
                estimate_reference_ratio=float(b.mean()/a.mean()) if a.mean()>0 else np.nan,
                r=float(np.corrcoef(a,b)[0,1]) if len(a)>1 and np.std(a)>0 and np.std(b)>0 else np.nan,
                RMSE=float(np.sqrt(np.mean(delta**2))),MAE=float(np.mean(abs(delta))),Bias=float(delta.mean()))

def calibrate_zone_holdout(g):
    output=pd.DataFrame(index=g.index,columns=['C_calibrated_OOF_tC_ha','calibration_factor'],dtype=float)
    output['heldout_zone']='';records=[]
    assert g.zone.nunique()==5
    for fold,(train,test) in enumerate(GroupKFold(n_splits=5).split(g,groups=g.zone)):
        tr=g.iloc[train];te=g.iloc[test]
        assert set(tr.zone).isdisjoint(set(te.zone)) and te.zone.nunique()==1
        factor=training_ratio(tr.C_reference_tC_ha,tr.C_estimated_tC_ha)
        output.loc[te.index,'calibration_factor']=factor
        output.loc[te.index,'C_calibrated_OOF_tC_ha']=te.C_estimated_tC_ha*factor
        output.loc[te.index,'heldout_zone']=te.zone.iloc[0]
        records.append(dict(outer_fold=fold,train_zones=';'.join(sorted(set(tr.zone))),test_zone=te.zone.iloc[0],
                            train_plots=';'.join(tr.plot_id),test_plots=';'.join(te.plot_id),factor=factor,
                            n_train=len(tr),n_test=len(te),fit='sum(train_reference_density)/sum(train_estimated_density)'))
    return output,records

def compare_plots(detected,reference):
    rows=[]
    scenarios=[('square20','CF0.5','height0','carbon_kg'),('circle30','CF0.5','height0','carbon_kg'),
               ('square20','species_CF','height0','carbon_cf_kg'),('square20','CF0.5','height_minus10','carbon_height_minus10_kg'),
               ('square20','CF0.5','height_plus10','carbon_height_plus10_kg')]
    for _,p in PLOTS.iterrows():
        refs=reference[reference.plot_id==p.plot_id]
        boxes=detected[detected.plot_id==p.plot_id]
        for area_kind,cf,height,column in scenarios:
            geometry=box(p.center_x-10,p.center_y-10,p.center_x+10,p.center_y+10) if area_kind=='square20' else p.geometry.buffer(30,quad_segs=128)
            area=400. if area_kind=='square20' else math.pi*30**2
            chosen=boxes[boxes.geometry.centroid.map(geometry.covers).astype(bool)].copy()
            for threshold,include in ((10,'included_ge10'),(6,'included_ge6')):
                fr=refs[refs[include]];refcol='carbon_cf_kg' if cf=='species_CF' else 'carbon_kg'
                cref=float(fr[refcol].sum())/1000/.04
                cest=float(chosen[column].sum())/1000*10000/area
                visible=fr[fr.visible_by_T1];hidden=fr[~fr.visible_by_T1]
                center=forest_at(p.geometry)
                formula=' | '.join(sorted(set(chosen.formula)))
                r=dict(plot_id=p.plot_id,zone=p.zone,area_kind=area_kind,area_m2=area,reference_area_m2=400,
                       CF_scenario=cf,height_scenario=height,field_DBH_min_cm=threshold,
                       C_reference_tC_ha=cref,C_estimated_tC_ha=cest,estimate_reference_ratio=cest/cref if cref>0 else np.nan,
                       n_detected=len(chosen),n_detected_per400m2=len(chosen)*400/area,
                       n_field=len(fr),n_visible=len(visible),n_occluded=len(hidden),
                       C_visible_tC_ha=float(visible[refcol].sum())/1000/.04,C_occluded_tC_ha=float(hidden[refcol].sum())/1000/.04,
                       visible_carbon_share=float(visible[refcol].sum()/fr[refcol].sum()) if fr[refcol].sum()>0 else np.nan,
                       occluded_carbon_share=float(hidden[refcol].sum()/fr[refcol].sum()) if fr[refcol].sum()>0 else np.nan,
                       forest_type=center['forest_type'],map_species=center['map_species'],used_formula=formula,
                       detected_range_extrapolation_count=int(chosen.range_flag.fillna('').ne('').sum()),
                       detected_below10_count=int((chosen.DBH_cm<10).sum()),
                       C_estimated_DBH10only_tC_ha=float(chosen.loc[chosen.DBH_cm>=10,column].sum())/1000*10000/area,
                       reference_proxy_for_circle=area_kind=='circle30',T1_min_m=p.T1_min_m)
                rows.append(r)
    result=pd.DataFrame(rows)
    assert len(result)==50*5*2 and result.groupby(['area_kind','CF_scenario','height_scenario','field_DBH_min_cm']).size().eq(50).all()
    result['C_calibrated_OOF_tC_ha']=np.nan;result['calibration_factor']=np.nan;result['heldout_zone']=''
    folds=[]
    keys=['area_kind','CF_scenario','height_scenario','field_DBH_min_cm']
    for key,g in result.groupby(keys,sort=False):
        calibrated,records=calibrate_zone_holdout(g)
        for col in calibrated:result.loc[g.index,col]=calibrated[col]
        folds.extend(dict(zip(keys,key),**record) for record in records)
    result.to_csv(OUT/'plot_comparison_sensitivity.csv',index=False)
    primary=baseline_rows(result)
    primary['circle30_reference_is_20m_proxy']=True
    primary['C_reference_ge6_tC_ha']=primary.plot_id.map(result[(result.area_kind=='square20')&(result.CF_scenario=='CF0.5')&(result.height_scenario=='height0')&(result.field_DBH_min_cm==6)].set_index('plot_id').C_reference_tC_ha)
    for label,area_kind,cf,height in [('circle30','circle30','CF0.5','height0'),('speciesCF','square20','species_CF','height0'),
                                     ('height_minus10','square20','CF0.5','height_minus10'),('height_plus10','square20','CF0.5','height_plus10')]:
        scenario=result[(result.area_kind==area_kind)&(result.CF_scenario==cf)&(result.height_scenario==height)&(result.field_DBH_min_cm==10)].set_index('plot_id')
        primary['C_estimated_'+label+'_tC_ha']=primary.plot_id.map(scenario.C_estimated_tC_ha)
    primary.to_csv(OUT/'plot_comparison.csv',index=False)
    pd.DataFrame(folds).to_csv(OUT/'calibration_folds.csv',index=False)
    summaries=[]
    for key,g in result.groupby(keys,sort=False):
        for zone,q in [('ALL',g)]+list(g.groupby('zone')):
            a=metric(q.C_reference_tC_ha,q.C_estimated_tC_ha);b=metric(q.C_reference_tC_ha,q.C_calibrated_OOF_tC_ha)
            summaries.append(dict(zip(keys,key),zone=zone,source='금정산',**a,
                calibrated_RMSE=b['RMSE'],calibrated_MAE=b['MAE'],calibrated_ratio=b['estimate_reference_ratio'],calibrated_r=b['r'],
                n_detected=int(q.n_detected.sum()),n_visible=int(q.n_visible.sum()),n_occluded=int(q.n_occluded.sum()),F1_visible=np.nan))
    summaries.append(dict(source='NEON25cm 기준선(보정 미사용)',zone='ALL',area_kind='NEON original footprint',CF_scenario='CF0.5',
        height_scenario='height0',field_DBH_min_cm=10,n=23,estimate_reference_ratio=NEON_25CM_BASELINE['estimate_reference_ratio'],
        F1_visible=NEON_25CM_BASELINE['F1_visible'],RMSE=NEON_25CM_BASELINE['carbon_RMSE_tC_ha'],
        reference_mean=NEON_25CM_BASELINE['reference_mean'],estimate_mean=NEON_25CM_BASELINE['estimate_mean'],
        r=NEON_25CM_BASELINE['r'],MAE=NEON_25CM_BASELINE['MAE']))
    summary=pd.DataFrame(summaries);summary.to_csv(OUT/'summary_sensitivity.csv',index=False)
    core=summary[(summary.source!='금정산')|((summary.area_kind=='square20')&(summary.CF_scenario=='CF0.5')&(summary.height_scenario=='height0')&(summary.field_DBH_min_cm==10))]
    assert len(core)==7
    core.to_csv(OUT/'summary.csv',index=False)
    return result,summary

def baseline_rows(comparison):
    return comparison[(comparison.area_kind=='square20')&(comparison.CF_scenario=='CF0.5')&(comparison.height_scenario=='height0')&(comparison.field_DBH_min_cm==10)].copy()

def figures(detected,comparison):
    directory=OUT/'figures';baseline=baseline_rows(comparison)
    colors=dict(zip(('중앙','동측','서측','남측','북측'),('#386cb0','#f6a53a','#6ca448','#cf6689','#8b6bb2')))
    # Full image thumbnails only: no complete 3GB mosaic is read into memory.
    fig,ax=plt.subplots(figsize=(10,12))
    for o in ORTHO:
        signature=hashlib.sha256(json.dumps(stat_signature(o['path']),sort_keys=True).encode()).hexdigest()[:16]
        thumb=CACHE/'thumbnails'/(o['sheet']+'_'+signature+'.png');thumb.parent.mkdir(exist_ok=True)
        if thumb.exists():image=plt.imread(thumb)
        else:
            with rasterio.open(o['vrt']) as ds:
                ratio=max(ds.width,ds.height)/850
                arr=ds.read([1,2,3],out_shape=(3,int(ds.height/ratio),int(ds.width/ratio)),resampling=Resampling.nearest)
            image=arr.transpose(1,2,0);plt.imsave(thumb,image)
        left,bottom,right,top=o['geometry'].bounds
        ax.imshow(image,extent=(left,right,bottom,top),zorder=0)
        x,y=o['geometry'].exterior.xy;ax.plot(x,y,c='yellow',lw=.6);ax.annotate(o['sheet'],(left,top),fontsize=7,color='blue')
    for zone,g in PLOTS.groupby('zone'):ax.scatter(g.center_x,g.center_y,s=18,c=colors[zone],edgecolors='white',linewidths=.3,label=ZONE_EN[zone],zorder=4)
    ax.legend();ax.set_aspect('equal');ax.set_title('2021 orthophotos: nine manually referenced sheets / 50 field plots')
    ax.ticklabel_format(useOffset=False,style='plain');fig.savefig(directory/'overview_map.png',dpi=160,bbox_inches='tight');plt.close(fig)
    # Five zone representatives plus one additional opposite carbon extreme.
    low=baseline.sort_values('C_reference_tC_ha').iloc[0];high=baseline.sort_values('C_reference_tC_ha').iloc[-1]
    selected=[low.plot_id,high.plot_id]
    for zone,g in baseline.groupby('zone'):
        if zone not in set(baseline[baseline.plot_id.isin(selected)].zone):selected.append(g.sort_values('C_reference_tC_ha').iloc[len(g)//2].plot_id)
    for pid in baseline.sort_values('C_reference_tC_ha').plot_id:
        if len(selected)>=6:break
        if pid not in selected:selected.append(pid)
    assert len(selected)==6 and set(PLOTS[PLOTS.plot_id.isin(selected)].zone)==set(colors)
    pd.DataFrame({'plot_id':selected}).to_csv(OUT/'selected_overlay_plots.csv',index=False)
    for i,pid in enumerate(selected):
        p=PLOTS.set_index('plot_id').loc[pid];low_path=RASTERS[pid][0]
        fig,ax=plt.subplots(figsize=(8,8))
        with rasterio.open(low_path) as ds:show(ds,ax=ax)
        q=detected[detected.plot_id==pid]
        if len(q):q.boundary.plot(ax=ax,color='cyan',lw=.8)
        for geom,color,style in [(box(p.center_x-10,p.center_y-10,p.center_x+10,p.center_y+10),'yellow','-'),(p.geometry.buffer(30),'yellow','--')]:
            x,y=geom.exterior.xy;ax.plot(x,y,c=color,ls=style,lw=1.5)
        ax.scatter([p.center_x],[p.center_y],c='red',s=20,marker='+')
        ax.set_title(f'{ZONE_EN[p.zone]} {pid.split("-")[-1]} | cyan: crowns; yellow:20m square /30m circle')
        ax.ticklabel_format(useOffset=False,style='plain');fig.savefig(directory/f'{ZONE_EN[p.zone]}_{pid.split("-")[-1]}_overlay.png',dpi=160,bbox_inches='tight');plt.close(fig)
        if i<2:
            fig,ax=plt.subplots(figsize=(9,9))
            with rasterio.open(low_path) as ds:show(ds,ax=ax)
            extent=box(p.center_x-60,p.center_y-60,p.center_x+60,p.center_y+60)
            f=FOREST[FOREST.geometry.intersects(extent)];t=TRAILS[TRAILS.geometry.intersects(extent)]
            if len(f):f.boundary.plot(ax=ax,color='orange',lw=.8)
            if len(t):t.plot(ax=ax,color='cyan',lw=1.2)
            ax.set_xlim(p.center_x-60,p.center_x+60);ax.set_ylim(p.center_y-60,p.center_y+60)
            ax.set_title('Visual georeferencing check | orange:2025 forest map / cyan:trails')
            fig.text(.12,.02,'Reference bbox test is numerical; boundary alignment still requires visual review.',fontsize=9)
            fig.savefig(directory/f'georef_overlay_{i+1}.png',dpi=160,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for zone,g in baseline.groupby('zone'):
        for ax,column in zip(axes,('C_estimated_tC_ha','C_calibrated_OOF_tC_ha')):ax.scatter(g.C_reference_tC_ha,g[column],c=colors[zone],label=ZONE_EN[zone],s=24)
    for ax,column,title in zip(axes,('C_estimated_tC_ha','C_calibrated_OOF_tC_ha'),('Uncalibrated','Leave-one-zone-out calibrated')):
        limit=max(baseline.C_reference_tC_ha.max(),baseline[column].max())*1.06
        ax.plot([0,limit],[0,limit],'k--',lw=1);ax.set_xlim(0,limit);ax.set_ylim(0,limit);ax.set_aspect('equal')
        ax.set(xlabel='Field-derived aboveground carbon (tC/ha)',ylabel='Estimated aboveground carbon (tC/ha)',title=title);ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(directory/'plot_carbon_scatter.png',dpi=170,bbox_inches='tight');plt.close(fig)
    g=baseline.groupby('zone')[['C_visible_tC_ha','C_occluded_tC_ha','n_visible','n_occluded','n_detected']].mean().reindex(colors)
    fig,axes=plt.subplots(1,2,figsize=(13,5));xx=np.arange(5)
    axes[0].bar(xx,g.C_visible_tC_ha,label='Visible by T1 minimum',color='#73af75');axes[0].bar(xx,g.C_occluded_tC_ha,bottom=g.C_visible_tC_ha,label='Occluded by T1 minimum',color='#d982a4')
    for j,(col,label,c) in enumerate([('n_visible','Field visible','#73af75'),('n_occluded','Field occluded','#d982a4'),('n_detected','Detected','#5a86b6')]):axes[1].bar(xx+(j-1)*.25,g[col],.24,label=label,color=c)
    for ax in axes:ax.set_xticks(xx,[ZONE_EN[z] for z in g.index]);ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Mean field aboveground C (tC/ha)');axes[1].set_ylabel('Mean tree count per400m2')
    fig.tight_layout();fig.savefig(directory/'visible_occluded_bar.png',dpi=170,bbox_inches='tight');plt.close(fig)

def write_report(comparison,summary,detected,reference):
    base=summary[(summary.source=='금정산')&(summary.area_kind=='square20')&(summary.CF_scenario=='CF0.5')&(summary.height_scenario=='height0')&(summary.field_DBH_min_cm==10)]
    total=base[base.zone=='ALL'].iloc[0]
    proxy=detected.mapping.fillna('').str.contains('平均|평균',regex=True) if len(detected) else pd.Series([],dtype=bool)
    source_unknown=int(detected.get('forest_fallback',pd.Series([],dtype=bool)).fillna(False).astype(bool).sum())
    table=base[['zone','n','reference_mean','estimate_mean','estimate_reference_ratio','r','RMSE','MAE','calibrated_RMSE','calibrated_ratio']].to_markdown(index=False,floatfmt='.3f')
    sensitivity=summary[(summary.source=='금정산')&(summary.zone=='ALL')][['area_kind','CF_scenario','height_scenario','field_DBH_min_cm','reference_mean','estimate_mean','estimate_reference_ratio','RMSE','calibrated_RMSE']].to_markdown(index=False,floatfmt='.3f')
    report=f'''# 금정산 사진·높이 기반 지상부 탄소 비교

**50개20m 조사구, DBH≥10cm: 실측 기반 {total.reference_mean:.2f} / 영상 추정 {total.estimate_mean:.2f}tC/ha, 비율 {total.estimate_reference_ratio:.3f}, RMSE {total.RMSE:.2f}, zone OOF 보정 후 RMSE {total.calibrated_RMSE:.2f}tC/ha.**

## 전체·zone별 결과

{table}

모든 탄소값은 지상부(줄기+가지+잎)의 tC/ha이다. 뿌리·토양·고사목은 포함하지 않는다.
이전 금정산 A/B/C의 지상부+뿌리 탄소값과 계산 범위가 다르다.
논문 발표는 2011년이지만 실제 측정은 2009년 5월~2010년 9월이며, 정사영상은 2021년이다.

## NEON25cm 참고 기준선

같은 HF 커밋·score 0.2/patch 600·25cm→10cm 입력 실험의 추정/참조 비율은 {NEON_25CM_BASELINE['estimate_reference_ratio']:.6f},
F1_visible은 {NEON_25CM_BASELINE['F1_visible']:.6f}, 탄소 RMSE는 {NEON_25CM_BASELINE['carbon_RMSE_tC_ha']:.3f}tC/ha이다.
이 값은 보정계수로 사용하지 않았다. NEON은 지역·공식·참조 모집단이 달라 국내 성능을 보장하지 않는다.
금정산은 개별 줄기 좌표가 없어 precision/recall/F1이나 1:1 개체 매칭을 계산하지 않고, 탐지 수와 보이는 나무 수만 비교했다.

## 좌표와 영상

- 지정된 1:5000 도엽 규칙으로 경위도 모서리를 EPSG:5187로 변환했다. 좌상단은 (minx−50m,maxy+50m), 픽셀은 25cm이며 회전은 없다고 가정했다.
- 원본 TIFF 9개는 수정하지 않았다. 캐시 VRT에 좌표를 할당한 뒤 조사구 crop만 GeoTIFF로 저장했다.
- 35913003의 영상 전체를 EPSG:5179로 변환한 bbox가 기준의 각 변 3m 이내인지 assert했다. 9개 영상의 폭·높이 차이는 georef_check.csv에 기록했다.
- 각 중심 ±60m crop의 원본 25cm와 bilinear 10cm 확대본을 저장했다. 확대해도 선명도가 복원되지 않는다.
- Meta CHM은 대상 zoom9 quad만 사용하고 원격 COG의 필요한 window만 읽었다. 제품의 명목 1m 해상도와 native CRS 픽셀 간격을 구분하여 원본 tags를 보존했다.
  H는 native Meta 셀의 최대값으로 추출했다. 별도 EPSG:5187 1m CHM은 max 리샘플한 표시용이며 추론에는 쓰지 않았다.
- CHM 관측 메타데이터는 raster_manifest.csv와 meta_chm_sources.json에 있다. 날짜를 찾지 못하면 unknown으로 남기며, 파일 생성일과 촬영일을 구분했다.
- 임상도·등산로 겹침 그림 2개를 저장했다. 수치 bbox 검산은 현장 정합의 증명이 아니므로 사람이 그림을 확인해야 한다.

## 추정과 실측 공식

- HF모델:{HF_MODEL_PROVENANCE['model_name']} / 고정커밋:{HF_MODEL_PROVENANCE['resolved_revision']},DeepForest2.1.0.
- NEON의predict_crowns를재사용:score.2,patch600,overlap.25,NMS.15,상자중심영역판정,H≥3m,H=상자내native CHM최대,CD=상자가로·세로평균.
- Jucker:DBH(cm)=a×(H(m)×CD(m))^b×exp(g²/2).침엽biome21:{JUCKER[21]},활엽biome20:{JUCKER[20]}.
  혼효림은두DBH예측평균이다.수종미상은임상별대표(침엽:소나무류2식평균,활엽:낙엽활엽수평균,혼효:두그룹각0.5)를사용했다.
- 임상도밖/비산림·죽림등임상미확인탐지 {source_unknown}개는혼효대표대용으로계산하고forest_fallback으로표시했다.
  지도미피복을탄소0으로단정하거나탐지를숨겨삭제하지않았다.이대용의기여는mapping_contributions.csv.
- 국내2014상대생장식aD^b의줄기·가지·잎합을사용했다.그룹대용은계수평균이아닌각식의예측질량평균이다.
- 소나무D>40cm는강원식,나머지는중부식.떡갈/기타참나무는4종평균,개서어나무/소사나무는서어나무.
  야장의등나무와D<6cm는제외했다.필드기본은D≥10cm,보조는D≥6cm.범위밖은외삽+flag이다.
- 영상탐지는NEON과같이H≥3m만으로기본계산하며추정DBH<10cm도포함한다.추정DBH≥10cm만의탄소도별도열로저장해모집단차이를확인할수있다.
- CF기본.5,민감도는국내2010수종별목재CF를모든지상부부위에공통적용했다.미측정수종CF는.5로유지하고대용가중치를기록했다.
  대용그룹은각수종의AGB와그수종CF를먼저짝지어탄소를계산한뒤평균했다.CF열은그결과의유효CF이다.
- 교목층T1최소수고이상은위에서보이는나무로대리분류했다.이는실제수관가시성정답이아니며개체별겹침을알수없다.
- 원문표기의심개어서나무·대팻집·고뢰쇠나무는임의정정없이기타활엽수대용을적용했다.

## 보정과 민감도

zone 5개를 GroupKFold 5개로 나누어 매번 한 zone을 통째로 제외했다. 학습 4개 zone의 합계 실측 밀도/합계 추정 밀도로 보정계수를 구했다.
시험 zone의 정답은 계수에 사용하지 않았으며 train/test ID와 계수는 calibration_folds.csv에 보존했다.
학습 추정 합계가 0이면 계수를 구할 수 없어 NaN으로 표시한다. 이는 이 50개 자료 안의 공간 보정 시험이며 새 지역 성능 검증이 아니다.

{sensitivity}

30m 원은 추정 범위의 민감도다. 원 안의 나무를 새로 조사한 것이 아니므로 20m 조사구의 참조 밀도를 대용했다.
면적이 달라 탐지 수는 n_detected_per400m2도 함께 확인해야 한다. 높이 ±10%는 상자·탐지 H≥3m 필터·분석 영역을 고정하고 H→DBH→AGB만 변경했다.

## 해석상 한계

- 현장2009~2010과영상2021사이에11~12년성장·고사·교란이있다.논문연도2011과영상사이는약10년이다.CHM촬영일은별도메타데이터로확인해야한다.
- 초 단위 GPS의 약25~30m 오차는20m 조사구보다 크고, 조사구 회전 방향도 모른다. 축정렬20m 정사각형은 비교를 위한 근사이다. 원 GPS 측지계가 명시되지 않아 경위도는 WGS84(EPSG:4326)로 가정했다.
- 수치bbox1개일치는위치정확도인증이아니며,영상9개전체의회전·기하왜곡·50m외곽여유가완벽히같다는보장은없다.
- 실제국내25cm영상은인위적으로축소한NEON영상과센서·압축·촬영각·색감·계절이다르다.10cm확대는정보를복원하지않는다.
- MetaCHM자체모델오차·촬영시점·지형/수평정합오차가있고최대높이가이웃큰나무를반영할수있다.
- 임상도2025의임분수종을2021개체수종으로대용했다.혼효림구성비와개별수종은모르며경계/중복의영향이있다.
- Jucker구북구온대식과국내수종식이금정산/상자CD에맞는지는검증되지않았다.모든범위외삽과평균대용을표시했다.
- 비교대상의실측탄소역시상대생장식참조값으로실제무게측정값이아니다.수고에따른가시성대리분류와DBH10cm기준도불확실성을만든다.

## 검산과 출처

- plot_comparison.csv는기본50행,별도민감도long표는50개×5시나리오×2직경기준=500행이다.기본summary.csv는전체·5zone·NEON참고7행이다.
  zone별검증분리,400㎡/ha·kg/t·CF환산,수종평균식·DBH40cm분기검사를수행했다.
- [Jucker원문구현]({JUCKER_URL}),계수와SHA256은provenance.json.
- 국내2014/2010보고서의제공CSV와README를사용했다.나무별formula열에수종·부위·계수·가중치·원문쪽을기록했다.
- README의 가지 예시0.069×24^1.995는39.6kg가 아니라39.117448kg이다. 예시 산술 불일치를 기록하고 제공 CSV 계수는 그대로 사용했다.
- [Meta/WRI공개자료](https://registry.opendata.aws/dataforgood-fb-forests/),출처Maxar/Meta/WRI;관련관측메타데이터보존.
- 국토지리정보원정사영상XML2021,산림청임상도2025·등산로,서정범2011동아대박사부록야장(측정2009~2010).
- 계산·그림은Colab실행결과이며Windows검증은표준라이브러리와노트북구조확인에한정한다.학위논문자료는공개업로드하지않는다.
'''
    (OUT/'report.md').write_text(report,encoding='utf-8')
    pd.DataFrame(AUDIT,columns=['stage','item','message']).to_csv(OUT/'audit.csv',index=False)
    pd.DataFrame(NEON.AUDIT,columns=['stage','item','reason']).to_csv(OUT/'neon_inference_audit.csv',index=False)
    contributions=[]
    for label,df in [('detected_60m_crops',detected),('field_ge10',reference[reference.included_ge10])]:
        total=float(df.carbon_kg.sum())
        for rule,g in df.groupby('mapping'):
            contributions.append(dict(population=label,mapping=rule,n=len(g),carbon_kg=float(g.carbon_kg.sum()),carbon_share_pct=100*float(g.carbon_kg.sum())/total if total else np.nan,
                note='detected_60m includes per-plot buffers and possible repeated real trees; not a landscape total' if label.startswith('detected') else '50 original20m field plots'))
    pd.DataFrame(contributions).to_csv(OUT/'mapping_contributions.csv',index=False)
    json_save(OUT/'run_validation.json',dict(plots=50,reference_rows=len(reference),comparison_rows=len(comparison),georef_control_passed=True,
        no_root_biomass=True,calibration='leave-one-zone-out; train-only ratio',no_individual_F1=True,
        synthetic_and_unit_checks='passed before real inference',georef_visual_review='two figures saved; requires human review'))

def runtime_checks():
    dimensional_checks()
    test={20:dict(a=.694,b=.73,g=0),21:dict(a=1.004,b=.73,g=0)}
    broad,_=jucker_diameter(20,5,'2',test);con,_=jucker_diameter(20,5,'1',test);mixed,_=jucker_diameter(20,5,'3',test)
    assert np.isclose(mixed,(broad+con)/2)
    # Real Rasterio crop/grid test, tiny data only; no source image mutation.
    with tempfile.TemporaryDirectory() as folder:
        file=Path(folder)/'test.tif';arr=np.arange(12*12,dtype='uint8').reshape(12,12)
        tr=from_origin(200000,300000,.25,.25)
        with rasterio.open(file,'w',driver='GTiff',height=12,width=12,count=3,dtype='uint8',crs=CRS_WORK,transform=tr) as ds:ds.write(np.stack([arr]*3))
        with rasterio.open(file) as ds:
            up=ds.read(out_shape=(3,30,30),resampling=Resampling.bilinear)
            transformed=ds.transform*rasterio.Affine.scale(.4,.4)
            assert up.shape==(3,30,30) and np.allclose((transformed.a,abs(transformed.e)),.1)
            assert np.allclose(rasterio.transform.array_bounds(12,12,tr),rasterio.transform.array_bounds(30,30,transformed))
    assert ALLO.mapping('소나무',40)[0][0][0]=='중부지방소나무'
    assert ALLO.mapping('소나무',40.1)[0][0][0]=='강원지방소나무'
    assert np.isclose(ALLO.predict('떡갈나무',20)['AGB_kg'],np.mean([ALLO.predict(s,20)['AGB_kg'] for s in ALLO.oaks]))
    # Perturb held-out targets and rerun the actual calibration implementation.
    synthetic=pd.DataFrame(dict(plot_id=list('ABCDE'),zone=list('ABCDE'),C_reference_tC_ha=[10.,20,30,40,50],C_estimated_tC_ha=[5.,10,15,20,25]))
    before,_=calibrate_zone_holdout(synthetic)
    changed=synthetic.copy();changed.loc[changed.zone=='A','C_reference_tC_ha']=1e9
    after,_=calibrate_zone_holdout(changed)
    assert before.loc[0,'calibration_factor']==after.loc[0,'calibration_factor']==2
    assert before.loc[0,'C_calibrated_OOF_tC_ha']==after.loc[0,'C_calibrated_OOF_tC_ha']==10
    print('Colab 실제 라이브러리 격자/단위/수종/혼효림 예측량평균 검증 통과')

ZONE_EN={'중앙':'Central','동측':'East','서측':'West','남측':'South','북측':'North'}
