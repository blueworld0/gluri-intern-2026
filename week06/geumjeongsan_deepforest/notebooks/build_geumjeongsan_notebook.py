"""Author a self-contained Colab notebook without importing scientific packages."""
import ast,hashlib,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]
NEON=WORK/'neon_grsm_ornl/notebooks/neon_pipeline.py'
MODEL=WORK/'neon_grsm_ornl/outputs_v4_resolution/resolution_model_provenance.json'
BASELINE=WORK/'neon_grsm_ornl/outputs_v4_resolution/resolution_comparison.csv'
import csv
with BASELINE.open(encoding='utf-8-sig',newline='') as stream:
    baseline=next(r for r in csv.DictReader(stream) if r['resolution_cm']=='25' and r['setting']=='fixed' and r['siteID']=='ALL')
baseline={new:float(baseline[old]) for old,new in (
    ('estimate_reference_ratio','estimate_reference_ratio'),('F1_visible','F1_visible'),('carbon_RMSE_tC_ha','carbon_RMSE_tC_ha'),
    ('C_ref_mean_tC_ha','reference_mean'),('C_est_mean_tC_ha','estimate_mean'),('carbon_r','r'),('carbon_MAE_tC_ha','MAE'))}
core=(HERE/'geumjeongsan_core.py').read_text(encoding='utf-8')
runtime=(HERE/'geumjeongsan_runtime.py').read_text(encoding='utf-8')
neon=NEON.read_text(encoding='utf-8');model=json.loads(MODEL.read_text(encoding='utf-8'))
fingerprint=hashlib.sha256((core+runtime+neon).encode()).hexdigest()
cells=[]
def md(text):cells.append(dict(cell_type='markdown',metadata={},source=text.strip().splitlines(True)))
def code(text):cells.append(dict(cell_type='code',metadata={},execution_count=None,outputs=[],source=text.strip().splitlines(True)))
def section(title,description,src):md('## '+title+'\n\n'+description);code(src)

md('''# 금정산 2021 정사영상·Meta CHM → 국내 공식 지상부 탄소

**Google Colab / T4 GPU용 새 노트북입니다. Windows에서 과학 패키지를 실행하지 않습니다.**
서정범 야장50개에서 사진·높이 추정과 국내 상대생장식 실측 기반 값을 비교합니다.
기본 참조는DBH≥10cm, 보조는≥6cm입니다. **줄기+가지+잎 지상부 탄소만 계산하며 뿌리는 제외합니다.**
개별 나무 위치가 없어F1/precision/recall·헝가리안 개체매칭을 수행하지 않습니다.
NEON의고정 탐지·CHM최대·상자CD 함수만재사용하며 한국수종공식으로교체합니다.

실행 순서: 아래 업로드 안내 → 패키지 설치 → 런타임 재시작이 필요하면 재시작 → Drive마운트부터 위에서 아래로 실행.
첫 실행은50개crop·COG접속·GPU·Drive속도에 따라 수십분 이상 걸릴 수 있습니다.
원본9장약3GB는수정하지않으며 crop과탐지를Drive에캐시합니다. 공개저장소로자료를전송하는코드는없습니다.
''')
md('''## 00. Drive 업로드 구조

이노트북하나만Colab에올리면됩니다.아래입력자료는회사승인범위내개인/조직Drive에직접업로드하세요.
Drive루트기본값은 `/content/drive/MyDrive/geumjeongsan`입니다.

```text
geumjeongsan/
  gedi_geumjeongsan/data/
    ortho/2021/       # 지정9개 (B060)…tif + 대응XML9개, 이름유지
    field/           # seo2011_plots.csv, seo2011_trees.csv, seo2011_layers.csv
    reference/       # README.md, nifos2014_allometry.csv,
                     # nifos2014_carbon_factors.csv, nifos2010_carbon_fraction.csv
    forest_map/
      2025_임상도/{금정구,북구,양산시}/  # .shp/.shx/.dbf/.prj/.cpg 등 함께
      등산로/{금정구,북구,양산시}/       # 동일
    raw/meta_chm/    # 이미받은native CHM quadkey TIFF가있으면선택업로드
```

자료를 `geumjeongsan/data/`에바로배치한경우도자동인식합니다. PDF원문이나임목축적KOSIS파일은필요없습니다.
NEON커밋·25cm기준선·재사용함수는노트북에내장되어추가모듈업로드가필요없습니다.
Meta CHM은이미받은quad파일을먼저찾고없으면해당quad의필요COG window만원격으로읽습니다.
HF모델과Jucker계수원문은공개서버에서받습니다.국내정사영상을공개서버로전송하지않습니다.
''')
section('01. Colab 패키지 설치','설치후재시작안내가나오면재시작하고02번부터실행하세요. NEON과같은DeepForest2.1.0을사용합니다.',
'''%pip -q install "deepforest==2.1.0" "scikit-learn==1.6.1" rasterio geopandas shapely pyproj numpy scipy scikit-image matplotlib requests tabulate''')
section('02. Drive 마운트·설정','회사승인입력만사용하세요.이노트북의출력과캐시는아래ROOT하위에저장됩니다. HF캐시는Drive와분리합니다.',f'''
from google.colab import drive
drive.mount('/content/drive')
from pathlib import Path
import os
ROOT=Path('/content/drive/MyDrive/geumjeongsan')
ROOT.mkdir(parents=True,exist_ok=True)
EXTRA_CHM_DIRS=[
    '/content/drive/MyDrive/gedi_geumjeongsan/data/raw/meta_chm',
    '/content/drive/MyDrive/gedi_hongcheon/data/raw/chm',
]
os.environ['HF_HOME']='/content/hf_cache'
os.environ['HF_HUB_CACHE']='/content/hf_cache/hub'
os.environ['TORCH_HOME']='/content/torch_cache'
Path('/content/hf_cache').mkdir(parents=True,exist_ok=True)
HF_MODEL_PROVENANCE={model!r}
NEON_25CM_BASELINE={baseline!r}
GEUM_CODE_SHA={fingerprint!r}
''')
section('03. 순수계산 함수','도엽·quadkey·국내수종대응·부위별공식입니다.실측과영상추정이같은국내공식을사용하며평균대용은예측질량을평균합니다.',core)
section('04. Colab 처리 함수','기존NEON모듈을내장합니다.새자료다운로드·3차결과재계산셀은없으며추론·수고·CD부분만호출합니다.',f'NEON_PIPELINE_SOURCE={neon!r}\n'+runtime)
section('05. 입력 확인·README·공식 원문 검증','필수파일과50개조사구를확인하고README를출력합니다.Jucker원문biome20/21을파싱하고계수·해시를기록합니다.지정숫자와다르면원문을사용하며경고합니다.',
'''initialize()
runtime_checks()
print('사용할 Jucker계수:',{k:JUCKER[k] for k in (20,21)})
display(PLOTS.drop(columns='geometry').head())''')
section('06. 도엽 좌표 할당·엄격한 범위 검증','원본TIFF를수정하지않고캐시VRT에5187좌표를할당합니다.35913003기준bbox각변3m이내와9개영상폭·높이차이를검사합니다.실패하면GPU추론전중단합니다.',
'''GEOREF_CHECK=georeference()
print('35913003 기준범위 대조:')
display(GEOREF_CHECK[GEOREF_CHECK.sheet=='35913003'])''')
section('07. 임상도·등산로 필드 확인','임상/수종필드의실제이름과예시를먼저출력합니다.코드필드만있으면명시한코드표를사용합니다.지도밖·비입목등은대용여부를표시합니다.',
'''read_vectors()
print('임상도/등산로:',len(FOREST),len(TRAILS))''')
section('08. 야장 참조 계산','기본≥10cm,보조≥6cm.줄기+가지+잎합×CF.5이며뿌리는제외합니다.T1최소수고이상/미만의개수와탄소를나눌수있도록보임표시를저장합니다.',
'''FIELD_REFERENCE=field_reference()
display(FIELD_REFERENCE[['plot_id','수종','DBH_cm','H_m','included_ge10','visible_by_T1','AGB_kg','carbon_kg','range_flag']].head())
print('원자료/6cm이상/10cm이상:',len(FIELD_REFERENCE),int(FIELD_REFERENCE.included_ge6.sum()),int(FIELD_REFERENCE.included_ge10.sum()))''')
section('09. 원본25cm·10cm확대본·CHM crop 캐시','각중심±60m,120m폭crop입니다.25cm사진은bilinear10cm로확대합니다.Meta CHMnative crop은H최대추출용,5187 1m CHM은별도표시용입니다.관측일을못찾으면unknown으로남깁니다.',
'''build_crops()
display(CHM_INFO[['plot_id','observation_dates','valid_fraction']])''')
section('10. 고정DeepForest 추론과국내공식','같은HF커밋,score.2/patch600,CHM≥3m.H·CD는NEON함수그대로,DBH는구북구20/21,AGB는국내2014지상부식으로교체합니다.원시탐지·처리결과를Drive에캐시합니다.',
'''DETECTED=estimate_boxes()
print('60m crop 안 탐지(비교영역필터전):',len(DETECTED))
display(DETECTED.drop(columns='geometry').head())''')
section('11. 조사구 비교·zone분리 보정·민감도','공식20m정사각형과30m원,CF,높이±10%,참조10/6cm를계산합니다.30m원의실측밀도는20m조사구대용입니다.보정계수는시험zone을뺀4개zone에서만적합합니다.NEON기준선은보정에쓰지않습니다.',
'''PLOT_COMPARISON,SUMMARY=compare_plots(DETECTED,FIELD_REFERENCE)
display(baseline_rows(PLOT_COMPARISON)[['plot_id','zone','C_reference_tC_ha','C_estimated_tC_ha','estimate_reference_ratio','C_calibrated_OOF_tC_ha']])''')
section('12. 지도·시각 정합·비교 그림','도엽9개와조사구50개overview,임상도·등산로겹침2장,zone별산점도,가시/가림막대,5zone를포함한저/고탄소6개overlay를저장합니다.실측개체줄기점은없으므로임의로나무위치를그리지않습니다.',
'''figures(DETECTED,PLOT_COMPARISON)
from IPython.display import display,Image
display(Image(filename=str(OUT/'figures/overview_map.png')))
for path in sorted((OUT/'figures').glob('georef_overlay_*.png')):
    display(Image(filename=str(path)))
print('임상도·등산로 정합은 위그림을사람이확인해야합니다.수치assert만으로정합을확정하지않습니다.')''')
section('13. 한국어 보고서·출처·검산 기록','맨위한줄요약에기본20m·≥10cm의결과를씁니다.좌표오차·시점·임상대용·원공식범위·NEON25cm차이를기록합니다.',
'''write_report(PLOT_COMPARISON,SUMMARY,DETECTED,FIELD_REFERENCE)
print('보고서:',OUT/'report.md')
print('생성파일:')
for p in sorted(OUT.rglob('*')):
    if p.is_file():print(p.relative_to(OUT))''')
section('14. 사람이 읽는 최종 요약','기본비교는400㎡·실측DBH≥10cm·CF.5·높이원값입니다.모든탄소와RMSE단위는tC/ha이며NEON기준선은참고행입니다.',
'''human=SUMMARY[((SUMMARY.source=='금정산')&(SUMMARY.area_kind=='square20')&(SUMMARY.CF_scenario=='CF0.5')&(SUMMARY.height_scenario=='height0')&(SUMMARY.field_DBH_min_cm==10))|(SUMMARY.source!='금정산')].copy()
human=human[['source','zone','n','reference_mean','estimate_mean','estimate_reference_ratio','r','RMSE','calibrated_RMSE','F1_visible']]
human=human.rename(columns={'source':'자료','zone':'구역','n':'조사구수','reference_mean':'실측기반평균(tC/ha)','estimate_mean':'추정평균(tC/ha)','estimate_reference_ratio':'추정/실측','r':'상관r','RMSE':'RMSE(tC/ha)','calibrated_RMSE':'보정후RMSE(tC/ha)','F1_visible':'NEON_F1참고'})
display(human.round(3))
print('NEON25cm 참고비율:',NEON_25CM_BASELINE['estimate_reference_ratio'],'(보정 미사용)')''')
section('15. 결과ZIP·Drive사본·다운로드','결과폴더만ZIP으로묶습니다.수GB원본정사영상과crop/탐지캐시는Drive에보존하고ZIP에중복포함하지않습니다.',
'''from google.colab import files
archive=shutil.make_archive('/content/geumjeongsan_results','zip',root_dir=str(OUT))
shutil.copy2(archive,ROOT/'outputs/geumjeongsan_results.zip')
print('ZIP:',archive,'Drive 사본:',ROOT/'outputs/geumjeongsan_results.zip')
files.download('/content/geumjeongsan_results.zip')''')

for i,c in enumerate(cells):
    c['id']=f'geum-{i:03}'
    if c['cell_type']=='code':
        assert cells[i-1]['cell_type']=='markdown'
        src=''.join(c['source'])
        if not src.startswith('%'):ast.parse(src,filename=f'geum_cell_{i}')
notebook=dict(nbformat=4,nbformat_minor=5,metadata=dict(kernelspec=dict(display_name='Python3',language='python',name='python3'),language_info=dict(name='python',version='3.11'),colab=dict(name='GEUMJEONGSAN_DeepForest_2021.ipynb',provenance=[]),accelerator='GPU'),cells=cells)
target=HERE/'GEUMJEONGSAN_DeepForest_2021.ipynb'
target.write_text(json.dumps(notebook,ensure_ascii=False,indent=1),encoding='utf-8')
print(target,'cells=',len(cells),'syntax OK')
