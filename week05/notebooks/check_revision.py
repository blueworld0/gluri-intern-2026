"""Regression/metadata checks using stdlib only; never import GIS/ML packages."""
import ast
import collections
import csv
import hashlib
import json
import math
import sqlite3
import statistics
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
CHECK=BASE/'outputs/revision_check'
CHECK.mkdir(exist_ok=True)
nb=json.loads((HERE/'NEON_DeepForest_AGB_2021.ipynb').read_text(encoding='utf-8'))
functions={}
for cell in nb['cells']:
    text=''.join(cell['source'])
    if cell['cell_type']=='code' and not text.startswith('%'):
        for node in ast.parse(text).body:
            if isinstance(node,ast.FunctionDef): functions[node.name]=node
namespace={'math':math,'Path':Path}
for name in ['product_site_entry','select_aop_pair','plot_biomass_totals','aop_tiff_files']:
    exec(compile(ast.Module(body=[functions[name]],type_ignores=[]),name,'exec'),namespace)
select=namespace['select_aop_pair']; totals=namespace['plot_biomass_totals']; entry=namespace['product_site_entry']

# Real product API evidence: ORNL 2021 exists but is provisional. Never replace it with 2022.
rgb=json.loads((HERE/'api_evidence/rgb_product.json').read_text(encoding='utf-8'))
chm=json.loads((HERE/'api_evidence/chm_product.json').read_text(encoding='utf-8'))
selected=select(entry(rgb,'ORNL'),entry(chm,'ORNL'),2021,True,[2022,2020],True)
assert selected['imagery_year']==2021 and selected['include_provisional']
assert selected['common_months']==['2021-06']
try: select(entry(rgb,'ORNL'),entry(chm,'ORNL'),2021,True,[2022,2020],False)
except ValueError: pass
else: raise AssertionError('PROVISIONAL disabled must not silently switch years')
synthetic={'availableMonths':['2020-06','2022-06'],'availableReleases':[]}
assert select(synthetic,synthetic,2021,True,[2022,2020])['imagery_year']==2022
try: select(synthetic,synthetic,2021,False)
except ValueError: pass
else: raise AssertionError('GRSM/nonfallback site must not change imagery year')

# Missing coordinates are valid for total reference. Missing biomass is not.
truth=[dict(individualID='A',ref_AGB_kg=100),dict(individualID='B',ref_AGB_kg=200)]
det=[dict(det_id='D',AGB_default_kg=80)]
pairs=[dict(det_id='D',individualID='A')]
plot=dict(plotID='T',siteID='TEST',field_year=2021,imagery_year=2022,
          area_m2=390,reported_area_m2=400,polygon_area_m2=390)
t=totals(det,truth,truth[:1],pairs,plot)
assert math.isclose(t['AGB_ref_t_ha'],7.5) and math.isclose(t['AGB_ref_coordonly_t_ha'],2.5)
assert t['area_m2']==400 and t['polygon_area_m2']==390
assert not t['coordinates_complete'] and t['reference_complete'] and not t['same_year']
assert t['n_without_matching_coordinate']==1
assert math.isclose(t['unmapped_truth_t_ha'],5)
assert math.isclose(t['AGB_est_t_ha']-t['AGB_ref_t_ha'],t['matched_error_t_ha']+t['unmatched_detection_t_ha']-t['missed_truth_t_ha'])
missing=totals(det,[truth[0],dict(individualID='B',ref_AGB_kg=None)],truth[:1],pairs,plot)
assert math.isnan(missing['AGB_ref_t_ha']) and missing['n_missing_ref_AGB']==1
no_coordinates=totals(det,truth,[],[],plot)
assert no_coordinates['reference_complete'] and math.isfinite(no_coordinates['AGB_ref_t_ha'])

baseline=json.loads((HERE/'preserved_pipeline_baseline.json').read_text())
for name,digest in baseline.items():
    if name=='sensitivity_cell':
        source=next(''.join(c['source']) for c in nb['cells'] if ''.join(c['source']).startswith('sensitivity=[]'))
        actual=hashlib.sha256(source.encode()).hexdigest()
    else: actual=hashlib.sha256(ast.dump(functions[name]).encode()).hexdigest()
    assert actual==digest,f'Protected pipeline changed: {name}'

def read_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def iter_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f: yield from csv.DictReader(f)
old=BASE/'outputs'
plots=read_csv(old/'plot_comparison.csv')
mapped=read_csv(old/'field_trees.csv')
detected=read_csv(old/'detected_trees.csv')
matches=read_csv(old/'matches_all.csv')
plotids={r['plotID'] for r in plots}
# Read only nongeometry GPKG attributes with sqlite3; no geopandas/rasterio on Windows.
con=sqlite3.connect(f'file:{(old/"sampled_footprints.gpkg").as_posix()}?mode=ro',uri=True)
table=con.execute("SELECT table_name FROM gpkg_contents WHERE data_type='features'").fetchone()[0]
safe='"'+table.replace('"','""')+'"'
areas=dict(con.execute(f'SELECT plotID,reported_area_m2 FROM {safe}'))
con.close()

# Reconstruct historical full reference population from bundled fixed AGBv2 for a cross-check.
# Actual rerun uses field_all from DP1, not this reconstruction.
best={}
for row in iter_csv(BASE/'data/raw/neonforestagb/NEONForestAGBv2.csv'):
    if row['plotID'] not in plotids or not row['date'].startswith('2021'): continue
    if not row['plantStatus'].startswith('Live') or 'tree' not in row['growthForm'].lower(): continue
    try:
        if float(row['stemDiameter'])<10 or not math.isfinite(float(row['AGB'])) or float(row['AGB'])<=0: continue
    except ValueError: continue
    rank={'AGBChojnacky':0,'AGBJenkins':1}.get(row['allometry'])
    if rank is None: continue
    key=row['individualID']; previous=best.get(key)
    if previous is None or rank<previous[0] or (rank==previous[0] and row['date']>previous[1]['date']):
        row['ref_AGB_kg']=float(row['AGB']); best[key]=(rank,row)
all_truth=[x[1] for x in best.values()]
rows=[]; mapped_est=[]; mapped_ref=[]
for oldp in plots:
    pid=oldp['plotID']; p=dict(oldp)
    p.update(field_year=2021,imagery_year=2021,reported_area_m2=areas[pid],polygon_area_m2=float(oldp['area_m2']))
    fa=[r for r in all_truth if r['plotID']==pid]
    fc=[r for r in mapped if r['plotID']==pid]
    ds=[r for r in detected if r['plotID']==pid]
    ms=[r for r in matches if r['plotID']==pid]
    if not {r['individualID'] for r in fc}<={r['individualID'] for r in fa}:
        raise AssertionError(f'{pid}: historical AGBv2 cannot reconstruct every mapped reference')
    result=totals(ds,fa,fc,ms,p)
    result['old_expected_field_count']=int(oldp['n_field_before_coordinate_filter'])
    result['historical_reconstruction_count_matches']=len(fa)==result['old_expected_field_count']
    rows.append(result)
    # Reproduce user's coordinates-only check with old polygon denominators before area fix.
    mapped_est.append(float(oldp['C_est_tC_ha']))
    mapped_ref.append(sum(float(r['ref_AGB_kg']) for r in fc)/1000*10000/float(oldp['area_m2'])*0.5)

def corr(a,b):
    return statistics.correlation(a,b) if len(a)>1 and len(set(a))>1 and len(set(b))>1 else float('nan')
def stats(a,b):
    errors=[p-r for r,p in zip(a,b)]
    return dict(n=len(a),MAE=statistics.mean(abs(e) for e in errors),
                RMSE=math.sqrt(statistics.mean(e*e for e in errors)),Bias=statistics.mean(errors),
                r=corr(a,b),mean_estimate=statistics.mean(b),mean_reference=statistics.mean(a),
                estimate_reference_ratio=statistics.mean(b)/statistics.mean(a))
before=stats(mapped_ref,mapped_est)
after=stats([r['C_ref_tC_ha'] for r in rows],[r['C_est_tC_ha'] for r in rows])
assert 75<before['mean_estimate']<81 and 130<before['mean_reference']<145
assert 0.5<before['estimate_reference_ratio']<0.65 and 0.35<before['r']<0.5
assert after['mean_reference']>statistics.mean(r['C_ref_coordonly_tC_ha'] for r in rows)
with (CHECK/'GRSM_plot_totals_crosscheck.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary=dict(ORNL_api_selection=selected,coordinates_only_old_area=before,
             full_reference_official_area=after,
             matched_reconstructed_counts=sum(r['historical_reconstruction_count_matches'] for r in rows),
             n_plots=len(rows),pipeline_unchanged=list(baseline))
(CHECK/'crosscheck.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
text=f'''# 1차 수정 검산

이 검산은 기존 GRSM 탐지를 다시 실행하지 않고 CSV·공식 API metadata·GPKG의 면적 속성과
고정 NEONForestAGBv2 파일만 표준 라이브러리로 읽어 수행했습니다. 원본 outputs/ 파일은 변경하지 않았습니다.

- ORNL: RGB·CHM 공통 2021-06 제공, PROVISIONAL=True. 2022로 바꾸지 않음.
- 좌표 누락이 있어도 총량 참조값 계산: 회귀 검사 통과.
- 실제 AGB 누락이면 총량 NaN: 회귀 검사 통과.
- 공식 조사면적 우선, 전체/좌표 한정 참조량 분리, 매칭 누락 분해 항등식: 통과.
- 탐지·Jucker·Chojnacky·개체 매칭·회귀·수고 민감도 코드 AST 보존: 통과.
- 전체 참조 모집단을 AGBv2로 복원한 개체 수가 이전 DP1 필터 전 개체 수와 일치: {summary['matched_reconstructed_counts']}/{len(rows)} 조사구.

| 검산 범위 | 평균 추정 tC/ha | 평균 참조 tC/ha | 추정/참조 | r |
|---|---:|---:|---:|---:|
| 기존 좌표 있는 나무만, 기존 폴리곤 면적 | {before['mean_estimate']:.3f} | {before['mean_reference']:.3f} | {before['estimate_reference_ratio']:.3f} | {before['r']:.3f} |
| 좌표 없는 나무 포함, 공식 조사면적 우선 | {after['mean_estimate']:.3f} | {after['mean_reference']:.3f} | {after['estimate_reference_ratio']:.3f} | {after['r']:.3f} |

두 번째 행은 과거 고정 AGBv2로 전체 참조를 재구성한 교차 검산입니다.
최종 노트북은 최신 DP1의 field_all을 기준으로 조인하므로 배포본 개정에 따라 값이 달라질 수 있습니다.
ORNL 실제 영상 다운로드와 GPU 탐지는 Colab 재실행이 필요합니다.
'''
(CHECK/'README.md').write_text(text,encoding='utf-8')
print(json.dumps(summary,indent=2))
print('PASS: revision regression checks; no pandas/geopandas/rasterio imports')
