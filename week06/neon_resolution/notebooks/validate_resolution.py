"""Scientific-package-free tests of the resolution extension and real v3 caches."""
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import statistics
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / 'outputs_v3'

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def main():
    nb = json.loads((HERE/'NEON_DeepForest_AGB_2021.ipynb').read_text(encoding='utf-8'))
    sources = [''.join(c['source']) for c in nb['cells'] if c['cell_type']=='code']
    combined = '\n'.join(s for s in sources if not s.startswith('%'))
    tree = ast.parse(combined)
    functions = {n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
    protected = ast.parse((HERE/'neon_pipeline.py').read_text(encoding='utf-8'))
    names = ('estimate_dbh','estimate_agb','predict_crowns','match_one_to_one',
             'nested_group_parameter_selection','evaluate_detection_method','plot_biomass_totals')
    for n in protected.body:
        if isinstance(n,ast.FunctionDef) and n.name in names:
            assert ast.dump(n)==ast.dump(functions[n.name]), n.name
    # Execute definition/guard cells with all scientific imports removed.
    # Any old inference/download action escaping its guard raises NameError.
    guarded = 0
    env = {'RESOLUTION_ONLY':True,'CFG':{'agb_default':{'beta0':-2.2118,'beta1':2.4133}}}
    for source in sources:
        if source.startswith('%') or 'if not RESOLUTION_ONLY:' not in source or 'drive.mount' in source:
            continue
        module = ast.parse(source)
        # The final summary intentionally runs; skip it for this guard-only test.
        if 'resolution_summary_display' in source:
            continue
        module.body = [n for n in module.body if not isinstance(n,(ast.Import,ast.ImportFrom))]
        exec(compile(module,'guard_test','exec'),env)
        guarded += 1
    assert guarded == 16
    assert 'RESOLUTION_ONLY = True' in combined and "TOKEN = ''" in combined
    assert sources[0].startswith('%pip') and 'deepforest==2.1.0' in sources[0]
    assert 'scikit-learn==1.6.1' in sources[0]
    assert "files.download('/content/neon_results.zip')" in ''.join(nb['cells'][-1]['source'])
    assert 'from rasterio.enums import Resampling' in combined
    assert 'resampling=Resampling.average' in combined and 'resampling=Resampling.bilinear' in combined
    assert 'model_input.transform == transform0' in combined
    assert 'model_input.shape == shape0' in combined
    assert 'no_imagery_download=True' in combined
    assert 'old_methods_recomputed=False' in combined
    # No downloader is called from the resolution runtime itself.
    runtime = ast.parse((HERE/'resolution_experiment.py').read_text(encoding='utf-8'))
    forbidden = {'cached_download','by_tile_aop','crop_mosaic','load_field_year','label_visibility','chm_watershed_crowns'}
    assert not any(isinstance(n,ast.Call) and (getattr(n.func,'id','') in forbidden or getattr(n.func,'attr','')=='by_tile_aop') for n in ast.walk(runtime))
    # Pure helpers are actually executed, not merely searched for text.
    pure = {'math':math,'RESOLUTION_THRESHOLDS':(.1,.2,.3,.4,.5),'Path':Path,'hashlib':hashlib}
    for name in ('resolution_grid','resolution_verify_cached_winners','resolution_sha'):
        exec(compile(ast.Module(body=[functions[name]],type_ignores=[]),name,'exec'),pure)
    grid = pure['resolution_grid']
    assert grid(1000,1000,.1,.1,.25)==(400,400)
    assert grid(1000,1000,.1,.1,.51)==(196,196)
    for width,height in ((1000,1000),(1001,999),(1,1)):
        for target in (.25,.51):
            w,h=grid(width,height,.1,.1,target)
            assert w>=1 and h>=1
            assert math.isclose(w*(width*.1/w),width*.1)
            assert math.isclose(h*(height*.1/h),height*.1)
    outer = read(BASE/'A_outer_groupkfold.csv')
    inner = read(BASE/'A_inner_grid_scores.csv')
    verified = pure['resolution_verify_cached_winners'](outer,inner)
    assert verified == {i:.2 for i in range(5)}
    wrong = [r.copy() for r in outer]; wrong[0]['patch_size']='400'
    try: pure['resolution_verify_cached_winners'](wrong,inner)
    except AssertionError: pass
    else: raise AssertionError('Must reject relabeling mixed-patch cached results')
    assignments = read(BASE/'A_fold_plot_assignments.csv')
    for fold in range(5):
        train={r['plotID'] for r in assignments if int(r['outer_fold'])==fold and r['role']=='train'}
        test={r['plotID'] for r in assignments if int(r['outer_fold'])==fold and r['role']=='test'}
        assert train.isdisjoint(test) and len(train|test)==23
    metrics=read(BASE/'methods/A_DF_nested_OOF/detection_metrics.csv')
    q=[r for r in metrics if r['population']=='visible']
    tp,fp,fn=(sum(int(r[k]) for r in q) for k in ('TP','FP','FN'))
    f1=2*tp/(2*tp+fp+fn)
    pc=read(BASE/'methods/A_DF_nested_OOF/plot_comparison.csv')
    ratio=statistics.mean(float(r['C_est_tC_ha']) for r in pc)/statistics.mean(float(r['C_ref_tC_ha']) for r in pc)
    saved=next(r for r in read(BASE/'method_comparison.csv') if r['method']=='A_DF_nested_OOF' and r['siteID']=='ALL')
    assert math.isclose(f1,float(saved['F1_visible']),abs_tol=1e-12)
    assert math.isclose(ratio,float(saved['estimate_reference_ratio']),abs_tol=1e-12)
    assert abs(f1-.55)<.02 and abs(ratio-.68)<.02
    fa=read(BASE/'field_trees_all.csv'); fg=read(BASE/'field_trees.csv')
    assert len(fa)==617 and len(fg)==532
    for r in pc:
        refs=[float(t['ref_AGB_kg']) for t in fa if t['plotID']==r['plotID']]
        carbon=sum(refs)/1000*10000/float(r['area_m2'])*.5
        assert math.isclose(carbon,float(r['C_ref_tC_ha']),rel_tol=1e-12)
    sites={s:sum(r['siteID']==s for r in pc) for s in ('GRSM','ORNL')}
    assert sites=={'GRSM':15,'ORNL':8}
    assert all(int(r['imagery_year'])==2021 for r in pc)
    assert all(r['imagery_provisional']=='True' for r in pc if r['siteID']=='ORNL')
    # Input stamp records and third-run package version remain identifiable.
    assert json.loads((BASE/'versions.json').read_text())['deepforest']=='2.1.0'
    assert not any(p in sys.modules for p in ('pandas','geopandas','rasterio','pyogrio','deepforest'))
    result=dict(status='PASS',cells=len(nb['cells']),guarded_legacy_cells=guarded,
        preserved_core_functions=list(names),cohort=sites,reference_trees=len(fa),coordinate_trees=len(fg),
        cached_10cm_F1_visible=f1,cached_10cm_carbon_ratio=ratio,cached_10cm_carbon_RMSE=float(saved['carbon_RMSE_tC_ha']),
        cached_patch600_subset_winners=verified,actual_51cm_gsd_for_100m_crop=100/196,
        scientific_imports_executed=False,GPU_inference_executed=False,
        remaining='실제 Rasterio 축소·DeepForest 추론·25/51cm 평가는 Colab 실행 필요')
    (HERE/'RESOLUTION_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    validation = HERE/'VALIDATION.md'
    previous = validation.read_text(encoding='utf-8').split('## 해상도 실험 추가 검증')[0].rstrip()
    validation.write_text(previous+f'''\n\n## 해상도 실험 추가 검증

- 50개 셀(코드24개), 전체 문법·기존 함수 AST 일치: 통과
- 기존 처리16개 셀을 캐시 전용 모드에서 실제로 건너뛰는 제어 흐름 검사: 통과
- 기존3차23개(GRSM15+ORNL8), 참조617그루/좌표532그루: 원본 CSV 검산
- 캐시A의 F1_visible={f1:.6f}, 탄소비율={ratio:.6f}, RMSE={float(saved['carbon_RMSE_tC_ha']):.6f}tC/ha: 재추론 없이 재집계로 대조
- 캐시A의5개 outer fold가 모두600/.2이며600 제한 안쪽 탐색에서도 모두.2가 최적: 통과
- 다른patch 캐시를600 고정으로 재표기하면 거부하는 음성 검사: 통과
- 실제GSD 계산:100m/196픽셀=0.5102040816m. 원본 범위 보존과 정수크기 규칙: 통과
- 원본10cm 네 방법 결과 재사용, 새RGB/CHM 다운로드 호출 없음, CHM고정, 임계값만 nested재선택, 완료 캐시 경로: 소스 검사
- 평균 축소·bilinear 확대·격자/nodata·캐시 재사용·실측 참조 총량·합성CHM/CV 검사는Colab 실행 셀에 포함
- Windows에서금지 패키지 import, Rasterio실행, GPU추론은 수행하지 않음.25/51cm 성능과 그림은아직미생성
- 마지막요약·/content/neon_results.zip 다운로드 셀 유지. 상세검산: RESOLUTION_VALIDATION.json
''',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
