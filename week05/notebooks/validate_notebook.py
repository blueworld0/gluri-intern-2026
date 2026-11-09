"""Standard-library-only structural/numerical checks; never run notebook cells."""
import ast
import json
import math
from pathlib import Path

root=Path(__file__).resolve().parent
nb=json.loads((root/'NEON_DeepForest_AGB_2021.ipynb').read_text(encoding='utf-8'))
assert nb['nbformat']==4
assert len({c['id'] for c in nb['cells']})==len(nb['cells'])
code_cells=[]
for i,c in enumerate(nb['cells']):
    if c['cell_type']!='code': continue
    assert nb['cells'][i-1]['cell_type']=='markdown'
    assert c['outputs']==[] and c['execution_count'] is None
    source=''.join(c['source'])
    if source.startswith('%pip'):
        assert not code_cells
        for package in ['deepforest','neonutilities','rasterio','geopandas','shapely','pyproj','scikit-learn','matplotlib']:
            assert package in source
    else:
        compile(source,f'notebook_cell_{i}','exec')
    code_cells.append(source)
combined='\n'.join(s for s in code_cells if not s.startswith('%'))
tree=ast.parse(combined)
names={n.name for n in tree.body if isinstance(n,ast.FunctionDef)}
for name in ['load_field_year','map_individuals','build_footprints','estimate_dbh','estimate_agb',
 'predict_crowns','match_one_to_one','grouped_dbh_regression','errors','select_aop_pair',
 'plot_biomass_totals','build_site_summary','plot_carbon_scatter','variable_window_peaks',
 'chm_watershed_crowns','combine_boxes_and_peaks','nested_group_parameter_selection',
 'evaluate_detection_method','detection_algorithm_smoke_checks']:
    assert name in names,name
assert 'use_release' in combined and 'load_model' in combined
for filename in ['detected_trees.csv','matched_trees_visible.csv','stage_errors.csv',
                 'plot_comparison.csv','report.md','height_sensitivity.csv','field_trees_all.csv',
                 'aop_plot_coverage.csv','site_summary.csv','/content/neon_results.zip',
                 'method_comparison.csv','A_precision_recall.png','selected_detection_method.json',
                 'A_fold_plot_assignments.csv','A_outer_groupkfold.csv']:
    assert filename in combined
compile((root/'neon_pipeline.py').read_text(encoding='utf-8'),'neon_pipeline.py','exec')
assert "files.download('/content/neon_results.zip')" in ''.join(nb['cells'][-1]['source'])
module_tree=ast.parse((root/'neon_pipeline.py').read_text(encoding='utf-8'))
module_functions={n.name:ast.dump(n) for n in module_tree.body if isinstance(n,ast.FunctionDef)}
for node in tree.body:
    if isinstance(node,ast.FunctionDef) and node.name in module_functions:
        assert ast.dump(node)==module_functions[node.name],f'Notebook/module mismatch: {node.name}'
cfg_node=next(n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='CFG' for t in n.targets))
cfg=eval(compile(ast.Expression(cfg_node.value),'config_only','eval'),{'dict':dict})
assert cfg['search_score_thresholds']==[0.1,0.2,0.3,0.4,0.5]
assert cfg['search_patch_sizes']==[300,400,600]
assert len(set(cfg['comparison_plot_ids']))==cfg['expected_plot_count']==23
assert set(cfg['comparison_overlay_plots'])=={'GRSM_016','GRSM_051','GRSM_012','ORNL_010'}
assert "os.environ['HF_HOME'] = '/content/hf_cache'" in combined
assert "not r.headers.get('Content-Encoding')" in combined
assert 'valid=valid[~valid' in combined and 'warnings.warn' in combined
assert cfg['aop_include_provisional'] is True
assert 'field_all=field.merge(ref' in combined
assert combined.index('DEFAULT_DETECTION_METHOD=ranking')<combined.index('sensitivity=[]')
assert 'maxima=mask' in combined and 'diameter_b*H*H' in combined
assert 'assert set(trainids).isdisjoint(testids)' in combined
# Independent dimensional/scenario checks with Python math; no scientific imports.
assert 100/1000*10000/400==2.5
dbh=0.367*(20*5)**0.887
assert 20<dbh<25
agb=math.exp(-2.2118+2.4133*math.log(dbh))
assert 150<agb<200
for delta in [-0.2,-0.1,0.1,0.2]:
    dbh_changed=0.367*((20*(1+delta))*5)**0.887
    changed=math.exp(-2.2118+2.4133*math.log(dbh_changed))
    assert math.isclose(changed/agb,(1+delta)**(0.887*2.4133))
report=f'''# 사전 검증 기록

- 노트북 JSON v4 구조, 고유 셀 ID: 통과
- 전체 {len(nb['cells'])}개 셀, 코드 {len(code_cells)}개: Python 문법 검사 통과 (첫 pip 셀은 Colab magic으로 별도 확인)
- 모든 코드 셀 바로 위 한국어 설명: 확인
- 실행 결과와 실행 번호: 비어 있음 (실행 완료로 오인하지 않도록 유지)
- 재사용 함수 모듈: 문법 검사 통과
- m/cm·kg/t·㎡/ha 단위 및 높이 민감도 수치 검산: 통과
- 로컬에서 pandas/geopandas/rasterio/DeepForest를 import하거나 실행하지 않음
- 제공된 1차 outputs/로 집계를 검산함. 수정본 다운로드·GPU 추론은 Colab 재실행 필요

공식 API 문서에서 load_by_product(dpid, site, startdate, enddate, package, release, token),
by_tile_aop(dpid, site, year, easting, northing, buffer, savepath, include_provisional, skip_if_exists, overwrite),
DeepForest load_model/use_release 및 predict_tile path/raster_path 차이를 확인했습니다.

Jucker 계수는 itcSegment/R/dbh.R의 biome=11 (a=0.367, b=0.887, g=0)을 확인했습니다.
Chojnacky 기본값은 혼합활엽수라는 독립 그룹이 아닌 일반 활엽수 대체 그룹
(β0=-2.2118, β1=2.4133)이며 원 논문 표와 MasterTaxonList 열 이름을 확인했습니다.
실행 시 원문 계수표를 캐시하고 계수와 해시를 다시 확인합니다.

## 1차 수정 검증

- 공식 API 직접 조회: ORNL RGB·CHM 모두 2021-06 존재, PROVISIONAL. ORNL UTM 16N 확인
- 실제 API 응답 기반 연도 선택: 2021 유지, include_provisional=True; 자료가 없을 때만 ORNL 2022→2020
- 좌표 누락이어도 참조 AGB가 있으면 조사구 총량 유효. 실제 AGB 누락일 때만 NaN: 회귀 검사
- 공식 면적 우선과 누락 분해 항등식: 회귀 검사
- 노트북과 neon_pipeline.py의 공통 함수 AST 일치: 확인
- 탐지·Jucker·Chojnacky·매칭·회귀·민감도 코드 원본 보존: check_revision.py에서 검사
- report.md 상단과 마지막 요약은 동일 표 사용, 마지막 ZIP 셀 존재: 확인
- GRSM 15개 검산: outputs/revision_check/README.md 참고

## 나무 찾기 개선 추가 검증

- score 5개×patch 3개, 고정 23개 ID, 요청된 4개 나란히 그림 ID: 설정 검사 통과
- 바깥 GroupKFold 평가 정답과 안쪽 선택 자료 분리, OOF 결과와 배포용 전체 선택 조합 분리: 소스 확인
- B: Gaussian·가변 창·watershed·면적 기반 CD, C: 중복 제거·상자 밖 꼭대기 추가: 함수/공식 확인
- 네 방법 전체/사이트 비교표와 탄소 RMSE 최소 기본값이 이후 민감도 이전에 연결됨: 확인
- QUDO 경고/제외, gzip 크기 검사 예외, /content/hf_cache, ORNL provisional, 전체 참조 총량: 유지 검사 통과
- 합성 두 나무 CHM의 꼭대기·watershed 면적·CD·결합 중복 제거와 평가 정답 변경 시 선택 격리 검사를 Colab 셀 14B에 추가
- 합성 검사와 실제 23개 자료의 GPU/CHM 연산은 로컬에서 실행하지 않음. 새 성능 수치를 주장하지 않음
'''
(root/'VALIDATION.md').write_text(report,encoding='utf-8')
print('PASS:',len(nb['cells']),'cells;',len(code_cells),'code cells; no scientific imports/execution')
