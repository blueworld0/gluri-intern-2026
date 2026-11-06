# Week05 모델 개선 탐색 실험 보고서

## 1. 기존 모델의 문제와 자료

Week04는 SCL 통과 GEDI L4A 90점의 `agbd`(Mg/ha)를 목표로 `ndvi`, `evi`, `ndmi`만 사용했다. RF는 `n_estimators=100`, `random_state=42`이며 깊이·잎 크기 제한이 없다(`week04/train_baseline.py:27-32,116-130`). Random 80/20 분할은 72/18점, Spatial은 동쪽 끝 1km 블록 20점을 평가하고 500m 버퍼 15점을 제외해 55점으로 학습한다(`week04/train_baseline.py:80-113`, `week04/compare_baselines.py:30-40`). Week04 RF의 R²는 두 분할 모두 음수였고 Spatial에서는 Dummy보다 세 지표가 나빴다. Sentinel-2 `/65535` 변환도 물리적으로 미검증이다(`week05/step3/report.md`).

## 2. 실행한 개선

**A — 입력 변수:** 기존 세 지수에 없는 `B03`, `B12` 원시 TIFF 값을 기존 90점의 행·열에서 직접 추출했다. B08과 B12로 `NBR=(B08−B12)/(B08+B12)`를 계산했다. 필요한 밴드 3개(B03/B08/B12)는 모두 존재하고 90점에서 마스크·분모가 유효함을 실행 중 확인했다. NBR은 동일 양의 공통 배율·영 offset에서 불변이지만 밴드별 offset·비선형 출력에서는 불변이 아니다. B03·B12 원시값의 절대 물리 단위도 검증되지 않았다. 따라서 이 후보는 **기존 파일값을 사용한 탐색적 입력**이며 검증된 새 반사율 자료로 취급하지 않는다. 후보는 `기존+NBR`, `기존+B03/B12`, `기존+NBR+B03/B12` 세 가지였고, 외부 평가점이 아닌 학습 내부 CV RMSE로 하나를 선택했다. B02/B04/B08/B11은 이미 지수 계산에 쓰이므로 '기존에 사용하지 않은 밴드' 후보로 추가하지 않았다.

**B — 복잡도:** 기본 3지수를 유지하고 `max_depth`, `min_samples_leaf`, `max_features`의 축소 후보 세 개만 비교했다. 각 외부 학습 집합에서 500m easting 그룹을 사용하는 3-fold GroupKFold의 pooled out-of-fold RMSE가 최소인 구성을 선택했다. 파라미터는 평가점에서 선택하지 않았다. 전체 검색공간과 fold별 크기는 [`experiment_design.json`](experiment_design.json)에 있다.

**C — GEDI 품질:** 최종 CSV에는 일반명 `quality_flag` 대신 V003의 `l2a_quality_flag_rel3`, `l4a_quality_flag_rel3`, `degrade_flag`, `degrade_include_flag`, `elev_highestreturn_outlier_flag`, `predictor_limit_flag`, `response_limit_flag`, `sensitivity`, `agbd_se` 등이 있다. 90점에서는 두 rel3 플래그가 모두 1, `degrade_flag=0`, `degrade_include_flag=True`, elevation outlier 및 predictor/response limit 플래그가 모두 0이다([`quality_flag_summary.csv`](quality_flag_summary.csv)). 기존 `week03/analyze_gedi.py:84-99` 필터와 SCL 선별을 이미 통과했으므로 이진 플래그를 다시 걸어도 차이가 없다. [GEDI L4A V3 공식 가이드](https://daac.ornl.gov/GEDI/guides/GEDI_L4A_AGB_Density_V3.html)는 sensitivity가 파형 SNR에 따른 최대 관통 가능 수관 피복을 나타내며 지역에 따라 품질·표본 수의 절충을 검토할 수 있다고 설명한다. **0.95를 추가 민감도 조건**으로 사전 지정해 학습점에만 적용하고 평가점은 바꾸지 않았다. 이보다 낮은 점을 무효라고 판정한 것은 아니다.

## 3. 성능 비교

전체 수치는 [모델 비교표](model_comparison.md)와 [`experiment_results.csv`](experiment_results.csv)에 있다.

- **A:** Random 평가 RMSE 52.34→52.78, MAE 42.87→43.03으로 악화했다. Spatial은 RMSE 70.61→69.52, MAE 57.71→56.27로 소폭 개선했으나 Dummy RMSE 57.16보다 여전히 높다. 두 분할에서 **서로 다른 후보**가 학습 CV로 선택됐다.
- **B:** Random 평가 RMSE 52.34→49.65, MAE 42.87→41.17, R² −0.139→−0.025였다. 기존 Random Dummy(RMSE 50.94, MAE 43.34)보다도 낮았다. Spatial은 RMSE 70.61→62.35, MAE 57.71→50.47, R² −2.165→−1.468로 개선했지만 **Spatial Dummy(RMSE 57.16, MAE 48.07)에는 못 미쳤다**. R²는 여전히 음수다.
- **C:** Random 평가점 18개는 그대로 두고 학습점을 72→39개로 줄였다. RF RMSE는 52.34→67.62, MAE는 42.87→54.35로 악화했다. 변경된 학습점으로 학습한 Dummy도 RMSE 64.15, MAE 52.59로 악화해, 단순히 RF의 문제로 해석할 수 없다. Spatial에서는 55개 학습점 중 10개만 조건에 남아 사전 최소 기준 30개 미만이므로 **실험하지 않았다**.

## 4. 품질 조건의 표본 차이

전체 90점 중 sensitivity≥0.95는 45점이다. 이 45점의 AGBD 평균은 132.40 Mg/ha이고 나머지 45점은 75.48 Mg/ha로 서로 다르다. Random 학습 전체 평균은 106.71이지만 필터 후 134.25이며, 평가점 전체 평균은 92.88이다. Spatial 학습 전체 55점 중 조건 통과는 10점(평균 154.15), 평가 20점은 모두 조건을 통과한다. 따라서 필터 실험의 점수 차이는 **품질, 학습 규모, AGBD 분포 변화가 함께 작용한 결과**로 분리할 수 없다. 원자료 분포는 [`quality_distribution.csv`](quality_distribution.csv)에 기록했다. Spatial 버퍼 제외 15점은 어떤 실험에도 학습·평가로 사용하지 않았다.

## 5. 한계와 추가 검증

- 외부 평가점은 Week04와 Step1에서 이미 반복적으로 관찰했다. 이번 선택 로직은 학습 내부에서만 동작하지만, 연구 과정 전체에는 기존 평가점에 대한 사전 지식이 있다. **독립 최종 일반화 성능·유의성 검정이 아니다.**
- 90점은 한 GEDI 파일의 제한된 위치·시점이고 외부 평가점은 Random 18, Spatial 20점뿐이다. 내부 CV 그룹도 불균형하며 500m 버퍼가 없다. 하나의 동쪽 Spatial 홀드아웃으로 다른 지역 성능을 주장할 수 없다.
- 추가 원시 밴드와 NBR은 방사 변환이 미검증인 동일 TIFF에서 왔다. 특히 EVI는 공통 배율에도 달라지며, 밴드별 offset·클리핑이 있다면 NBR·NDVI·NDMI도 영향을 받을 수 있다. 물리적 자료 신뢰성 검증이 선행돼야 한다.
- AGBD는 GEDI L4A 모델 추정값이며 독립 현장 실측값이 아니다. `agbd_se`도 GEDI 예측의 불확실성 지표이지 현장 오차가 아니다.

## 6. 후속 연구 우선순위

1. 2025 Sentinel-2 출력 설정/evalscript와 동일 장면의 검증된 FLOAT32 반사율을 확보해 방사 변환을 확정한다(`week05/step3/report.md`). 확인된 자료 버전을 따로 만든 뒤 현재 결과와 비교한다.
2. 독립 장면·GEDI 궤도를 늘려 위치와 시기의 대표성을 확보한다. 모델 선택은 학습 내부의 **버퍼를 둔** 공간 CV로 진행하고, 마지막까지 사용하지 않은 지역·기간에서 RF와 Dummy를 다시 비교한다.
3. GEDI sensitivity와 AGBD·위치의 관계를 별도로 분석해 품질 조건의 선택 편향을 평가한다. 충분한 표본을 확보하기 전에는 강한 품질 필터로 성능이 개선된다고 주장하지 않는다.

## 재현과 변경 범위

프로젝트 루트에서 `week04/.venv/Scripts/python.exe week05/improvement/improve_model.py`를 실행하면 결과 CSV, 예측값 CSV, 내부 CV 기록, 품질 분포 CSV, JSON 설계 기록, PNG 그래프를 이 디렉터리에 다시 생성한다. Week03–Week05 기존 코드와 원본 데이터는 수정하지 않는다. 나무 수·난수 시드, 외부 분할과 평가점 ID는 Week04와 동일하며 기존 RF 지표·검증 예측값의 재현을 1e-8 허용오차로 검사한다.
