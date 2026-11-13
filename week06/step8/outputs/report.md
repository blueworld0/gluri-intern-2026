# Week06 Step 8 — 확대 데이터 기반 모델 성능 평가

## 목적과 기본 결과

산림으로 독립 검증되지 않은 Sentinel-2 SCL 식생 분류 통과 GEDI L4A 표본 441점만 평가했다. 기존 90점을 합치지 않았다. 기본 결과는 EVI 이상점도 포함한다. 전체 OOF RMSE가 가장 낮은 모델은 **Week05 Spatial RF (44.3820 Mg/ha)**이다. 최선이라는 표현은 이번 고정 비교의 관측 지표에 한정한다.

| model | n | rmse | mae | r2 |
| --- | --- | --- | --- | --- |
| week05_spatial_B_fixed | 441 | 44.3820 | 29.0806 | 0.0722 |
| week04_rf_baseline | 441 | 46.5713 | 30.0932 | -0.0216 |
| dummy_mean | 441 | 47.4352 | 36.5294 | -0.0598 |

Dummy 대비 RMSE 감소량은 Week04 RF **0.8639 Mg/ha (1.82%)**, Week05 Spatial RF **3.0532 Mg/ha (6.44%)**이다. 음수는 악화다. 기존 RF 대비 개선 RF의 RMSE 변화(개선−기존)는 **-2.1893 Mg/ha**이다. pooled OOF 지표는 각 점의 한 번의 평가 예측을 모아 직접 계산했으며 fold 지표 평균과 다르다.

## 고정 실험 설정과 검증

입력: `week06/step7/outputs/new441_evaluation_input.csv`. NDVI·EVI·NDMI → AGBD(Mg/ha). EPSG:32652의 원점 0 기준 5km 격자, 기존 5-fold 배정, 평가점 중심 500m 미만 학습점 제외를 그대로 재사용했다. 차단된 점을 다른 fold에서 임의로 재배정하지 않았다. 점·shot·좌표·공통 10m 셀 중복은 없고, 모든 441점이 정확히 한 번 평가된다. 학습/평가 블록은 겹치지 않고 모든 학습/평가점 거리가 500m 이상이다. 세 모델의 평가 ID는 fold 및 실험 조건마다 같다.

| fold | train_n | test_n | buffer_excluded_n | minimum_train_test_distance_m | outlier_role |
| --- | --- | --- | --- | --- | --- |
| 1 | 321 | 119 | 1 | 513.5205 | train |
| 2 | 349 | 79 | 13 | 513.4145 | test |
| 3 | 353 | 87 | 1 | 634.9766 | train |
| 4 | 344 | 80 | 17 | 513.6963 | train |
| 5 | 362 | 76 | 3 | 571.2422 | train |

- Dummy: 각 fold 학습 AGBD의 평균만 예측한다.
- Week04 RF: n_estimators=100, max_depth=None, min_samples_leaf=1, max_features=1.0.
- Week05 개선 RF: **Spatial B**의 max_depth=6, min_samples_leaf=3, max_features=0.67을 고정했다. Week05 Random B(depth=3, leaf=4, max_features=1.0)는 적용하지 않았다. 근거는 `week05/improvement/model_comparison.md`, `experiment_design.json`, Step7 `evaluation_plan.json`이다.
- RF 공통 random_state=42, n_jobs=1. 나머지 sklearn 기본값은 `experiment_config.json`의 전체 get_params로 기록했다. 전체 데이터에 기반한 스케일링·결측 대치·변수 선택·신규 탐색은 없다. 저장된 기존 모델을 불러온 평가가 아니라 고정 train ID에서의 재학습이다.
- Python 3.12.2, numpy 2.5.4, pandas 3.1.0, sklearn 1.9.1. 30회 본/민감도 적합과 fold1 RF 두 모델의 결정적 재실행 검증 2회, 총 32회 적합했다. 재실행 예측 차이는 0이다. 원본/분할 파일 SHA-256은 실행 전후 동일했다.

## fold별 편차

| model | fold | train_n | n | rmse | mae | r2 |
| --- | --- | --- | --- | --- | --- | --- |
| dummy_mean | 1 | 321 | 119 | 44.7348 | 34.1687 | -0.0120 |
| week04_rf_baseline | 1 | 321 | 119 | 45.8131 | 27.5502 | -0.0614 |
| week05_spatial_B_fixed | 1 | 321 | 119 | 41.4221 | 25.6412 | 0.1323 |
| dummy_mean | 2 | 349 | 79 | 62.8822 | 45.4689 | -0.3082 |
| week04_rf_baseline | 2 | 349 | 79 | 57.1367 | 39.9934 | -0.0801 |
| week05_spatial_B_fixed | 2 | 349 | 79 | 56.6692 | 40.2307 | -0.0625 |
| dummy_mean | 3 | 353 | 87 | 43.5762 | 36.0604 | -0.0415 |
| week04_rf_baseline | 3 | 353 | 87 | 39.5782 | 25.5526 | 0.1408 |
| week05_spatial_B_fixed | 3 | 353 | 87 | 38.7081 | 25.7474 | 0.1782 |
| dummy_mean | 4 | 344 | 80 | 40.6196 | 30.0312 | -0.1254 |
| week04_rf_baseline | 4 | 344 | 80 | 45.7046 | 32.8668 | -0.4249 |
| week05_spatial_B_fixed | 4 | 344 | 80 | 43.2449 | 31.0316 | -0.2756 |
| dummy_mean | 5 | 362 | 76 | 43.6150 | 38.3101 | -0.5259 |
| week04_rf_baseline | 5 | 362 | 76 | 43.7428 | 26.0624 | -0.5349 |
| week05_spatial_B_fixed | 5 | 362 | 76 | 41.4888 | 24.6374 | -0.3808 |

| model | rmse_min | rmse_max | rmse_sample_sd | r2_min | r2_max |
| --- | --- | --- | --- | --- | --- |
| dummy_mean | 40.6196 | 62.8822 | 8.9614 | -0.5259 | -0.0120 |
| week04_rf_baseline | 39.5782 | 57.1367 | 6.5135 | -0.5349 | 0.1408 |
| week05_spatial_B_fixed | 38.7081 | 56.6692 | 7.0987 | -0.3808 | 0.1782 |

표본 수와 AGBD 분포가 다른 공간 fold 간 편차를 함께 본다. `rmse_sample_sd`는 5개 fold RMSE의 기술적 표준편차(ddof=1)이며 독립 반복 실험의 표준오차나 신뢰구간이 아니다. 음의 R²는 해당 평가집합 평균 기준 제곱오차보다 더 큰 오차를 뜻한다. Dummy 자체는 학습 평균을 쓰므로 R²=0일 필요가 없다.

## EVI 이상점 민감도

RMSE 기준 기존 RF는 5개 중 **2개**, 개선 RF는 **4개** fold에서 Dummy보다 좋았다. 개선 RF는 기존 RF보다 모든 fold의 RMSE가 낮지만, fold4에서는 Dummy보다 나쁘다. 개선 RF의 fold2·4·5 R²는 음수다. 개선 RF의 MAE도 pooled 기준 낮지만 fold2·3·4에서는 기존 RF보다 높다. 따라서 지역 및 오차 지표 전반에서 안정적인 개선을 확보했다고 판단할 수 없다.

관측 ID: `GEDI04_A_2025115163530_O36071_02_T10090_02_004_02_V003.h5:BEAM0110:360710600200387796`. EVI=-41.465836, AGBD=9.350689 Mg/ha, EVI 분모=-0.013300. SCL=4가 지수 안정성을 보장하지 않는다. 값이 극단적이라는 이유만으로 물리적 오류라고 단정하지 않았다.

이 점은 fold2에서는 평가점이고 fold1/3/4/5에서는 학습점이다. 기본 실험은 포함, 민감도만 해당 ID를 모든 역할에서 제외했다. 원래 fold 및 buffer_excluded 배정은 유지했다. 따라서 fold2는 학습 349점/평가 78점, 다른 fold는 학습 1점 감소/평가 수 유지다. 500m 버퍼를 다시 산정하거나 이전 제외점을 복원하지 않았다.

기본 실험에서 이상점의 예측:

| model | fold | agbd | prediction | residual |
| --- | --- | --- | --- | --- |
| dummy_mean | 2 | 9.3507 | 36.5327 | -27.1820 |
| week04_rf_baseline | 2 | 9.3507 | 11.7757 | -2.4250 |
| week05_spatial_B_fixed | 2 | 9.3507 | 14.7716 | -5.4209 |

| model | n_included | n_excluded | included_rmse | base_common_rmse | excluded_rmse | evaluation_removal_delta_rmse | refit_common_delta_rmse | total_delta_rmse |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dummy_mean | 441 | 440 | 47.4352 | 47.4714 | 47.4837 | 0.0362 | 0.0122 | 0.0484 |
| week04_rf_baseline | 441 | 440 | 46.5713 | 46.6241 | 46.7940 | 0.0527 | 0.1700 | 0.2227 |
| week05_spatial_B_fixed | 441 | 440 | 44.3820 | 44.4316 | 44.3776 | 0.0497 | -0.0540 | -0.0044 |

`base_common_rmse`는 **재학습 없이** 기본 예측에서 해당 평가점만 뺀 공통 440점의 지표다. `evaluation_removal_delta_rmse`는 평가집합 변화 효과, `refit_common_delta_rmse`는 공통 440점에서 학습점 제외 후 재학습 효과다. 두 차이의 합은 전체 RMSE 변화와 같다. RMSE는 비선형 지표이므로 이 차이를 개별 관측점의 인과 기여로 해석하지 않는다. fold2의 공통 78점에서는 세 모델 모두 예측 변화가 정확히 0이며, 지표 변화는 평가점 제외 때문이다. 다른 fold의 변화는 학습에서 이 점을 제외한 영향이며 RF 부트스트랩 표본의 구성 변화도 함께 포함한다. 따라서 변화 전부를 EVI 수치 자체의 효과로 단정하지 않는다. 전체 및 fold별 RMSE·MAE·R² 분해는 `sensitivity_results.csv`에 저장했다. 민감도 결과로 기본 실험을 대체하거나 설정을 다시 선택하지 않았다.

## AGBD 구간별 편향

이상점 제외 후 기존 RF RMSE는 **46.7940 Mg/ha**로 +0.2227 변하고, 개선 RF는 **44.3776 Mg/ha**로 -0.0044 변한다. 개선 RF의 평가점 제거 효과는 +0.0497, 공통점 재학습 효과는 -0.0540 Mg/ha로 일부 상쇄된다. 이상점 제거만으로 성능 문제가 해결되지 않으며 기본 모델 순위도 유지된다.

잔차는 **관측 GEDI AGBD − 예측**이다. 양수는 과소 예측, 음수는 과대 예측이다. 다음은 기본 OOF 결과이며, 작은 구간 표본 수를 함께 확인해야 한다.

| model | agbd_bin | n | mean_agbd | mean_prediction | mean_residual | rmse |
| --- | --- | --- | --- | --- | --- | --- |
| dummy_mean | [0,50) | 311 | 17.8458 | 43.1069 | -25.2611 | 28.1619 |
| dummy_mean | [50,100) | 74 | 71.2897 | 39.9464 | 31.3433 | 35.1348 |
| dummy_mean | [150,200) | 17 | 171.5559 | 41.0650 | 130.4908 | 131.5799 |
| dummy_mean | [100,150) | 35 | 120.9554 | 40.1843 | 80.7712 | 81.9592 |
| dummy_mean | [200,inf) | 4 | 219.3076 | 43.1533 | 176.1543 | 176.6788 |
| week04_rf_baseline | [0,50) | 311 | 17.8458 | 34.8697 | -17.0239 | 31.7431 |
| week04_rf_baseline | [50,100) | 74 | 71.2897 | 60.4417 | 10.8480 | 27.8401 |
| week04_rf_baseline | [150,200) | 17 | 171.5559 | 51.0933 | 120.4626 | 123.2874 |
| week04_rf_baseline | [100,150) | 35 | 120.9554 | 52.0104 | 68.9450 | 74.4593 |
| week04_rf_baseline | [200,inf) | 4 | 219.3076 | 36.7954 | 182.5122 | 182.5583 |
| week05_spatial_B_fixed | [0,50) | 311 | 17.8458 | 34.8723 | -17.0264 | 29.0684 |
| week05_spatial_B_fixed | [50,100) | 74 | 71.2897 | 59.7106 | 11.5791 | 25.2590 |
| week05_spatial_B_fixed | [150,200) | 17 | 171.5559 | 51.4332 | 120.1227 | 122.1522 |
| week05_spatial_B_fixed | [100,150) | 35 | 120.9554 | 54.5011 | 66.4543 | 71.2725 |
| week05_spatial_B_fixed | [200,inf) | 4 | 219.3076 | 41.0108 | 178.2967 | 178.3351 |

고 AGBD 구간의 표본 부족과 평균 회귀 경향을 확인하는 진단표다. 구간 경계는 성능 선택에 사용하지 않았으며 `[하한, 상한)`이다. 모든 조건의 상세 지표는 `agbd_bin_metrics.csv`에 있다.

## 해석의 한계와 다음 단계

실제로 50 Mg/ha 미만 표본이 311/441점(70.5%)이다. 개선 RF는 이 구간에서 평균 약 17.03 Mg/ha 과대 예측하며, 100–150 구간에서는 약 66.45, 150–200에서는 약 120.12, 200 이상에서는 약 178.30 Mg/ha 과소 예측한다. 200 이상은 4점뿐이며 이 구간의 RF RMSE는 두 모델 모두 Dummy보다 크다. 높은 AGBD를 구분하는 능력이 부족하고 평균 쪽으로 예측이 수축되는 문제가 남는다.

기존 90점은 Browser TIFF의 /65535 반사율 가정이 미검증이며, 신규 441점은 SAFE XML의 offset=-1000, quantification=10000을 확인해 처리했다. 품질 마스크 적용 기준도 다르다. 따라서 두 집합을 합치지 않았고 Week04/05의 기존 90점 성능 수치와 직접 개선률을 계산하지 않았다.

신규 표본도 기존과 동일한 2025-04-25·O36071 GEDI 관측에 속한다. 이번 평가는 독립 궤도 검증이 아니다. AGBD는 GEDI L4A 모델 추정치이며 독립 현장 실측값이 아니다. 대상은 산림 확정 표본이 아닌 SCL 식생 후보이며 농경지 등 다른 식생이 포함될 수 있다. 5km 격자·500m 버퍼는 공간 독립성을 보증하는 거리나 자기상관 범위 추정치가 아니다. 이 공간 CV를 독립 최종 일반화 성능이나 통계적 유의성의 증거로 해석하지 않는다.

다음 단계에서는 (1) EVI 분모·원본 픽셀과 주변 품질을 점검해 지수 안정성 처리 규칙을 사전에 정하고, (2) 고 AGBD 표본 부족 및 지역별 잔차를 점검하고, (3) 허용 시 독립 관측일·궤도와 검증된 산림 자료를 추가하며, (4) 기존 90점의 반사율/품질 조건을 통일한 뒤 결합 가능성을 재검토한다. RF가 Dummy보다 좋더라도 이 한계와 fold 편차는 남는다. 이번 단계에서는 추가 자료 확보나 튜닝을 수행하지 않았다.

## 재현 및 파일

저장소 루트의 PowerShell에서 다음을 실행한다. 출력 폴더에 파일이 있으면 덮어쓰지 않고 중단하므로 재현 시 새로운 경로를 지정한다.

```powershell
& .\week04\.venv\Scripts\python.exe .\week06\step8\evaluate_models.py --output-dir .\week06\step8\reproduction_outputs
```

의존성은 위 버전 및 기존 Week04 환경을 사용한다. 설치·네트워크 연결은 코드에서 하지 않는다. `fold_metrics.csv`/`overall_metrics.csv`/`predictions.csv`의 case는 기본 `main_all441`, 민감도 `sensitivity_exclude_evi`를 구분한다. 지표 단위는 RMSE·MAE·잔차 Mg/ha, R² 무차원이다. `validation.json`은 입력 해시, 분할 검증, 재현성 및 지표 독립 재계산 결과를 담는다. `experiment_config.json`은 전체 모델 파라미터와 버전을, `fold_validation.csv`는 학습/평가 수와 최소 거리를 기록한다.
