# Week05 Step1: Week04 모델 성능 및 예측 오차 분석

## 자료와 방법

- Week03 SCL 통과 90점의 GEDI L4A 추정 AGBD(Mg/ha)를 목표값으로, Sentinel-2 NDVI·EVI·NDMI를 입력으로 사용했다. GEDI AGBD는 현장 실측값이 아니다.
- Week04의 RF 지표·검증점 예측 CSV를 재사용했다. Dummy는 Week04의 분할 함수와 학습 집합 평균으로 예측값을 재현했다. 네 실험의 지표를 기존 CSV와 1e-8 이내로 대조했다.
- 잔차 = 실제 GEDI 추정 AGBD − 모델 예측값. 양수는 과소 예측이다. 원본 데이터·Week04 결과는 변경하지 않았다.

## Week04 실험 결과

| 분할 | 모델 | 학습 n | 평가 n | R² | RMSE | MAE |
|---|---|---:|---:|---:|---:|---:|
| random_80_20 | RandomForest | 72 | 18 | -0.139 | 52.34 | 42.87 |
| random_80_20 | DummyMean | 72 | 18 | -0.080 | 50.94 | 43.34 |
| spatial_block_500m_buffer | RandomForest | 55 | 20 | -2.165 | 70.61 | 57.71 |
| spatial_block_500m_buffer | DummyMean | 55 | 20 | -1.074 | 57.16 | 48.07 |

RMSE·MAE 단위는 Mg/ha이다. Random 분할에서 RF의 MAE는 Dummy보다 0.47 낮지만 RMSE는 1.39 높다. Spatial 분할에서 RF는 Dummy보다 MAE 9.64, RMSE 13.45 높다. 두 분할 모두 RF R²가 음수이며 단일 분할 결과다.

## Random Split과 Spatial Split

- Random: 학습 72점, 평가 18점, 최소 학습·평가 거리 55.2m.
- Spatial: 동쪽 끝 1km easting 블록 20점을 평가하고, 500m 버퍼의 15점을 제외해 55점으로 학습했다. 최소 거리 500.6m. 두 분할의 평가점·학습 규모가 다르므로 점수 차이 전부를 공간 효과로 단정할 수 없다.
- 공간 분할은 좌표만으로 정해지며 학습·평가 point ID 및 raster pixel 중복은 0개다. 검증점 AGBD를 분할 선택에 사용하지 않았다. 동일 GEDI 궤도와 공간 자기상관, 분석 과정의 반복적 모델 선택에 따른 평가 편향 가능성은 남아 있다. 이 코드는 500m 밖의 자기상관 소멸을 증명하지 않는다.

## 예측 오차

- random_80_20, low_actual_le_50: n=3, 평균 잔차 -71.11 Mg/ha, 과소 예측 0점, 과대 예측 3점.
- random_80_20, high_actual_ge_150: n=2, 평균 잔차 +76.21 Mg/ha, 과소 예측 2점, 과대 예측 0점.
- random_80_20, very_high_actual_ge_200: n=1, 평균 잔차 +102.60 Mg/ha, 과소 예측 1점, 과대 예측 0점.
- spatial_block_500m_buffer, low_actual_le_50: n=1, 평균 잔차 -36.30 Mg/ha, 과소 예측 0점, 과대 예측 1점.
- spatial_block_500m_buffer, high_actual_ge_150: n=7, 평균 잔차 +98.98 Mg/ha, 과소 예측 7점, 과대 예측 0점.
- 높은/낮은 AGBD 구간은 사전에 정한 기술적 기준(150 이상 / 50 이하 Mg/ha)이며, 표본이 작아 경향의 일반화는 어렵다. 잔차와 예측값의 관계는 별도 산점도에서 확인할 수 있다.

## 데이터 품질과 결합 점검

- 최종 자료 90행·70열. 결측이 있는 컬럼: agbd_pi_lower_Mg_ha 9개, scl_class 90개. 학습 입력·목표·좌표 결측은 0개다.
- AGBD: 최소 11.00, 최대 301.97, 평균 103.94, 중앙값 103.08, 표본 표준편차 60.16 Mg/ha. IQR 1.5배 상한 290.55 초과 1점; 200 이상 4점, 250 이상 1점. 제거하지 않았다.
- NDVI: 범위 0.301–0.909, 중앙값 0.810, IQR 이상치 표시 7점. 분포는 `feature_distributions.png` 참조.
- EVI: 범위 0.145–0.772, 중앙값 0.509, IQR 이상치 표시 6점. 분포는 `feature_distributions.png` 참조.
- NDMI: 범위 0.050–0.493, 중앙값 0.387, IQR 이상치 표시 11점. 분포는 `feature_distributions.png` 참조.
- NDVI·NDMI의 [-1, 1] 밖 값 0개, |EVI|>2 값 0개. 최종점의 raw 65535 표시 0개, raw 0 표시 1개, SCL 경계 판정 불일치 0개. 단순 범위 검사만으로 반사율의 물리적 단위는 확인되지 않는다.
- 최종점의 GEDI–Sentinel 시차는 3–3일, GEDI 원본 파일 1개, 중복 point ID 0개, 중복 raster pixel 0개다. 원시 결합 과정은 좌표를 UTM 52N으로 변환해 지수 raster 픽셀을 추출했고, SCL 필터로 최종 90점을 선택했다.
- **반사율 단위 미확정:** Week03은 Copernicus Browser Analytical Raw uint16 TIFF를 `/65535`로 변환했다. Week04 TIFF 메타데이터에는 실질적인 radiometric scale/offset 근거가 없고, `/10000`은 민감도 비교 시나리오일 뿐이다. NDVI·NDMI는 동일 배율이면 값이 유지되지만 EVI는 상수항 때문에 달라진다. Week04 90점 비교에서 `/65535` 재계산값과 저장 지수는 거의 일치했으나 이는 물리적 단위의 검증이 아니다.
- 2025 TIFF 원시 중앙값은 B02 1848, B04 2549, B08 21830이다. `/65535` 가정 시 각각 0.028, 0.039, 0.333이지만 물리 반사율임이 입증되지는 않았다. 2025 TIFF 유효 픽셀의 65535도 존재한다. 90점 EVI 평균은 `/65535` 가정 0.511, `/10000` 민감도 가정 1.335으로 차이가 크다.

## 성능 저하 원인: 확인 사실과 가설

**확인 사실:** 표본은 단일 GEDI 파일의 90점이며 세 식생지수만 입력했다. Spatial 학습은 55점으로 줄고 평가 위치는 동쪽 끝에 집중된다. RF는 Spatial에서 Dummy보다 성능이 낮다. 반사율 변환은 미확정이다.

**가설:** 제한된 입력이 수종·수관 구조·지형 차이를 설명하지 못하거나, 적은 표본과 위치별 분포 차이로 RF가 과적합했을 수 있다. EVI의 반사율 배율 가정도 영향을 줄 수 있다. 각 요인의 기여도는 현재 결과만으로 분리할 수 없다.

## Week05 우선 작업

1. 원본 TIFF 생성 설정/evalscript 또는 검증된 FLOAT32 반사율로 배율·offset을 확인한다. 확인 전에는 기존 지수를 수정하지 않는다.
2. 버퍼를 둔 여러 공간 블록과 반복 분할로 성능 변동성을 평가하고, 각 fold에서 Dummy를 동일 학습점으로 비교한다.
3. 평가점의 AGBD·입력 지수·공간 분포 차이와 극단값의 영향도를 점검한다. 극단값을 임의 제거하지 않는다.
4. 검증 가능한 시공간 일치 자료와 수관 구조·지형 등 추가 변수를 검토하되, 독립 현장 자료가 없으면 현장 정확도로 해석하지 않는다.

## 산출물

`model_comparison.csv`, `validation_residuals.csv`, `residual_summary.csv`, `data_quality.json`과 PNG 그래프 5개. `data_quality.json`은 전체 컬럼 결측과 분할 검사 수치를 포함한다.
