# Random Forest와 학습 평균 Dummy 비교

## 방법과 재현 확인

- 입력: Week03 SCL 통과 90개 점의 `ndvi`, `evi`, `ndmi`; 목표: GEDI L4A 추정 `agbd`(Mg/ha).
- 기존 `train_baseline.py`의 `SEED=42`, 무작위 80/20 분할과 `spatial_split` 함수를 그대로 사용했다. 각 검증에서 두 모델에 동일한 학습·검증 **인덱스와 순서**를 전달했다.
- 공간 분할의 500m 버퍼 제외점은 두 모델 모두 학습·검증에 사용하지 않았다. Dummy 평균은 전체 90개가 아닌 해당 **학습 집합의 AGBD 평균**으로 학습했다.
- 실행 환경: scikit-learn 1.9.1; Random Forest는 `n_estimators=100, random_state=42, n_jobs=1`, Dummy는 `strategy='mean'`.
- 기존 `baseline_model/model_metrics.csv` 및 `validation_predictions.csv`와 RF의 검증점 ID·순서·실제값·예측값·지표를 대조했다.
- **기존 RF 결과 재현: 성공** (지표·예측값 허용오차 1e-8). 기존 결과 파일은 변경하지 않았다.

## 검증 데이터 성능

| 분할 | 모델 | 학습 | 검증 | MAE (Mg/ha) | RMSE (Mg/ha) | R² |
|---|---|---:|---:|---:|---:|---:|
| random_80_20 | RandomForest | 72 | 18 | 42.87 | 52.34 | -0.139 |
| random_80_20 | DummyMean | 72 | 18 | 43.34 | 50.94 | -0.080 |
| spatial_block_500m_buffer | RandomForest | 55 | 20 | 57.71 | 70.61 | -2.165 |
| spatial_block_500m_buffer | DummyMean | 55 | 20 | 48.07 | 57.16 | -1.074 |

## 해석

**결론:** 무작위 분할에서는 RF의 MAE만 근소하게 낮고 RMSE·R²는 Dummy보다 나쁘다. 공간 분리에서는 RF가 세 지표 모두 Dummy보다 나쁘다. 현재 표본과 단일 분할에서는 RF가 단순 학습 평균보다 일관되게 낫다는 근거가 없다.

- **random_80_20:** RF − Dummy MAE -0.47 Mg/ha, RMSE +1.39 Mg/ha, R² -0.060. 학습 AGBD 평균 106.71 Mg/ha; 검증 AGBD 평균 92.88 Mg/ha. 버퍼 제외 0개.
- **spatial_block_500m_buffer:** RF − Dummy MAE +9.64 Mg/ha, RMSE +13.45 Mg/ha, R² -1.091. 학습 AGBD 평균 89.79 Mg/ha; 검증 AGBD 평균 130.92 Mg/ha. 버퍼 제외 15개.
- Dummy의 R²가 반드시 0인 것은 아니다. `strategy='mean'`은 **학습** 평균을 예측하지만 R²의 분모 기준은 **검증** 집합 평균이기 때문이다.
- 이 한 번의 분할 결과는 통계적 유의성 검정이 아니다. RF와 Dummy의 우열을 다른 위치·기간에 일반화할 수 없다.
- 가능한 원인: 90개 중 한 GEDI 궤도의 점만 사용해 공간 대표성이 낮고, NDVI·EVI·NDMI 세 변수만으로 수고·수종·밀도·지형 차이를 담기 어렵다. 공간 홀드아웃에서는 학습 55개와 동쪽 검증 20개의 분포가 다를 수 있다. 반사도 `/65535` 변환은 미확정이며 특히 EVI에 영향을 줄 수 있다. 이 설명들은 진단 가능한 가설이지 인과적 결론은 아니다.
- GEDI L4A AGBD는 기존 모델의 추정값이며 현장 실측값이 아니다. 이 점수는 독립 현장 검증 성능이나 금정산 전체 정확도가 아니다.
- CSV에는 네 실험의 검증 점수와 분할 크기, 최소 학습·검증 거리, 버퍼 제외 수가 있다. 그래프 위쪽은 **동일 축의 검증점 실제값·두 모델 예측값**, 아래쪽은 지표 비교다.
