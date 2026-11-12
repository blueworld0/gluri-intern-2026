내부 비교 23개: 새 방법 / 방법 A의 RMSE 137.72 / 162.47 Mg/ha, 평균 편향 -14.04 / -87.02 Mg/ha, 평균 추정÷평균 기준 0.948 / 0.676 (독립 검증 아님).

# NEON 조사구 CHM AGB 회귀 — 2026-11-12

## 평가 비교

| 평가 범위      | 방법                 |   n |   RMSE (Mg/ha) |   평균 편향 (Mg/ha) |   MAE (Mg/ha) |   평균 추정÷평균 기준 |   Pearson r |
|:---------------|:---------------------|----:|---------------:|--------------------:|--------------:|----------------------:|------------:|
| ORNL 평가      | 방법 A (기존 OOF)    |   8 |        111.369 |             -88.639 |        95.143 |                 0.601 |       0.441 |
| ORNL 평가      | 새 방법 (조사구 CHM) |   8 |         89.469 |              66.986 |        71.692 |                 1.301 |       0.547 |
| GRSM 평가      | 방법 A (기존 OOF)    |  15 |        184.014 |             -86.162 |       119.041 |                 0.707 |       0.504 |
| GRSM 평가      | 새 방법 (조사구 CHM) |  15 |        157.518 |             -57.248 |        80.824 |                 0.805 |       0.591 |
| 통합 23개 평가 | 방법 A (기존 OOF)    |  23 |        162.473 |             -87.023 |       110.729 |                 0.676 |       0.530 |
| 통합 23개 평가 | 새 방법 (조사구 CHM) |  23 |        137.717 |             -14.036 |        77.648 |                 0.948 |       0.462 |

단위는 지상부 AGB Mg/ha (= t/ha)입니다. 기준은 좌표 없는 나무까지 포함한 AGB_ref_t_ha입니다.
편향 = mean(추정−기준), 비율 = mean(추정)/mean(기준). 통합 23개는 지역이 평가일 때의 예측만 사용했습니다.
방법 A는 기존 바깥 GroupKFold OOF 예측이므로 새 방법의 지역 분리 평가와 평가 분할 방식이 다릅니다.
따라서 이 표는 같은 조사구에 대한 내부 비교이며, 두 방법의 동일한 지역 외 일반화 성능 비교는 아닙니다.

## 학습 지역 안 적합 (in-sample, 참고용)

아래 행은 새 모델이 학습에 사용한 지역에서의 적합도이며 비교 결론에 사용하지 않습니다.

| 평가 범위                                | 방법                 |   n |   RMSE (Mg/ha) |   평균 편향 (Mg/ha) |   MAE (Mg/ha) |   평균 추정÷평균 기준 |   Pearson r |
|:-----------------------------------------|:---------------------|----:|---------------:|--------------------:|--------------:|----------------------:|------------:|
| GRSM 학습 지역 안 적합 (in-sample, 참고) | 새 방법 (조사구 CHM) |  15 |        143.508 |              -0.000 |        92.152 |                 1.000 |       0.599 |
| ORNL 학습 지역 안 적합 (in-sample, 참고) | 새 방법 (조사구 CHM) |   8 |         59.244 |               0.000 |        53.205 |                 1.000 |       0.547 |

## 계수와 학습 R²

| direction       | train_site   | test_site   |          b0 |       b1 |         b2 |   train_n |   train_R2 |
|:----------------|:-------------|:------------|------------:|---------:|-----------:|----------:|-----------:|
| D1_GRSM_to_ORNL | GRSM         | ORNL        | -385.565565 | 9.118904 | 456.608151 |        15 |   0.359390 |
| D2_ORNL_to_GRSM | ORNL         | GRSM        |   18.746800 | 8.505411 |   0.000000 |         8 |   0.299187 |

고정식: AGB = b0 + b1 × mean_height_m + b2 × cover_2m (0~1).
절편 포함, 정규화·변환·규제·설정 탐색 없음. 평가 지역 정답은 학습에 사용하지 않았습니다.
음수 예측은 보정하지 않았고 predictions.csv의 negative_flag로 표시했습니다.

## 특징 및 입력 점검

- 경계와 CHM의 좌표계를 맞추고 픽셀 중심이 경계 안에 있는 1m 픽셀만 사용(all_touched=False).
- n_pixels는 결측을 포함한 전체 경계 내부 픽셀, n_valid_pixels는 유효 픽셀 수입니다.
- mean_height_m은 유효 픽셀의 평균(0m 포함); cover_2m은 유효 픽셀 중 CHM ≥ 2m 비율입니다.
- 결측률 >10%는 nodata_warning_gt10pct에만 표시하며, 조사구를 제외하지 않았습니다.
- 정수 CHM도 float32로 읽고 nodata=None이면 창 밖 채움값을 NaN으로 명시했습니다.

- 결측률 >10%: 0개; 음수 평가 예측: 0개; 제외: 0개.
- AGB_ref × 0.5 = C_ref 확인만: True; 최대 절대 차이 5.6843419e-14.
- 방법 A 입력: `/content/drive/MyDrive/neon_grsm_ornl/outputs/f55cc4636e94/methods/A_DF_nested_OOF/plot_comparison.csv` (SHA-256 `0109466981fd99a47e59d037713ac31b362e8ae0dbfb873a38011212924cfbff`).
- 경계 입력: `/content/drive/MyDrive/neon_grsm_ornl/outputs/f55cc4636e94/sampled_footprints.gpkg`, 레이어 `sampled_footprints` (SHA-256 `1ddc3e36a1cf1febe6eca3008f81567aea904a6418a1aecc9e4e232a9010d08c`).
- 경계 원본 CRS: `EPSG:32617`; 아래 각 CHM CRS에서 폴리곤 면적 계산.

| plotID   | chm_crs    |   resolution_m | dtype   |   source_nodata |   polygon_area_m2 | chm_path                                                                            |
|:---------|:-----------|---------------:|:--------|----------------:|------------------:|:------------------------------------------------------------------------------------|
| GRSM_001 | EPSG:32617 |              1 | float32 |           -9999 |           388.466 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_001_2021_chm_f55cc4636e94.tif |
| GRSM_009 | EPSG:32617 |              1 | float32 |           -9999 |           394.801 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_009_2021_chm_f55cc4636e94.tif |
| GRSM_012 | EPSG:32617 |              1 | float32 |           -9999 |           382.248 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_012_2021_chm_f55cc4636e94.tif |
| GRSM_014 | EPSG:32617 |              1 | float32 |           -9999 |           389.608 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_014_2021_chm_f55cc4636e94.tif |
| GRSM_016 | EPSG:32617 |              1 | float32 |           -9999 |           390.75  | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_016_2021_chm_f55cc4636e94.tif |
| GRSM_017 | EPSG:32617 |              1 | float32 |           -9999 |           392.465 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_017_2021_chm_f55cc4636e94.tif |
| GRSM_047 | EPSG:32617 |              1 | float32 |           -9999 |           800     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_047_2021_chm_f55cc4636e94.tif |
| GRSM_049 | EPSG:32617 |              1 | float32 |           -9999 |           790.75  | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_049_2021_chm_f55cc4636e94.tif |
| GRSM_050 | EPSG:32617 |              1 | float32 |           -9999 |           797.795 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_050_2021_chm_f55cc4636e94.tif |
| GRSM_051 | EPSG:32617 |              1 | float32 |           -9999 |           780.852 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_051_2021_chm_f55cc4636e94.tif |
| GRSM_055 | EPSG:32617 |              1 | float32 |           -9999 |           800     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_055_2021_chm_f55cc4636e94.tif |
| GRSM_057 | EPSG:32617 |              1 | float32 |           -9999 |           800     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_057_2021_chm_f55cc4636e94.tif |
| GRSM_058 | EPSG:32617 |              1 | float32 |           -9999 |           799.583 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_058_2021_chm_f55cc4636e94.tif |
| GRSM_059 | EPSG:32617 |              1 | float32 |           -9999 |           788.563 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_059_2021_chm_f55cc4636e94.tif |
| GRSM_060 | EPSG:32617 |              1 | float32 |           -9999 |           800     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/GRSM_060_2021_chm_f55cc4636e94.tif |
| ORNL_001 | EPSG:32616 |              1 | float32 |           -9999 |           400     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_001_2021_chm_f55cc4636e94.tif |
| ORNL_004 | EPSG:32616 |              1 | float32 |           -9999 |           381.942 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_004_2021_chm_f55cc4636e94.tif |
| ORNL_006 | EPSG:32616 |              1 | float32 |           -9999 |           390.326 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_006_2021_chm_f55cc4636e94.tif |
| ORNL_009 | EPSG:32616 |              1 | float32 |           -9999 |           400     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_009_2021_chm_f55cc4636e94.tif |
| ORNL_010 | EPSG:32616 |              1 | float32 |           -9999 |           393.498 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_010_2021_chm_f55cc4636e94.tif |
| ORNL_027 | EPSG:32616 |              1 | float32 |           -9999 |           390.269 | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_027_2021_chm_f55cc4636e94.tif |
| ORNL_049 | EPSG:32616 |              1 | float32 |           -9999 |           800     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_049_2021_chm_f55cc4636e94.tif |
| ORNL_057 | EPSG:32616 |              1 | float32 |           -9999 |           800     | /content/drive/MyDrive/neon_grsm_ornl/data/crops/ORNL_057_2021_chm_f55cc4636e94.tif |

실행 시각(UTC): 2026-11-12T04:46:27.239869+00:00

패키지 버전: `{"numpy": "2.4.6", "pandas": "2.2.3", "geopandas": "1.2.0", "rasterio": "1.5.2", "scikit-learn": "1.6.1", "matplotlib": "3.10.0"}`

## 해석상 한계

- 23개는 앞선 실험에서 이미 본 자료로 내부 비교이며 독립 검증이 아닙니다.
- 방법 A의 설정(score 0.2, patch 600)도 같은 23개에서 선택되었습니다. 여기서는 배포 설정으로 다시 예측하지 않고 기존 바깥 GroupKFold OOF 열을 사용했습니다.
- 학습 표본은 8~15개로 작고 지역은 2곳뿐입니다.
- NEON LiDAR CHM과 국내 Meta 추정 높이는 다른 자료이므로 국내 성능으로 해석할 수 없습니다.
- 평균 추정÷평균 기준이 1에 가까워도 조사구별 오차와 RMSE·MAE·편향·상관을 함께 봐야 합니다.
- 방법 A와 새 방법의 평가 분할이 다릅니다. 비교표만으로 우열이나 지역 외 성능을 확정할 수 없습니다.
- 기준 AGB도 기존 실측 기반 공식의 산출값이며 벌채 실측의 독립 정답은 아닙니다.

## 다음 단계 판단 (직접 작성)

- 추가 검증 진행 여부: ____________________
- 판단 근거: ____________________
- 다음 자료·실험: ____________________
