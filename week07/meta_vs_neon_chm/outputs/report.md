내부 비교 23개: Meta−NEON 조사구 평균 높이 편향 -7.12m, RMSE 9.05m; 고정 계수 AGB RMSE NEON 137.72 → Meta 164.86 Mg/ha임.

# Meta v1 vs NEON CHM — 2026-11-13

## 조사구 특징 비교 (주 결과)

| scope   | feature       | unit     |   n |     bias |    RMSE |     MAE |   pearson_r |   mean_ratio |
|:--------|:--------------|:---------|----:|---------:|--------:|--------:|------------:|-------------:|
| GRSM    | mean_height_m | m        |  15 |  -8.0625 | 10.2235 |  8.4217 |      0.8347 |       0.6852 |
| GRSM    | cover_2m      | fraction |  15 |   0.0186 |  0.0444 |  0.0186 |      0.9697 |       1.0191 |
| GRSM    | p95_height_m  | m        |  15 | -10.9330 | 12.5564 | 10.9330 |      0.7665 |       0.6574 |
| ORNL    | mean_height_m | m        |   8 |  -5.3529 |  6.2848 |  5.4911 |      0.6914 |       0.7765 |
| ORNL    | cover_2m      | fraction |   8 |   0.0000 |  0.0000 |  0.0000 |    nan      |       1.0000 |
| ORNL    | p95_height_m  | m        |   8 |  -8.4438 |  9.4305 |  8.4438 |      0.3477 |       0.7043 |
| ALL_23  | mean_height_m | m        |  23 |  -7.1200 |  9.0501 |  7.4024 |      0.7796 |       0.7156 |
| ALL_23  | cover_2m      | fraction |  23 |   0.0121 |  0.0359 |  0.0121 |      0.9682 |       1.0123 |
| ALL_23  | p95_height_m  | m        |  23 | -10.0672 | 11.5653 | 10.0672 |      0.6925 |       0.6726 |

## 픽셀 쌍 비교 (보조 결과)

| scope   |   n_plots |   n_pairs |   bias_m |   RMSE_m |   pearson_r |
|:--------|----------:|----------:|---------:|---------:|------------:|
| GRSM    |        15 |      9506 |  -7.7828 |  11.1087 |      0.7756 |
| ORNL    |         8 |      3954 |  -5.7124 |   7.1160 |      0.6113 |
| ALL_23  |        23 |     13460 |  -7.1746 |  10.1009 |      0.7432 |

## NEON 높이 구간별 편향

| scope   | neon_height_bin_m   |   n_pairs |   bias_m |
|:--------|:--------------------|----------:|---------:|
| GRSM    | 0-10                |      1420 |   3.6509 |
| GRSM    | 10-20               |      1548 |  -1.1598 |
| GRSM    | 20-30               |      2590 |  -7.3730 |
| GRSM    | 30+                 |      3948 | -14.7609 |
| ORNL    | 0-10                |        73 |   5.6531 |
| ORNL    | 10-20               |       768 |  -1.2138 |
| ORNL    | 20-30               |      2675 |  -6.4898 |
| ORNL    | 30+                 |       438 | -10.7464 |
| ALL_23  | 0-10                |      1493 |   3.7488 |
| ALL_23  | 10-20               |      2316 |  -1.1777 |
| ALL_23  | 20-30               |      5265 |  -6.9243 |
| ALL_23  | 30+                 |      4386 | -14.3600 |

## 고정 계수 AGB 비교

| scope   | input_source   |   n |   RMSE_Mg_ha |   bias_Mg_ha |   mean_ratio |   pearson_r |
|:--------|:---------------|----:|-------------:|-------------:|-------------:|------------:|
| GRSM    | NEON           |  15 |     157.5182 |     -57.2477 |       0.8052 |      0.5915 |
| GRSM    | META           |  15 |     199.2864 |    -125.8225 |       0.5718 |      0.7133 |
| ORNL    | NEON           |   8 |      89.4692 |      66.9861 |       1.3012 |      0.5470 |
| ORNL    | META           |   8 |      60.5950 |      18.1737 |       1.0817 |      0.6021 |
| ALL_23  | NEON           |  23 |     137.7170 |     -14.0359 |       0.9478 |      0.4625 |
| ALL_23  | META           |  23 |     164.8583 |     -75.7368 |       0.7185 |      0.3187 |

편향은 Meta−NEON, AGB에서는 예측−기준으로 정의했음. 비율은 평균 추정÷평균 기준임.
조사구 특징은 각각의 유효 픽셀로 계산했으며 픽셀 비교는 양쪽 모두 유효한 쌍만 사용했음.
통합 픽셀 지표는 모든 쌍을 합친 픽셀 가중 결과이며, 조사구 지표는 조사구마다 같은 가중치임.
0m는 유효값으로 포함했음. 결측률 10% 초과도 경고만 남겼으며 조사구는 제외하지 않았음.
분모가 0인 비율과 상수 벡터의 r은 NaN으로 기록했음. 음수 높이·AGB는 임의 보정하지 않았음.
높이 구간은 NEON 기준 [0,10), [10,20), [20,30), [30,∞)m이며 음수가 있으면 별도 점검 구간을 기록했음.
AGB 단위는 지상부 Mg/ha (= t/ha)이며 탄소가 아님. 기준은 좌표 없는 나무도 포함했음.
6주차 계수를 그대로 썼으며 GRSM은 D2, ORNL은 D1 평가 방향만 적용했음. 재학습·보정하지 않았음.
예시 지도는 부호 있는 조사구 평균 높이 편향의 최솟값·최댓값 조사구를 선택했음.

## 관측일 및 단위 확인

| plotID   | siteID   | meta_observation_date   | meta_unit_warning   |   neon_nodata_frac |   meta_nodata_frac | neon_nodata_warning_gt10pct   | meta_nodata_warning_gt10pct   |
|:---------|:---------|:------------------------|:--------------------|-------------------:|-------------------:|:------------------------------|:------------------------------|
| GRSM_001 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_009 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_012 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_014 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_016 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_017 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_047 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_049 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_050 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_051 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_055 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_057 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_058 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_059 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| GRSM_060 | GRSM     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_001 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_004 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_006 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_009 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_010 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_027 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_049 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |
| ORNL_057 | ORNL     | 미확인                  | False               |                  0 |                  0 | False                         | False                         |

원본 CRS·픽셀 간격·dtype·nodata·tags·unit·읽은 window는 meta_sources.json에 기록했음. 높이 단위가 원본 tags에서 확인되지 않으면 경고를 남기고 m로 가정했으며 단위 변환은 하지 않았음. 파일 생성·수정일은 관측일로 쓰지 않았음. 관측일 미확인은 영상 시점이 일치함을 뜻하지 않음.

실제 TIFF 전송량: 23.06 MiB (반복 Range 포함, 선택적 JSON/XML 제외)였음.

## 6주차 고정 계수

| direction       | train_site   | test_site   |            b0 |         b1 |           b2 |   train_n |   train_R2 |
|:----------------|:-------------|:------------|--------------:|-----------:|-------------:|----------:|-----------:|
| D1_GRSM_to_ORNL | GRSM         | ORNL        | -385.56556528 | 9.11890419 | 456.60815052 |        15 | 0.35939010 |
| D2_ORNL_to_GRSM | ORNL         | GRSM        |   18.74680025 | 8.50541149 |   0.00000000 |         8 | 0.29918662 |

원본 coefficients.csv SHA-256: `0a5746926b1a90cc35c7bf2af09af9d93e1f4d2ea679c8d2548823b817fc3b28`

원본 features.csv SHA-256: `41ff333270d5b5acf0e086106cf3080f8b471d6434daaebf0a97310d8a714c41`

기준 CSV: `/content/drive/MyDrive/neon_grsm_ornl/outputs/f55cc4636e94/methods/A_DF_nested_OOF/plot_comparison.csv`; 경계: `/content/drive/MyDrive/neon_grsm_ornl/outputs/f55cc4636e94/sampled_footprints.gpkg`; 실행 시각 UTC: 2026-11-13T04:33:00.297743+00:00

## 확인한 내용과 남은 한계

- 시점 차이: NEON LiDAR는 2021년, Meta 영상은 관측일 기록값(없으면 미확인, 명목상 2018–2020 전후). 차이에 실제 생장·교란이 섞일 수 있음
- 학습 자료 겹침: Meta 모델은 NEON 라이다로 학습했으므로 GRSM·ORNL이 학습에 쓰였을 수 있음. 그러면 결과가 실제보다 좋게 나옴. 국내(금정산) 정확도는 이보다 나쁠 수 있음
- 23개 조사구의 내부 비교이며 독립 검증 아님. 미국 동부 활엽수림 결과라 국내 숲에 그대로 적용하지 않음
- bilinear 재투영의 평활화·위치 정합 오차와 원자료의 시점·센서 차이가 함께 섞임. 이번 비교만으로 순수한 모델 오차를 분리할 수 없음
- 픽셀은 공간적으로 상관되어 있어 유효 픽셀 수를 독립 표본 수로 해석하지 않음
- 6주차 계수도 8~15개 학습 표본과 2개 지역에 한정됨. AGB 영향은 고정 모델의 입력 민감도 비교이며 Meta로 새 모델을 학습한 성능이 아님
- 평균비가 1에 가까워도 조사구별 오차·RMSE·편향·상관을 함께 봐야 함
- Meta 높이 단위·관측일이 미확인이면 해당 가정이 남음. 원본 tags와 관측일 기록을 먼저 확인해야 함

## 다음 단계 판단 (직접 채움)

- 국내 적용 판단: ____________________
- 근거 및 추가 검증: ____________________
