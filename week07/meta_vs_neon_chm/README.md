# NEON LiDAR CHM vs Meta/WRI v1 높이 지도 — 7주차 작은 비교 실험

- 작성일: 2026-11-13
- 노트북: [Meta_vs_NEON_CHM.ipynb](Meta_vs_NEON_CHM.ipynb)
- 작성자: 이연지
- 상태: Colab 실행 완료 (2026-11-13, 실행 약 45초, Meta 원격 전송 23.1 MiB). 결과는 `outputs/`
- 결과 보고서: [20261113_Meta_NEON_높이비교_결과.md](../docs/20261113_Meta_NEON_높이비교_결과.md)

## 확인할 질문

Meta/WRI v1 높이 지도는 미국 NEON LiDAR 높이 대비 얼마나 틀리는지 확인함. 이 높이 입력 차이가 6주차 면적 회귀 AGB 추정에 미치는 영향도 확인함.

GRSM 15개와 ORNL 8개, 총 23개 조사구에서 조사구 평균 높이·2m 이상 피복률·p95 높이를 주 결과로 비교함. 양쪽 모두 유효한 픽셀 쌍 비교는 보조 결과로 제시함. AGB는 6주차 계수를 그대로 적용하며 재학습·변수 추가·보정하지 않음.

## 입력자료

NEON Drive 루트는 `/content/drive/MyDrive/neon_grsm_ornl/`임. 기존 원본·결과는 읽기만 함. 새로 읽는 외부 자료는 아래 Meta **v1** COG 부분 범위와 선택적 관측일 메타데이터뿐임.

| 입력 | 위치·버전 | 확인 내용 |
| --- | --- | --- |
| NEON LiDAR CHM | Drive `data/crops/<plotID>_2021_chm_f55cc4636e94.tif` | 2021년·1m·단일 밴드·EPSG:32617 또는 32616 확인함 |
| 조사구 경계 | Drive `outputs/f55cc4636e94/sampled_footprints.gpkg` | 6주차의 23개 실제 조사구 경계를 사용함 |
| 기준 AGB | Drive `outputs/f55cc4636e94/methods/A_DF_nested_OOF/plot_comparison.csv` | 좌표 없는 나무도 포함한 `AGB_ref_t_ha`를 사용함. 단위는 지상부 Mg/ha (= t/ha)임 |
| 6주차 특징 | [features.csv](../../week06/neon_area_chm/outputs/features.csv) | 원문 23행을 노트북에 포함했음. mean_height_m·cover_2m 재현 확인용임 |
| 6주차 고정 계수 | [coefficients.csv](../../week06/neon_area_chm/outputs/coefficients.csv) | 원문 2행을 노트북에 포함했음. 별도 Drive 업로드는 불필요함 |
| Meta/WRI Global Canopy Height v1 | `https://dataforgood-fb-data.s3.amazonaws.com/forests/v1/alsgedi_global_v6_float/chm/<quadkey>.tif` | 조사구 경계+30m 버퍼가 걸치는 zoom9 quadkey의 필요한 window만 읽음 |
| Meta 관측일 메타데이터 | 같은 v1 경로의 `tiles.geojson`, `chm/tiles.geojson`, 목록에서 찾은 관련 JSON/GeoJSON | 관측·영상 획득을 명시한 필드만 사용함. 없거나 접근·해석 실패 시 `미확인`으로 남김 |

6주차 CSV 원문의 SHA-256을 노트북과 `meta_sources.json`에 기록함. 고정 계수는 다음과 같음.

| 방향 | 적용할 평가 지역 | b0 | b1 (평균 높이) | b2 (피복률) |
| --- | --- | ---: | ---: | ---: |
| D1_GRSM_to_ORNL | ORNL 8개 | -385.5655652750414 | 9.11890418658884 | 456.60815051757095 |
| D2_ORNL_to_GRSM | GRSM 15개 | 18.746800248381135 | 8.505411493461185 | 0.0 |

quadkey 함수 3개는 [geumjeongsan_core.py](../../week06/geumjeongsan_deepforest/notebooks/geumjeongsan_core.py)에서 그대로 가져왔음. 원격 COG window 읽기·원래 격자에서의 병합·tags 기록 방식은 [geumjeongsan_runtime.py](../../week06/geumjeongsan_deepforest/notebooks/geumjeongsan_runtime.py)를 재사용했음. NEON 픽셀 읽기·결측 처리·저장 보호 구조는 [6주차 노트북](../../week06/neon_area_chm/NEON_area_CHM_regression.ipynb)을 재사용했음.

## 실행 환경과 순서

Google Colab CPU 런타임에서 위에서 아래로 실행함. GPU·NEON 토큰은 필요하지 않음. 코드 셀 12개이며 각 셀 위에 한국어 설명을 둠.

1. **설치 → Drive 마운트:** 패키지를 설치하고 입력·출력 경로를 지정함. 기존 출력 폴더나 ZIP이 있으면 중단함.
2. **입력 검증:** 6주차와 동일한 23개 ID·GRSM 15/ORNL 8·파일 존재·좌표계·1m 해상도를 확인함. 기준값과 내장 CSV도 검증함.
3. **NEON 재현 확인:** 6주차와 같은 float32 읽기·픽셀 중심 포함(`all_touched=False`)·0m와 결측 구분 규칙을 적용함. mean_height_m·cover_2m이 6주차 값과 절대오차 `1e-6` 이내인지 검사함(`rtol=0`). 다르면 Meta 읽기 전에 전체 중단함. p95는 유효 높이의 linear 백분위수로 추가 계산함.
4. **Meta 읽기:** 실제 조사구 30m 버퍼와 교차하는 zoom9 quadkey만 선택함. 관측일 메타데이터를 선택적으로 확인한 뒤 필요한 window를 읽음. strip TIFF도 허용함. 먼저 23개 조사구 모두에 대해 `겹치는 블록 수 × block_h × block_w × 원본 dtype 바이트 수`를 계산하고, 조사구에 걸친 모든 타일의 값을 합산함. 조사구당 **300 MB(300,000,000바이트, 압축 전)** 를 초과하면 이유를 출력하고 픽셀 읽기 전에 전체 중단함. 실행 시 조사구별 예상 MB와 통과 수를 먼저 표시함. Range 없는 GET·206이 아닌 응답·전체 타일 window는 계속 거절함. 전체 다운로드로 자동 전환하지 않음.
5. **정렬:** 작은 원본 창을 같은 native 격자에서 병합한 뒤 `bilinear` 한 번으로 NEON CRS·1m 격자에 재투영함. 목적 격자는 NEON 원래 격자에 정렬한 30m 버퍼 window임. 원래 조사구 window를 잘라 동일 경계에서 비교함. nodata는 명시적으로 NaN으로 처리함.
6. **조사구 비교:** 각 자료 자체의 유효 픽셀로 평균 높이·피복률·p95를 계산함. 특징마다 GRSM·ORNL·통합 23개의 Meta−NEON 편향·RMSE·MAE·r·평균비를 계산함.
7. **픽셀 비교:** 두 자료가 모두 유효한 쌍으로 사이트별·통합 편향·RMSE·r을 계산함. NEON 기준 [0,10), [10,20), [20,30), [30,∞)m 구간별 편향도 계산함.
8. **AGB 영향:** GRSM 조사구는 D2, ORNL 조사구는 D1의 고정 계수만 적용함. NEON 입력·Meta 입력 각각을 같은 기준 AGB와 비교해 RMSE·편향·평균비·r을 계산함.
9. **그림 → 저장:** 영어 라벨 그림 5개·표·출처·보고서를 저장하고 ZIP을 내려받음. 마지막에 조사구 평균 높이 표와 AGB 비교표를 출력함.

평균 높이에는 유효한 0m를 포함함. 결측률은 경계 안 전체 픽셀을 분모로 계산하고, 높이·피복률·p95는 각 자료의 유효 픽셀을 사용함. 결측률 >10%는 경고 열만 남기며 조사구를 제외하지 않음. 어느 자료든 유효 픽셀이 0개이거나 유효 픽셀 쌍이 0개인 조사구가 있으면 전체 중단함. 음수 높이·음수 AGB 예측은 임의 보정하지 않고 기록함. 음수 NEON 픽셀이 있으면 전체 픽셀 지표에 유지하고 높이 구간 표에 별도 점검 행을 추가함.

평균비는 `mean(Meta)/mean(NEON)`, AGB에서는 `mean(예측)/mean(기준)`임. AGB 편향은 예측−기준임. 평균 기준이 0이면 비율, 상수 벡터이면 Pearson r을 NaN으로 기록함. 피복률 오차 단위는 0~1 비율이며 높이 오차는 m, AGB 오차는 Mg/ha임.

실행 전 다음 사항을 확인함.

- Drive에 위 NEON 원본 23개와 경계·기준 CSV가 실제로 있는지 확인함. 6주차 CSV는 내장되어 별도 복사하지 않아도 됨.
- Drive `outputs/meta_vs_neon_20261113/`와 Colab `/content/meta_vs_neon_results.zip`이 없는지 확인함. 이미 있으면 이전 결과를 별도 위치에 보관한 후 실행함.
- 예상 실행 시간은 설치 포함 **약 5~20분**으로 잡음. 미실행 상태의 계획용 추정이며 서버·Drive 응답 속도에 따라 달라짐.
- Meta TIFF 원격 읽기는 **약 10~200MiB**를 계획용 범위로 잡았음. 실제 측정값은 아님. 작은 window도 내부 블록 전체를 읽으므로 픽셀 배열 크기보다 전송량이 커질 수 있음. 실제 TIFF 응답 본문 전송량은 실행 중 출력하고 `meta_sources.json`에 기록함. 이 수치에는 중복 Range와 TIFF 헤더가 포함되고 선택적 JSON/XML은 제외됨. 선택적 메타데이터는 파일당 8MiB 이하로 제한하며 한도 초과 시 미확인으로 남김. [Rasterio window 읽기 설명](https://rasterio.readthedocs.io/en/stable/topics/windowed-rw.html)을 참고함.

로컬 Windows 검증은 표준 라이브러리로 노트북 JSON·전체 코드 셀 문법·CSV 23개 ID와 2개 계수 방향·내장 원문 해시 등을 확인함. pandas·rasterio·geopandas로 실제 지리 처리나 회귀를 실행하지 않았음. 실제 원격 window 읽기·정렬·특징 재현은 Colab에서 검증함. 재투영은 [Rasterio 공식 설명](https://rasterio.readthedocs.io/en/stable/topics/reproject.html)을 참고함.

## 주요 결과

NEON 재현 확인 통과(6주차 features와 일치), 23개 조사구 모두 결측 0, 제외 조사구 0개임.

**조사구 평균 높이 (Meta − NEON)**

| 범위 | 편향 (m) | RMSE (m) | r | Meta ÷ NEON |
| --- | ---: | ---: | ---: | ---: |
| GRSM 15 | −8.1 | 10.2 | 0.83 | 0.69 |
| ORNL 8 | −5.4 | 6.3 | 0.69 | 0.78 |
| 통합 23 | **−7.1** | **9.1** | 0.78 | **0.72** |

- p95 높이는 통합 편향 −10.1m, 비율 0.67로 더 낮게 나옴. 2m 이상 피복률은 둘 다 거의 1.0이라 차이 판단 불가
- 픽셀 1.3만 쌍 기준 편향 −7.2m, RMSE 10.1m, r 0.74

**NEON 높이 구간별 편향 (픽셀, 통합)**

| NEON 높이 | 0–10m | 10–20m | 20–30m | 30m 이상 |
| --- | ---: | ---: | ---: | ---: |
| Meta 편향 (m) | +3.7 | **−1.2** | −6.9 | **−14.4** |

- 20m까지는 잘 맞고, 그 위에서는 Meta 값이 20~25m에서 더 올라가지 않음(포화). 수관 단위 모양 없이 평활화됨

**6주차 고정 계수 AGB (Mg/ha)**

| 범위 | NEON 입력 RMSE | Meta 입력 RMSE | Meta 입력 편향 | Meta 입력 추정÷기준 |
| --- | ---: | ---: | ---: | ---: |
| GRSM 15 | 158 | 199 | −126 | 0.57 |
| ORNL 8 | 89 | 61 | +18 | 1.08 |
| 통합 23 | **138** | **165** | **−76** | **0.72** |

- 통합 기준으로 Meta 입력 시 오차가 커지고 AGB를 약 28% 적게 추정함
- ORNL 개선은 원래 모델의 과대추정(+67)을 Meta의 낮은 높이가 상쇄한 결과임. Meta가 정확해서가 아님

**Meta 원본 확인:** EPSG:3857, 픽셀 1.19m, uint8(1m 단위 정수), 단위 METERS 확인함. 관측일은 `metadata/<quadkey>.geojson`이 있으나 8MiB 한도 초과로 읽지 않아 미확인임

| 결과 | 기록 내용 |
| --- | --- |
| `plot_features.csv` | 조사구·지역·면적·NEON/Meta 평균 높이·피복률·p95·픽셀 수·결측률·0m 수·음수 수·결측 경고·Meta 관측일·단위 경고·quadkey·평가 방향·고정 계수 AGB 예측과 오차를 기록함 |
| `plot_metrics.csv` | 특징 3개 × GRSM/ORNL/통합 23개에서 Meta−NEON 편향·RMSE·MAE·Pearson r·평균비를 기록함 (9행) |
| `pixel_metrics.csv` | 양쪽 유효 픽셀 쌍의 사이트별·통합 n_pairs·n_plots·편향·RMSE·r을 기록함 (3행) |
| `height_bin_bias.csv` | NEON 기준 높이 구간별 유효 쌍 수와 Meta−NEON 평균 편향을 기록함. 빈 구간은 n=0·NaN으로 남김 |
| `agb_comparison.csv` | 지역별·통합 고정 계수 AGB의 NEON 입력 vs Meta 입력 지표를 기록함 (6행) |
| `meta_sources.json` | quadkey·URL·원본 CRS·픽셀 간격과 좌표 단위·dtype·nodata·tags·높이 단위 확인 상태·is_tiled·block_shapes·compression·width·height·읽은 window·블록 수·압축 전 예상 바이트/MB·조사구별 합계 MB와 300 MB 통과 여부·실제 Range 전송량·관측일 조회 상태·고정 계수와 원문 해시를 기록함 |
| `report.md` | 수치 결과·비교표·관측일 및 단위 점검·필수 한계 문구·다음 단계 판단 빈칸을 기록함 |

영어 라벨 그림을 `figures/`에 저장함.

- `scatter_plot_mean_height.png`: NEON(x) vs Meta(y) 조사구 평균 높이·1:1선·사이트별 색을 표시함.
- `pixel_hexbin.png`: 모든 유효 픽셀 쌍의 NEON(x) vs Meta(y)·1:1선을 표시함.
- `bias_by_height_bin.png`: NEON 기준 높이 구간별 편향을 사이트별·통합으로 표시함.
- `agb_neon_vs_meta.png`: 기준 AGB(x) vs 고정 계수 예측(y), NEON 입력=원·Meta 입력=삼각형으로 표시함.
- `examples.png`: **부호 있는 조사구 평균 높이 편향(Meta−NEON)** 이 가장 작은 곳과 가장 큰 곳 2개를 선택함. NEON·Meta·차이 지도를 나란히 표시함. 절대오차 순위와는 다름.

Drive `/content/drive/MyDrive/neon_grsm_ornl/outputs/meta_vs_neon_20261113/`에 표 5개·출처 JSON·보고서·그림 5개와 ZIP 사본을 저장함. `/content/meta_vs_neon_results.zip`에는 산출물 12개를 포함함. 큰 원본 자료는 저장소에 넣지 않음.

팀 Drive 결과: [02_분석결과](https://drive.google.com/drive/folders/13dHSpKDKMAVI7MbmWNuPDVocu5R5LGR2) / `20261113_Meta_NEON높이비교_v01`

## 확인한 내용과 남은 한계

- 시점 차이: NEON LiDAR는 2021년, Meta 영상은 관측일 기록값(없으면 미확인, 명목상 2018–2020 전후). 차이에 실제 생장·교란이 섞일 수 있음
- 학습 자료 겹침: Meta 모델은 NEON 라이다로 학습했으므로 GRSM·ORNL이 학습에 쓰였을 수 있음. 그러면 결과가 실제보다 좋게 나옴. 국내(금정산) 정확도는 이보다 나쁠 수 있음
- 23개 조사구의 내부 비교이며 독립 검증 아님. 미국 동부 활엽수림 결과라 국내 숲에 그대로 적용하지 않음
- bilinear 재투영의 평활화·위치 정합·시점·센서 차이가 함께 섞임. 이번 비교만으로 순수한 모델 오차를 분리할 수 없음
- 픽셀은 공간적으로 상관되어 있어 유효 쌍 수를 독립 표본 수로 해석하지 않음. 통합 픽셀 지표는 픽셀 수로 가중되고 조사구 지표는 조사구마다 같은 가중치를 가짐
- 6주차 계수는 8~15개 학습 표본과 2개 지역에 한정됨. AGB 영향은 고정 모델의 입력 민감도 비교이며 Meta로 새 회귀를 학습한 성능이 아님
- 원본 tags에서 높이 단위 m을 확인하지 못하면 경고를 기록하고 m로 가정해 계속함. 날짜를 확인하지 못하면 `미확인`으로 남김. 파일 생성·수정일은 관측일로 대체하지 않음
- 평균비가 1에 가까워도 조사구별 오차·RMSE·편향·상관을 함께 봐야 함
- 확인: Meta v1은 NEON 대비 높이를 평균 7.1m(28%) 낮게 잡았고, 20m 이상에서 과소 추정이 커짐. 6주차 회귀에 그대로 넣으면 AGB를 과소 추정함
- 다음(후속 범위): Meta 관측일(`metadata/<quadkey>.geojson`) 확인, 금정산 조사구 높이 분포 확인 후 20m 이하 범위인지 판단, 높이 구간별 보정 가능성 검토
