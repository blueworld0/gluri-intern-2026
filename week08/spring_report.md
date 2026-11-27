# Week08 후속 분석 — 봄철 검증 후보와 고도·경사 비교

학습·검증 분포와 구간별 예측은 [고 AGBD 진단](high_agbd_report.md)에 정리했다.

## 2020년 5월 봄철 검증 결과

2025년 4월 산림 184점으로 학습하고 2020년 5월 10일 O07991의 산림 339점에서 검증했다. **10밴드에 고도·경사를 추가한 랜덤 포레스트의 RMSE는 63.87 Mg/ha로, 학습 평균 66.47 Mg/ha보다 2.60 Mg/ha 낮았다.** 10밴드만 사용한 모델은 66.50 Mg/ha였다. [예측값](outputs/multiyear_spring_O07991_predictions.csv) · [상세 지표](outputs/multiyear_spring_O07991_summary.json)

| 범위와 지표 | 학습 평균 | RF 10밴드 | RF 10밴드 + 고도·경사 |
| --- | ---: | ---: | ---: |
| 전체 339점 RMSE (Mg/ha) | 66.47 | 66.50 | 63.87 |
| 전체 339점 MAE (Mg/ha) | 50.92 | 50.26 | 49.58 |
| 전체 339점 R² | −0.081 | −0.082 | 0.002 |
| 학습점에서 500m 이상 떨어진 282점 RMSE (Mg/ha) | 65.55 | 65.93 | 63.89 |
| 고 AGBD 129점 평균 과소 예측 (Mg/ha) | 85.10 | 82.61 | 74.43 |

NASA 원본의 크기와 SHA-512를 확인했다. 구역 내 원본 4,286점 중 GEDI 품질 기준을 923점이 통과했고, EGIS 산림 점은 375점(고 AGBD 147점)이었다. 2020년 5월 4일 T52SDD 재처리 Sentinel-2 제품은 SCL 식생 기준을 통과한 점이 343점으로 가장 많았다. 이후 12.5m 원의 SCL·CLASSI·10밴드 원본 품질 조건을 적용해 339점(고 AGBD 129점)을 평가했다. [GEDI 품질 집계](outputs/multiyear_spring_gedi_quality_audit.json) · [영상 선택](outputs/multiyear_spring_O07991_scl_counts.csv) · [원본 영상 품질 기록](outputs/multiyear_spring_O07991_footprint_audit.csv)

고도·경사 추가 후에도 고 AGBD 129점은 모두 과소 예측됐다. 이 모델의 평균 과소 예측량은 10밴드 모델보다 8.18 Mg/ha 작았다. 두 궤도의 GEDI ID는 겹치지 않으며, 평가점 57점은 학습점과 500m 이내였다. 500m 이상 떨어진 평가점까지 같은 방향으로 개선됐지만, 이번 결과는 단일 과거 궤도 평가에 해당하며 시점에 따른 산림 변화와 모델 일반화 효과를 따로 추정할 수 없다.

## 다른 해의 4~5월 후보

동일 연구 구역에서 2019~2026년 4~5월 NASA CMR 목록을 확인했다. 2025년 외에는 2019년 1개, 2020년 7개, 2022년 2개, 2024년 4개 등 14개 관측 파일이 조회됐다. 2021·2023·2026년에는 이 기간과 구역에 해당하는 V003 목록이 조회되지 않았다. [다년 GEDI 목록](outputs/multiyear_spring_gedi_inventory.csv) · [Sentinel-2 후보](outputs/multiyear_spring_sentinel_candidates.csv)

GFZ 공개 GEDI V003에서 조회되는 표본을 대상으로 사용 가능한 GEDI 품질 플래그, EGIS 산림 기준, 168시간 이내 Sentinel-2 재처리 제품의 중심 SCL 식생 기준을 적용했다. 사전 점검상 **2020년 5월 10일 O07991**이 가장 유력했다. 5월 4일 T52SDD 재처리 제품에서 677점이 SCL을 통과했고, 그중 318점이 고 AGBD였다. [GFZ 공개 데이터베이스](https://gedidb.readthedocs.io/en/latest/user/tiledb_database.html) · [사전 점검표](outputs/multiyear_spring_gfz_precheck.csv)

| GEDI 날짜 (UTC) | 궤도 | GFZ 산림 후보 | SCL 통과 | 그중 고 AGBD |
| --- | --- | ---: | ---: | ---: |
| 2019-05-12 | O02345 | 46 | 39 | 0 |
| 2020-04-05 | O07442 | 184 | 117 | 29 |
| 2020-05-06 | O07930 | 145 | 145 | 54 |
| 2020-05-10 | O07991 | 742 | 677 | 318 |
| 2020-05-28 | O08261 | 142 | 137 | 16 |

공개 GFZ 자료에는 `degrade_include_flag`와 `elev_highestreturn_outlier_flag`가 없으며, 일부 NASA 궤도의 표본도 조회되지 않았다. 위 수치는 원본 다운로드 우선순위를 정하는 데 썼다. O07991 NASA HDF5는 받아서 `week08/source/`에 저장했다. [O07991 NASA 원본](https://data.ornldaac.earthdata.nasa.gov/ornl-cumulus-prod-protected/gedi/GEDI_L4A_AGB_Density_V3/data/GEDI04_A_2020131202827_O07991_03_T00715_02_004_01_V003.h5)

```powershell
& .\week04\.venv\Scripts\python.exe -B .\week08\inventory_spring.py --extended
& .\week04\.venv\Scripts\python.exe -B .\week08\evaluate_spring.py --extended --orbits O07991
```

## 2025년 봄철 검증

**봄철 새 궤도 두 건에서 week08의 GEDI 품질 기준을 모두 통과한 표본은 0점이었다.** 연구 구역(128.75–129.00°E, 35.00–35.25°N)에서 4월 16일 2,160점, 4월 29일 297점을 확인했다. 봄철 Sentinel-2 픽셀 매칭과 모델 평가는 입력 표본이 없어 수행할 수 없었다. [품질 단계별 집계 및 NASA SHA-512 검증](outputs/spring_gedi_quality_audit.json)

| GEDI 날짜 (UTC) | 궤도 | 구역 내 원본 관측점 | GEDI 품질 통과 | 그중 EGIS 산림 | 그중 고 AGBD | 인근 Sentinel-2 제품 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 4월 16일 | O35923 | 2,160 | 0 | 0 | 0 | 12 |
| 4월 25일 (학습) | O36071 | 4,265 | 1,055 | 187 | 52 | 12 |
| 4월 29일 | O36132 | 297 | 0 | 0 | 0 | 12 |

4월 16일은 L4A 품질 단계까지 3점이 남았고, 세 점 모두 `degrade_include_flag=1` 조건에서 제외됐다. 4월 29일은 AGBD·표준오차·알고리즘 실행 조건까지 11점이 남았고, 모두 `l2a_quality_flag_rel3=1` 단계에서 제외됐다. 4월 25일의 산림 187점은 Sentinel-2 선별 전 집계이며, 기존 학습 자료는 184점이다. 고 AGBD 기준은 100 Mg/ha 이상이다. [GEDI 목록](outputs/spring_gedi_inventory.csv) · [4월 16일 결과](outputs/spring_O35923_summary.json) · [4월 29일 결과](outputs/spring_O36132_summary.json)

NASA CMR에서 같은 구역과 2025년 4~5월에 겹치는 V003 파일은 위 세 건이었고, 이 검색에서 5월 파일은 조회되지 않았다. Sentinel-2 제품 수는 Earth Search에서 각 GEDI 파일 관측 구간의 앞뒤 7일과 연구 구역에 겹치는 제품을 센 값이다. [NASA 컬렉션](https://doi.org/10.3334/ORNLDAAC/2508) · [CMR 응답](outputs/gedi_april_may_cmr.json) · [Sentinel-2 제품 목록](outputs/spring_sentinel_candidates.csv)

## 확보된 2월 표본의 고도·경사 비교

4월 학습 184점과 기존 2월 평가 35점, RF 설정을 고정했다. Copernicus GLO-30 DSM 타일 두 장을 30m UTM 격자로 변환했다. 중심 셀의 고도와 동서·남북 인접 셀의 중앙차분으로 계산한 경사를 10개 반사율 변수에 추가했다. 이 자료는 지표면 고도에 식생과 구조물의 높이도 포함할 수 있다. 기존 10밴드 RF와 학습 평균 예측을 점별로 재현했다. [Copernicus GLO-30 자료](https://registry.opendata.aws/copernicus-dem/) · [지형 표본](outputs/topography_samples.csv)

| 모델 | RMSE (Mg/ha) | MAE (Mg/ha) | 고 AGBD 평균 과소 예측 (Mg/ha) |
| --- | ---: | ---: | ---: |
| 4월 학습 평균 | 42.77 | 33.83 | 57.38 |
| RF 10밴드 | 46.35 | 38.15 | 62.34 |
| RF 10밴드 + 고도·경사 | 46.01 | 35.57 | 60.34 |

고도·경사를 추가하자 기존 RF의 RMSE는 0.34 Mg/ha 낮아졌고, 학습 평균보다 3.23 Mg/ha 높았다. 고 AGBD 10점은 세 모델 모두 과소 예측했다. 이번 비교는 이미 확인한 2월 35점을 사용한 탐색 결과다. [점별 예측](outputs/topography_predictions.csv) · [지표와 원본 체크섬](outputs/topography_comparison.json)

저장소 루트에서 기존 환경으로 재현할 수 있다. NASA HDF5는 `week08/source/`에 직접 두거나 그 아래 Earthdata 다운로드 폴더에 둔다. 원본 DSM은 공개 AWS 버킷 `copernicus-dem-30m`의 `Copernicus_DSM_COG_10_N35_00_E128_00_DEM`과 `Copernicus_DSM_COG_10_N35_00_E129_00_DEM` 폴더의 동명 `.tif`를 `week08/source/`에 저장한다.

```powershell
& .\week04\.venv\Scripts\python.exe -B .\week08\inventory_spring.py
& .\week04\.venv\Scripts\python.exe -B .\week08\compare_topography.py
& .\week04\.venv\Scripts\python.exe -B .\week08\evaluate_spring.py
```

`evaluate_spring.py`는 원본을 확보한 궤도별로 GEDI 시각에서 168시간 이내인 제품을 비교한다. SCL 식생 통과점이 가장 많은 제품 한 개를 선택하고, 동률은 시간 차로 정한다. 원본 품질 마스크를 적용한 뒤 같은 표본에서 4월 학습 평균, 10밴드 RF, 10밴드에 고도·경사를 추가한 RF를 비교한다.
