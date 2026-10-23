# GLURI Week03 — GEDI04_A V003 분석

기존 Sentinel-2 분석과 별도로 `GEDI04_A/`의 원본 HDF5 5개를 읽어 Sentinel-2 영상의 **사각형 범위에 포함되는 GEDI 관측점**을 분석합니다. 이 범위는 금정산 산림 경계가 아닙니다. GEDI의 AGBD는 제품 모델의 예측값이며 현장 실측 정답값이 아닙니다.

## 실행

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe analyze_gedi.py
```

`outputs/gedi/`에 결과가 이미 있으면 덮어쓰지 않고 종료합니다. 새 결과를 만들 때는 별도 경로를 지정합니다.

```powershell
.\.venv\Scripts\python.exe analyze_gedi.py --output-dir outputs/gedi_run2
```

`gedi_all_points.csv`는 6,356,011개 전체 관측점을 담으므로 약 1.25GB입니다. 저장 공간과 실행 시간을 고려하세요.

## 좌표, 시간, 품질 기준

- 기준 영상은 Sentinel-2 B04 GeoTIFF입니다. EPSG:32652의 범위는 easting 501553.909–507332.449m, northing 3897622.903–3905957.956m입니다.
- 위치는 각 BEAM의 `lat_lowestmode`, `lon_lowestmode`를 사용합니다. 후보점을 WGS84에서 EPSG:32652로 변환한 후 정확한 사각형 범위로 판별합니다.
- 관측 시각은 실제 `delta_time`의 속성인 `seconds since 2018-01-01`을 사용하여 UTC로 변환합니다.
- [NASA ORNL DAAC GEDI L4A V3 가이드](https://daac.ornl.gov/GEDI/guides/GEDI_L4A_AGB_Density_V3.html)에 맞춰 `l2_algrunflag=1`, `l2a_quality_flag_rel3=1`, `l4a_quality_flag_rel3=1`, `degrade_include_flag=1`, `elev_highestreturn_outlier_flag=0`을 적용합니다. `degrade_flag` 원시 코드는 CSV에 보존합니다. V3에서는 데이터 융합용 위치 품질에 `degrade_include_flag`와 고도 이상치 플래그를 함께 사용하도록 권고합니다.
- `agbd`와 `agbd_se`가 유한하고 0 이상이어야 합니다. 실제 데이터에서 관찰된 `-9999` 결측 표시는 비웁니다. HDF5의 기본 할당 fill 값 0과 관찰된 `-9999`는 구분합니다. `predictor_limit_flag`와 `response_limit_flag`는 기록하되 NASA가 보편적 제외 기준을 제시하지 않아 자동 제외하지 않습니다.

## 결과

| 범위/단계 | 관측점 수 |
|---|---:|
| 5개 파일 전체 | 6,356,011 |
| Sentinel-2 사각형 안 | 2,090 |
| 유효 AGBD | 864 |
| L2A 품질 통과 | 814 |
| L4A 품질 포함 최종 통과 | 97 |

품질 유효점 97개는 모두 **2025-04-25** 관측으로, 2026-09-12 Sentinel-2 영상보다 약 16개월 전입니다. AGBD는 9.37–301.97 Mg/ha, 평균 103.19, 중앙값 102.58, 표준편차 59.26 Mg/ha입니다. `agbd_se` 평균은 8.29 Mg/ha입니다. 97개 중 예측구간 하한은 87개에서만 유효합니다. 관측점 평균을 금정산 전체 평균 또는 탄소 저장량으로 해석하지 않습니다.

## 저장 파일

- `outputs/gedi/gedi_all_points.csv`: 파일 전체 관측점, 위치·시각·AGBD·품질·불확실성 정보.
- `outputs/gedi/gedi_geumjeongsan_points.csv`: Sentinel-2 사각형 안의 모든 관측점. 실제 산림 경계로 필터링한 결과가 아닙니다.
- `outputs/gedi/gedi_valid_points.csv`: V3 품질 필터를 통과한 관측점.
- `outputs/gedi/gedi_statistics.csv`: AGBD 및 불확실성 기초 통계.
- `outputs/gedi/gedi_counts_by_file_date.csv`: 파일·관측일별 전체, 영역 안, 유효 관측점 수.
- `outputs/gedi/gedi_quality_counts.csv`: 단계별 필터 통과 수.
- `outputs/gedi/figures/`: 전체점, 유효점, AGBD 색상 지도와 히스토그램. 지도 배경은 Sentinel-2 RGB입니다.

## 다음 단계

유효점 97개는 서로 다른 Sentinel-2 픽셀 97개에 대응하고 NDVI NoData는 없습니다. 따라서 위치 기반 결합은 가능합니다. 다만 GEDI와 Sentinel-2의 촬영 시점 차이, 적은 표본 수, 관측점의 선형·편중된 분포를 고려해 학습 데이터로 사용할지 판단해야 합니다. 산림 경계와 구름·그림자·피복 품질을 추가한 뒤 공간적 검증을 설계하는 것이 필요합니다.
