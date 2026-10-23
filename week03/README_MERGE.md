# GLURI Week03 — Sentinel-2 + GEDI04_A 공간 결합

`merge_sentinel_gedi.py`는 기존 [GEDI 유효점 CSV](outputs/gedi/gedi_valid_points.csv)의 각 위경도에서 NDVI·EVI·NDMI GeoTIFF 픽셀값을 추출합니다. AGBD 모델은 학습하지 않습니다. 기존 Sentinel-2·GEDI 코드, 원본, 출력은 변경하지 않습니다.

## 실행

```powershell
.\.venv\Scripts\python.exe merge_sentinel_gedi.py
```

기본 출력은 `outputs/merged/`입니다. 결과 파일이 이미 있으면 덮어쓰지 않고 중단합니다. 새 실행 결과가 필요하면 `--output-dir outputs/merged_run2`처럼 다른 폴더를 지정하세요.

## 방법

1. GEDI CSV의 실제 컬럼을 확인하고 위도·경도, AGBD, `agbd_se`, 관측시각, 파일명, BEAM, `shot_number`, 품질 플래그를 명시적으로 매핑합니다. 필수 컬럼이 없거나 매핑이 모호하면 중단합니다.
2. NDVI·EVI·NDMI GeoTIFF의 CRS, 크기, Transform, 범위가 서로 일치하는지 확인합니다. 셋 모두 EPSG:32652, 578×834, 약 9.9975×9.9941m, FLOAT32, NoData=-9999이며 실제 유효 픽셀 비율은 100%입니다. 촬영일은 GeoTIFF 태그의 2026-09-12를 사용합니다.
3. `rasterio.warp.transform`에 **경도, 위도 순서**로 입력해 UTM 52N 좌표를 얻고, Transform으로 픽셀 행·열을 계산합니다. 이후 별도 `project_and_sample` 함수가 좌표별 값을 읽습니다. NoData, NaN, 범위 밖, NDVI·NDMI의 [-1,1] 이탈을 검사하고 실패 사유를 남깁니다. EVI는 [-1,1]로 자르지 않습니다.
4. `point_id`는 GEDI 파일명·BEAM·`shot_number`를 결합합니다. 동일 픽셀 중복과 가까운 관측점 수를 별도로 확인합니다.
5. 날짜 차이는 Sentinel-2 촬영일에서 GEDI UTC 관측일을 뺀 일수입니다. 날짜 차이를 이유로 자동 삭제하지 않습니다.

## 실제 결합 결과

- GEDI 입력 97개 중 범위 안 97개, NDVI·EVI·NDMI 추출 각각 97개, 세 지수 모두 유효한 표본 97개입니다. 실패 0개, 동일 픽셀 중복 0개입니다.
- 최종 `sentinel_gedi_clean.csv`는 **97행 × 41열**입니다. 주요 모델 입력인 NDVI·EVI·NDMI·AGBD·`agbd_se`에는 결측이 없습니다. 선택적 예측구간 하한 `agbd_pi_lower_Mg_ha`에는 결측 10개가 있습니다.
- GEDI 관측일은 전부 2025-04-25이며 Sentinel-2와 날짜 차이는 전부 **505일**입니다.
- 96개 점의 최근접 다른 점까지 거리가 60m 미만입니다. 최소 최근접 거리는 약 54.43m입니다. 무작위 train/test 분할은 공간적으로 가까운 점을 양쪽에 넣어 데이터 누수를 일으킬 수 있습니다.

| 변수 | 평균 | 중앙값 | 표준편차 | AGBD Pearson | AGBD Spearman | 상관 표본 수 |
|---|---:|---:|---:|---:|---:|---:|
| NDVI | 0.8257 | 0.8773 | 0.1586 | 0.1374 | 0.2725 | 97 |
| EVI | 0.4990 | 0.5215 | 0.1059 | 0.0693 | 0.0611 | 97 |
| NDMI | 0.3867 | 0.4220 | 0.1200 | 0.2280 | 0.3105 | 97 |
| AGBD (Mg/ha) | 103.1914 | 102.5772 | 59.2637 | — | — | 97 |

상관관계는 탐색적 기술통계이며 인과관계나 모델 성능을 뜻하지 않습니다. GEDI AGBD는 제품 모델의 예측값이지 독립적인 현장 실측값이 아닙니다. 이 표본은 Sentinel-2의 사각형 범위 안의 일부 관측점이며 금정산 전체를 대표하지 않습니다.

## 출력

- `outputs/merged/sentinel_gedi_merged.csv`: 97개 관측점의 결합 결과와 픽셀 행·열, 추출 상태, 실패 사유.
- `outputs/merged/sentinel_gedi_clean.csv`: 주요 입력과 AGBD가 모두 유효한 표본.
- `outputs/merged/merge_statistics.csv`: 기초 통계, Pearson·Spearman 상관계수와 표본 수.
- `outputs/merged/merge_quality_report.json`: 사용 파일, 래스터 메타데이터, 결합 성공·실패, 중복·결측, 날짜 차이, 이상값 진단.
- `outputs/merged/figures/`: NDVI 지도 위 GEDI AGBD 점, 지수별 산점도, Pearson 히트맵, 날짜 차이 분포.

## 다음 단계 판단

공간 조인 자체는 재현 가능하고 97개 행을 가진 탐색용 데이터셋을 확보했습니다. 현재 데이터만으로 AGB 모델 재현이나 신뢰할 만한 검증을 진행하기에는 **505일의 시차, 한 관측일에 집중된 표본, 선형 공간 배치, 작은 유효 표본 수**가 큰 제약입니다. GEDI 관측 시기와 가까운 Sentinel-2 영상, 추가 GEDI 날짜·궤도, 산림 경계 및 구름·그림자 품질 자료를 확보하고 공간적으로 분리된 검증 계획을 세워야 합니다.
