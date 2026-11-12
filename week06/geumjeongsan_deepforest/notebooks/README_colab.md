# 금정산 Colab 노트북

파일: `GEUMJEONGSAN_DeepForest_2021.ipynb` — 32개 셀, 기존 NEON 노트북과 분리된 새 파일입니다.

## 실행

1. Colab에 이 노트북 하나를 올리고 T4 GPU를 선택합니다. 보조 Python 파일을 함께 올릴 필요는 없습니다.
2. 노트북 00번 안내대로 Drive의 `MyDrive/geumjeongsan/gedi_geumjeongsan/data/`에 입력자료를 업로드합니다.
   정사영상 9장·XML 9개, 야장 CSV 3개, 계수 CSV 3개·README, 임상도와 등산로의 shapefile 부속 파일을 함께 유지하세요.
3. 패키지 설치 후 재시작 안내가 나오면 재시작하고 Drive 마운트부터 실행합니다.
4. 06번에서 도엽 좌표와 기준 bbox를 검사합니다. 오차가 기준을 넘거나 50개 crop을 덮지 못하면 추론 전에 중단합니다.
5. 입력 필드·원문 계수 출력을 확인한 뒤 crop·추론·비교 셀을 실행합니다.
6. 12번의 임상도·등산로 겹침 그림 2장을 눈으로 확인하세요. 수치 bbox 통과만으로 사진 정합이 확정되지는 않습니다.
7. 마지막 요약과 보고서를 확인한 뒤 ZIP 셀을 실행하면 `/content/geumjeongsan_results.zip`을 다운로드합니다.

Drive 자료를 `MyDrive/geumjeongsan/data/`에 바로 놓은 경우도 인식합니다. 원본 TIFF는 수정하지 않습니다.

## 재사용과 캐시

- NEON 3차 `neon_pipeline.py`를 노트북에 내장하고 `predict_crowns`의 탐지·수고·상자 CD를 재사용합니다.
- DeepForest 2.1.0, score 0.2, patch 600, overlap 0.25, NMS 0.15, 내부 score floor 0.05, H≥3m입니다.
- 확인한 NEON 해상도 실험 HF 커밋은 `cc21436bc5d572dde8ff5f93c1e71a32f563cace`입니다.
- HF 캐시는 `/content/hf_cache`입니다. 정사영상 crop, 10cm 확대본, CHM crop과 탐지는 Drive에 저장합니다.
- Meta CHM은 zoom9 quadkey를 AOI로 계산하고, 기존 native 타일이 있으면 먼저 사용합니다.
  없으면 해당 COG의 필요한 window만 읽어 crop을 캐시합니다. 전 세계 타일이나 거대한 quad 전체를 받지 않습니다.
- H는 재투영 전 native CHM 픽셀의 최대입니다. EPSG:5187 1m CHM도 별도로 저장하며 이 파일은 표시용입니다.
- HTTP 다운로드 재시도, gzip 등 인코딩 응답의 크기 검사 예외, 원시 탐지 완료 표식과 입력 변경에 따른 캐시 분리를 포함합니다.

## 비교 정의

- **지상부만** 계산합니다. 국내 식의 줄기·가지·잎 합계 × CF 0.5이며 뿌리는 제외합니다.
- 공식 비교는 중심 기준 축정렬 20m×20m, 실측 DBH≥10cm입니다. DBH≥6cm 결과도 저장합니다.
- 영상 탐지는 NEON과 같은 H≥3m 기준을 유지합니다. 추정 DBH≥10cm만 남긴 보조 추정치도 따로 기록합니다.
- 30m 원은 위치 오차 민감도이며, 원 안의 실측을 새로 관측한 값은 없습니다. 20m 조사구의 참조 밀도를 대용한 것으로 표시합니다.
- T1 최소 수고 이상을 보이는 나무로 대리 분류합니다. 개별 줄기 위치가 없어 개체 매칭과 F1을 계산하지 않습니다.
- 임상도 미피복·비산림·죽림 등에서 임상을 정하지 못한 탐지는 혼효 대표식 대용과 `forest_fallback`으로 표시합니다.
  이 정책은 보고서에 명시하며 지도 밖을 탄소 0으로 간주하지 않습니다.
- zone 5개 중 1개를 통째로 제외하고 나머지 4개에서 실측/추정 비율을 적합합니다. NEON 기준선은 보정에 사용하지 않습니다.
- NEON 25cm 고정 설정의 실제 비율 0.52799277, F1 0.44061303, RMSE 94.1765tC/ha를 참고행으로 기록합니다.

## 생성 결과

- `georef_check.csv`: 도엽 9개 범위·크기 차이·35913003 기준 오차
- `plot_comparison.csv`: 기본 비교 50행, 참조 ≥6cm와 민감도 추정치 보조 열
- `plot_comparison_sensitivity.csv`: 50개 × 5시나리오 × 2직경 기준 = 500행
- `summary.csv`: 전체·5zone·NEON 참고 7행, `summary_sensitivity.csv`: 전체 민감도 집계
- `detected_trees.csv`, `field_reference_trees.csv`: 사용 수종·부위 식·CF·적용범위·보임 여부
- `calibration_folds.csv`: 각 시나리오의 train/test zone·조사구·보정계수
- `provenance.json`, `raster_manifest.csv`, `meta_chm_sources.json`, `vector_fields.csv`, 각종 검산·대용 기여 표
- `figures/`: overview, 정합 확인 2장, 탄소 산점도, 가시/가림 막대, 5zone와 낮은/높은 탄소를 포함한 조사구 6개 overlay
- `report.md`: 결과 한 줄 요약·방법·공식·불확실성·출처

## 로컬 검증 범위

Windows에서는 pandas·rasterio·geopandas 등을 import하지 않았습니다.
표준 라이브러리로 TIFF 헤더·XML·독립 TM 급수·국내 공식과 노트북 문법/순서를 확인했습니다.
실제 PyProj bbox assert, 영상 crop, Meta CHM 읽기와 GPU 추론은 Colab 셀에 있으며 아직 실행하지 않았습니다.

사전 검산: 35913003 기준의 각 변 오차 최대 약1.66m, 영상 9개의 폭·높이 차이 5m 이내, 50개 중심±60m 모서리 범위 포함.
실측 기반 평균 지상부 탄소는 ≥10cm **61.7732tC/ha**(1,422그루), ≥6cm **62.9178tC/ha**(1,694그루)입니다.
이 값은 영상 추론 성능 결과가 아닙니다. 세부 검산은 `GEUMJEONGSAN_VALIDATION.json`에 있습니다.
