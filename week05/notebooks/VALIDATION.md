# 사전 검증 기록

- 노트북 JSON v4 구조, 고유 셀 ID: 통과
- 전체 46개 셀, 코드 22개: Python 문법 검사 통과 (첫 pip 셀은 Colab magic으로 별도 확인)
- 모든 코드 셀 바로 위 한국어 설명: 확인
- 실행 결과와 실행 번호: 비어 있음 (실행 완료로 오인하지 않도록 유지)
- 재사용 함수 모듈: 문법 검사 통과
- m/cm·kg/t·㎡/ha 단위 및 높이 민감도 수치 검산: 통과
- 로컬에서 pandas/geopandas/rasterio/DeepForest를 import하거나 실행하지 않음
- 제공된 1차 outputs/로 집계를 검산함. 수정본 다운로드·GPU 추론은 Colab 재실행 필요

공식 API 문서에서 load_by_product(dpid, site, startdate, enddate, package, release, token),
by_tile_aop(dpid, site, year, easting, northing, buffer, savepath, include_provisional, skip_if_exists, overwrite),
DeepForest load_model/use_release 및 predict_tile path/raster_path 차이를 확인했습니다.

Jucker 계수는 itcSegment/R/dbh.R의 biome=11 (a=0.367, b=0.887, g=0)을 확인했습니다.
Chojnacky 기본값은 혼합활엽수라는 독립 그룹이 아닌 일반 활엽수 대체 그룹
(β0=-2.2118, β1=2.4133)이며 원 논문 표와 MasterTaxonList 열 이름을 확인했습니다.
실행 시 원문 계수표를 캐시하고 계수와 해시를 다시 확인합니다.

## 1차 수정 검증

- 공식 API 직접 조회: ORNL RGB·CHM 모두 2021-06 존재, PROVISIONAL. ORNL UTM 16N 확인
- 실제 API 응답 기반 연도 선택: 2021 유지, include_provisional=True; 자료가 없을 때만 ORNL 2022→2020
- 좌표 누락이어도 참조 AGB가 있으면 조사구 총량 유효. 실제 AGB 누락일 때만 NaN: 회귀 검사
- 공식 면적 우선과 누락 분해 항등식: 회귀 검사
- 노트북과 neon_pipeline.py의 공통 함수 AST 일치: 확인
- 탐지·Jucker·Chojnacky·매칭·회귀·민감도 코드 원본 보존: check_revision.py에서 검사
- report.md 상단과 마지막 요약은 동일 표 사용, 마지막 ZIP 셀 존재: 확인
- GRSM 15개 검산: outputs/revision_check/README.md 참고

## 나무 찾기 개선 추가 검증

- score 5개×patch 3개, 고정 23개 ID, 요청된 4개 나란히 그림 ID: 설정 검사 통과
- 바깥 GroupKFold 평가 정답과 안쪽 선택 자료 분리, OOF 결과와 배포용 전체 선택 조합 분리: 소스 확인
- B: Gaussian·가변 창·watershed·면적 기반 CD, C: 중복 제거·상자 밖 꼭대기 추가: 함수/공식 확인
- 네 방법 전체/사이트 비교표와 탄소 RMSE 최소 기본값이 이후 민감도 이전에 연결됨: 확인
- QUDO 경고/제외, gzip 크기 검사 예외, /content/hf_cache, ORNL provisional, 전체 참조 총량: 유지 검사 통과
- 합성 두 나무 CHM의 꼭대기·watershed 면적·CD·결합 중복 제거와 평가 정답 변경 시 선택 격리 검사를 Colab 셀 14B에 추가
- 합성 검사와 실제 23개 자료의 GPU/CHM 연산은 로컬에서 실행하지 않음. 새 성능 수치를 주장하지 않음
