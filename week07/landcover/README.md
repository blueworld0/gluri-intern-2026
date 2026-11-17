# EGIS 토지피복지도 준비

사용자가 받은 2025년 중분류 6개 도엽을 `source/<도엽번호>/`에 보존했다. 각 도엽의 `.shp`, `.shx`, `.dbf`, `.prj`는 모두 필요하므로 24개 파일을 남겼다. ZIP 외에 버린 내부 파일은 없다.

ZIP CRC, 추출 전후 바이트 일치, SHP/SHX 레코드·오프셋·좌표 유한성·링 폐합, DBF 레코드 수와 필드, PRJ 좌표계를 검사했다. 좌표계와 속성 목록·코드 샘플은 `provenance.json`에 기록했다. 도형 위상 유효성을 검사한 것은 아니다.

441점 전체가 6개 파일의 결합 bbox 안에 포함됨을 확인했다. 이는 실제 polygon 포함이나 산림 판정을 의미하지 않는다. 피복 분류와 경계 검토는 다음 분석에서 수행한다. 촬영·갱신 시기는 ZIP 이름만으로 확정하지 않았다. DBF 문자 인코딩은 확인 없이 새 CPG를 만들지 않고 원본을 보존했다.

사용자 요청에 따라 검증 완료 후 이번 ZIP 6개를 삭제했다. 원본 ZIP의 파일명·크기·SHA-256과 보존 파일의 SHA-256·CRC를 `provenance.json`에 남겼다. 기존 주차 데이터·분석 결과는 변경하지 않았다.

## 후속 피복 결합

`classify_points.py`로 실제 polygon 포함과 주변 피복을 분석했다. 결과는 [outputs/report.md](outputs/report.md)와 `outputs/landcover_points.csv`에 있다. 중심점 기준 산림 184점·비산림 257점·판정 불가 0점이며 기존 예측을 수정하거나 재학습하지 않았다. 산림은 EGIS 토지피복 코드에 따른 판정이고 현장 검증을 의미하지 않는다.

`verify_classification.py`는 원본 SHP 링에 대한 별도의 even/odd 교차 알고리즘으로 전체 441점 결합을 독립 대조했다. 도형의 홀·다중 외곽 사례와 실제 64,663개 레코드를 확인해 Shapely 공간 결합과 모두 일치했다. 결과는 `outputs/independent_verification.json`에 있다.

기존 Python 가상환경에 Shapely 2.2.0을 추가했다. CP949로 원본 속성 문자열을 엄격하게 읽고 전체 분류 코드·명칭 쌍의 일치를 확인했다. 준비 단계의 `provenance.json`은 당시 검사 기록으로 보존하며, 최신 도형 유효성·좌표계·주변 면적 검증은 `outputs/validation.json`을 참조한다.

## 산림 표본 전용 재학습

`evaluate_forest.py`는 EGIS 산림 중심점 184점으로 기존 fold·버퍼 역할과 RF 설정을 유지해 재학습한다. 주변 혼합 22점은 유지한다. [forest_evaluation/report.md](forest_evaluation/report.md)에 같은 산림 평가점에서 기존 441점 학습 모델과 비교한 결과를 기록했다. 산림 전용 RF도 평균 Dummy보다 성능이 낮아 현재 입력·모델의 산림 내 예측력은 확보되지 않았다. fold5 평가점은 4점뿐이라는 제한이 있다. 원본과 기존 결과는 변경하지 않았다.

## 원형 영역 평균 반사율 비교

`compare_footprint.py`는 같은 산림 184점의 반경 12.5m와 10m 픽셀이 겹치는 면적으로 여섯 밴드를 평균한 뒤 지수를 계산한다. 184점 모두 원 전체의 공통 밴드·SCL·CLASSI 유효성을 통과했다. 원 내부를 산림 피복으로 자르지 않았다. 동일 fold·파라미터의 결과는 [forest_evaluation/footprint_report.md](forest_evaluation/footprint_report.md)에 정리했다. 밴드6 RMSE는 약 2.18% 낮아졌지만 지수3은 악화되고 모든 RF가 Dummy보다 나빴다. 기존 결과는 보존하고 `footprint_*` 파일로 구분했다.

`verify_footprint.py`는 저장 가중치로 전체 평균·지수를 재계산하고, 별도 원 방정식 수치 적분으로 21점의 픽셀 면적을 대조한다. 원형 평균은 실제 GEDI footprint나 에너지·위치 오차를 정확히 복원한 값이 아니다.

## 기존 red-edge·B8A 추가

`compare_rededge.py`는 원본 SAFE B05·B06·B07·B8A의 XML 반사율 보정·품질 마스크·보간 유효성을 적용하고 동일 원형 면적 평균으로 기존 밴드6에 추가한다. 같은 184점·분할·RF 설정에서 밴드10 RMSE는 54.55로 밴드6의 55.51보다 약 1.74% 낮지만, Dummy 51.53보다 나쁘다. 모든 fold에서 Dummy보다 RMSE가 높다. 결과는 [forest_evaluation/rededge_report.md](forest_evaluation/rededge_report.md)에 있다.

`verify_rededge.py`는 원본 20m 픽셀 4개의 수동 bilinear 계산과 원형 가중 평균으로 새 밴드 추출을 독립 대조한다. 전체 184점이 신규 네 밴드의 원 전체 유효성을 통과했다. 외부 다운로드·튜닝·기존 데이터 수정은 하지 않았다.
