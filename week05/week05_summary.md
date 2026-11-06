# Week05 공유 요약 — GLURI × 부산대학교 R&D 인턴십

## 1. 지난주 결과 및 현재 상황

Week04는 GEDI L4A AGBD 90점과 Sentinel-2 세 지수로 RF를 학습했다. 기존 RMSE는 Random **52.34**, Spatial **70.61 Mg/ha**였고, Spatial Dummy **57.16 Mg/ha**보다 나빴다([Step1 성능 CSV](step1/outputs/model_comparison.csv)).

## 2. 이번 주 수행한 작업

예측 잔차·데이터 품질을 분석하고, Sentinel-2 TIFF·마스크·좌표 결합과 `/65535` 변환의 출처를 추적했다([Step1](step1/outputs/report.md), [Step2](step2/report.md), [Step3](step3/report.md)). 입력 변수, RF 복잡도, GEDI 품질 조건을 탐색적으로 비교했으며([개선 실험](improvement/improvement_report.md)), GEDI L2A RH95 수고 연구는 문헌·자료 수준에서 검토했다([타당성 검토](gedi_l2a_feasibility.md)). L2A 모델 학습은 수행하지 않았다.

## 3. 주요 성과

RF 복잡도 축소 후 RMSE는 Random **49.65**, Spatial **62.35 Mg/ha**로 기존 RF보다 각각 **2.69**, **8.26** 낮아졌다. Random에서는 Dummy **50.94**보다 낮지만, Spatial에서는 Dummy **57.16**보다 높고 R²도 **−1.468**이다([개선 결과 CSV](improvement/experiment_results.csv)). 좌표·마스크·SCL 처리는 확인 범위에서 일관됐으며, 확인된 반사율 수치 오류는 없었다([Step2](step2/report.md)).

## 4. 어려운 점 및 지원이 필요한 사항

2025 TIFF 생성 요청·evalscript·제품 ID와 독립 FLOAT32 반사율이 없어 `/65535`의 물리적 타당성을 **판단할 수 없다**([Step3](step3/report.md)). 표본은 한 GEDI 궤도 90점으로 적고, 기존 평가점을 반복 확인했으므로 개선 수치는 탐색적이다([개선 보고서](improvement/improvement_report.md)). 원본 출력 설정, GEDI L2A RH95, 추가 궤도·독립 수고 검증자료의 확보 경로에 대한 지원이 필요하다.

## 5. 다음 연구 방향

먼저 같은 장면·픽셀의 검증 반사율로 TIFF 변환을 확인한다. 이후 추가 시기·지역에서 버퍼를 둔 공간 검증과 Dummy 비교를 진행하고, L2A RH95 표본 수·품질을 확인한 뒤 별도 수고 연구의 실행 여부를 결정한다([최종 보고서](week05_final_report.md), [L2A 검토](gedi_l2a_feasibility.md)).
