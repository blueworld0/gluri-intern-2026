# GLURI × 부산대학교 R&D 인턴십 Week05 최종 보고서

> 이 보고서는 [Step1](step1/outputs/report.md), [Step2](step2/report.md), [Step3](step3/report.md), [모델 개선 실험](improvement/improvement_report.md)의 실제 산출물과 CSV를 종합한다. 이번 최종 보고서 작성 중 추가 학습·외부 API 호출·기존 자료 변경은 하지 않았다.

## 1. 연구 개요 및 목적

GEDI L4A 추정 지상부 바이오매스 밀도(AGBD, Mg/ha)와 Sentinel-2 지수를 결합한 Week04 RF 기준선을 진단하고, 입력 변수·모델 복잡도·GEDI 품질 조건의 탐색적 개선 효과를 비교했다. 동시에 Sentinel-2 반사율 계보를 점검하고, 후속 GEDI L2A 수관 높이 연구의 실행 조건을 문헌으로 검토했다. **AGBD는 GEDI 제품의 모델 추정치**이며 현장 실측 정답이 아니다([Step1](step1/outputs/report.md)).

## 2. Week04 결과와 문제점

SCL 통과 **90점**에 `ndvi,evi,ndmi`를 입력으로 사용했다. Random 80/20 분할은 학습 72/평가 18점, Spatial은 동쪽 1km 블록 평가 20점과 500m 버퍼 제외 15점으로 학습 55점이다([Step1 비교 CSV](step1/outputs/model_comparison.csv), [Week04 학습 코드](../week04/train_baseline.py)). 기존 RF의 Random RMSE **52.34**, Spatial RMSE **70.61** Mg/ha이고, R²는 각각 **−0.139**, **−2.165**다. Spatial에서 Dummy보다 오차가 크며, 단일 궤도·소표본·세 지수라는 제약이 있다.

## 3. Week05 모델 성능·오차 분석

Step1은 Week04 RF 및 Dummy 지표를 재현하고 예측 잔차를 계산했다. 기존 RF의 Random 평가점 중 AGBD ≤50 Mg/ha는 3/3점 과대 예측, ≥150 Mg/ha는 2/2점 과소 예측했다. Spatial ≥150 Mg/ha는 7/7점 과소 예측했고 평균 잔차는 **+98.98 Mg/ha**였다([잔차 집계 CSV](step1/outputs/residual_summary.csv)). 구간별 점수가 작아 일반적 편향으로 확정할 수 없다. 최종 AGBD는 평균 **103.94**, 중앙값 **103.08**, 범위 **11.00–301.97 Mg/ha**다([Step1 보고서](step1/outputs/report.md)).

## 4. Sentinel-2 데이터 품질 검증

Step2는 2025 TIFF의 B02/B03/B04/B08/B11/B12 첫 밴드가 `uint16`, 둘째 밴드가 0/1 마스크임을 확인했다. 각 밴드의 마스크 유효 픽셀은 **909,212개**, 제외 픽셀은 **152개**다([밴드 통계 CSV](step2/band_quality.csv)). 최종 90점의 NDVI·EVI·NDMI는 각각 **0.301–0.909**, **0.145–0.772**, **0.050–0.493**이며 단순 범위 이탈은 없다([지수 품질 CSV](step2/index_quality.csv)). 90점은 원본·정렬 SCL에서 모두 식생 클래스 4이고, 여섯 TIFF의 EPSG:32652 격자와 점–픽셀 행·열 매핑은 확인 범위에서 일치한다([Step2 보고서](step2/report.md)). 이는 자동 SCL 오류나 GEDI footprint 전체의 공간 대표성이 없음을 증명하지는 않는다.

## 5. 반사율 변환 검증 결과

Week03의 `/65535`는 저장된 지수를 **계산상 재현**하지만, TIFF의 `scale=1`, `offset=0`, 비어 있는 NoData·단위·제품 ID만으로 물리 반사율 변환식은 확정할 수 없다. 2025 원본 요청의 `units`, `sampleType`, `harmonizeValues`, evalscript 및 독립 FLOAT32 밴드 반사율이 현재 저장소에서 확인되지 않았다([Step3 보고서](step3/report.md), [추적 문서](step3/source_trace.md)). TIFF 파일명의 2025-04-28은 요청 날짜 범위이며 정확한 촬영 시각으로 확정하지 않았다. **`/65535` 판정은 근거 부족(C)**이고, 확인된 반사율 수치 오류도 없다. `/10000`을 검증된 대체 계수로 적용하지 않는다. EVI는 공통 배율에도 달라질 수 있고, 밴드별 offset·클리핑은 NDVI·NDMI·NBR에도 영향을 줄 수 있다([변환 평가](step3/conversion_assessment.md)).

## 6. 실제 모델 개선 실험

추가 입력(A)은 원시 B03·B12와 B08/B12의 NBR 후보를 학습 내부 CV에서 선택했다. 복잡도 축소(B)는 `max_depth`, `min_samples_leaf`, `max_features`의 제한된 후보를 각 외부 학습 집합의 500m easting 그룹 3-fold CV로 선택했다. Random 선택값은 **깊이 3·잎 4·max_features 1.0**, Spatial은 **깊이 6·잎 3·max_features 0.67**이다([모델 비교](improvement/model_comparison.md), [내부 CV](improvement/inner_cv_results.csv)). 품질 조건(C)은 학습점에만 `sensitivity≥0.95`를 적용했다. 전체 90점 중 45점이 조건을 충족했고, Spatial 학습점은 55점에서 10점으로 줄어 사전 최소 30점 기준에 못 미쳐 **Spatial C는 실행하지 않았다**([품질 분포](improvement/quality_distribution.csv), [개선 보고서](improvement/improvement_report.md)).

## 7. 개선 전후 성능 비교

아래 A·B와 기존 RF·Dummy는 분할별 **동일 평가점**을 사용한다. RMSE·MAE 단위는 Mg/ha이며, 원값은 [개선 실험 CSV](improvement/experiment_results.csv)에 있다.

| 분할 | 모델·실험 | 학습/평가 n | R² | RMSE | MAE |
|---|---|---:|---:|---:|---:|
| Random | 기존 RF | 72/18 | −0.139 | 52.34 | 42.87 |
| Random | Dummy | 72/18 | −0.080 | 50.94 | 43.34 |
| Random | 추가 입력 A | 72/18 | −0.159 | 52.78 | 43.03 |
| Random | 복잡도 축소 B | 72/18 | −0.025 | **49.65** | **41.17** |
| Spatial | 기존 RF | 55/20 | −2.165 | 70.61 | 57.71 |
| Spatial | Dummy | 55/20 | −1.074 | **57.16** | **48.07** |
| Spatial | 추가 입력 A | 55/20 | −2.069 | 69.52 | 56.27 |
| Spatial | 복잡도 축소 B | 55/20 | **−1.468** | 62.35 | 50.47 |

복잡도 축소로 기존 RF 대비 RMSE가 Random **2.69**, Spatial **8.26 Mg/ha 감소**했다. Random에서는 Dummy보다 RMSE **1.29 Mg/ha 낮지만**, Spatial에서는 Dummy보다 **5.19 Mg/ha 높다**. 두 분할의 RF R²는 여전히 음수다. A는 Random에서 RMSE가 **0.44 악화**, Spatial에서는 **1.08 개선**됐지만 Dummy를 넘지 못했다([비교표](improvement/model_comparison.md)).

C의 Random은 학습점 **72→39**, 평가점 18개 유지 상태에서 RF RMSE **67.62**, 변경된 학습점의 Dummy RMSE **64.15**였다. 학습 표본과 AGBD 분포가 바뀌었으므로 A·B처럼 동일 학습 조건의 개선 효과로 해석할 수 없다. Spatial C의 성능 수치는 없다([실험 CSV](improvement/experiment_results.csv), [품질 분포](improvement/quality_distribution.csv)).

## 8. GEDI L2A 수관 높이 연구 검토

[별도 실행 가능성 검토](gedi_l2a_feasibility.md)는 프로젝트의 [논문 PDF](../week04/research/paper1.pdf)와 [재현 계획](../week04/research/paper1_reproduction_plan.md)을 근거로 작성했다. 논문은 **GEDI L2A RH95**를 수고 라벨로 쓰고 `quality_flag=1`, `degrade_flag=0`, `sensitivity≥0.95`를 적용하며, Sentinel-2·Sentinel-1·경사·토지피복을 10m 격자에서 결합한다(PDF pp. 5–6). 같은 픽셀의 RH95 복수 샷은 평균한다. 수고 평가는 홀드아웃 RH95와 독립 CHM·LiDAR·현장 수고를 구분해야 한다. **현재 프로젝트에는 L2A RH95 원자료와 수고 모델 학습·검증 결과가 없다.** 논문의 RMSE(m)와 이번 AGBD RMSE(Mg/ha)는 직접 비교하지 않는다.

## 9. 연구의 한계

- 90점은 한 GEDI 파일·한 관측일의 제한된 궤도이며 Random 평가 18점, Spatial 평가 20점이다([Step1](step1/outputs/report.md)). 한 동쪽 Spatial 홀드아웃은 다른 지역 일반화의 충분한 증거가 아니다.
- Week04·Step1에서 평가점을 이미 여러 번 확인했다. 개선 후보 선택은 학습 내부에서 했지만, 이번 점수는 **독립적인 최종 일반화 성능이나 통계적 유의성 증거가 아니다**([개선 보고서](improvement/improvement_report.md)). 내부 CV의 공간 그룹에는 500m 버퍼가 없다.
- Sentinel-2 방사 변환과 정확한 TIFF 촬영·제품 연결은 미검증이며, GEDI AGBD도 독립 현장값이 아니다([Step3](step3/report.md)).
- GEDI L2A 연구는 문헌 검토 단계다. L2A 샷 수, RH95 품질, 외부 수고 참조자료를 아직 평가하지 않았다([실행 가능성 검토](gedi_l2a_feasibility.md)).

## 10. 후속 연구 계획

1. 2025 TIFF 생성 요청·evalscript·제품 ID와 같은 장면의 검증된 FLOAT32 반사율을 확보해 동일 픽셀·밴드에서 `/65535`와 비교한다. 제품·촬영 시각과 SCL 장면 연결도 확인한다([Step3 교차 검증 절차](step3/reflectance_crosscheck.md)).
2. 새 장면·GEDI 궤도와 독립 지역·시기를 확보하고, 학습 내부에는 버퍼가 있는 공간 CV, 최종에는 미사용 홀드아웃을 두어 RF와 Dummy를 다시 비교한다([개선 보고서](improvement/improvement_report.md)).
3. L2A RH95의 접근 가능성·품질·표본 수부터 확인하고, 충분할 때 별도 수고 학습표와 공간·시간 정합, 독립 검증 설계를 만든다([L2A 검토](gedi_l2a_feasibility.md)). 이 계획은 이번 주에 수행한 실험이 아니다.

## 11. 대표자 피드백 요청 사항

1. **자료 접근:** 2025 Sentinel-2 Browser 다운로드 설정·evalscript·장면/제품 ID 또는 동일 장면 FLOAT32 반사율과 L2A 원자료 접근 경로를 확보할 수 있는가?
2. **연구 우선순위:** 방사 값 검증과 AGBD 데이터 확장을 먼저 진행할지, L2A RH95 수고 타당성 조사에 별도 자원을 배정할지 결정이 필요한가?
3. **검증 자료·범위:** 금정산 경계, 추가 궤도·기간, 독립 수고 기준(CHM·LiDAR·현장 조사)과 NFI/임상도 자료의 사용 가능 범위를 확인해 줄 수 있는가?
4. **성과 해석:** 이번 복잡도 축소 결과를 탐색적 내부 성과로만 보고하고, 공간 Dummy 미달 및 반사율 미확정을 최종 성능 주장에 명시하는 방향에 의견이 있는가?
