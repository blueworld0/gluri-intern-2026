# Week04: AGBD 예측 베이스라인

Week03의 SCL 품질 검증을 통과한 식생 후보 90개로 탐색용 Random Forest AGBD 베이스라인을 평가한다.

## 데이터 연결

- 입력: [`../week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv`](../week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv)
- `load_data.py`는 스크립트 위치를 기준으로 위 상대경로를 구성한다. 실행 위치와 관계없이 같은 파일을 읽는다.
- 원본 CSV는 Week03에 그대로 두며 Week04로 복사하지 않는다.
- SCL 품질 보고서의 최종 식생 후보 수는 90개다. CSV에는 `agbd`, `ndvi`, `evi`, `ndmi`, `latitude`, `longitude`, `easting_m`, `northing_m` 컬럼이 있다.
- 실행 확인 결과 CSV는 90행이다. AGBD·세 식생지수·좌표에는 결측치가 없다. 다른 컬럼 중 `agbd_pi_lower_Mg_ha`는 9개, `scl_class`는 90개가 비어 있다.

## 실행

프로젝트 루트에서 PowerShell로 실행한다.

```powershell
cd week04
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe load_data.py
```

스크립트는 전체 행 수, 컬럼 목록, AGBD·식생지수 기초 통계, 각 컬럼의 결측치 수를 출력한다.

SCL 클래스 4는 식생 분류이며 산림 여부를 보증하지 않는다. Week03의 2025년 지수 계산에 사용된 반사도 계수는 잠정적이다.

## Sentinel-2 반사도 변환 검증

`validate_reflectance.py`는 Week03의 2025년·2026년 Analytical Raw 밴드 TIFF 각 6개를 읽기 전용으로 검사한다. 파일별 메타데이터와 유효 픽셀 통계를 기록하고, 2025년 SCL 통과 점 90개에서 두 **가정**을 비교한다.

```powershell
.\.venv\Scripts\python.exe validate_reflectance.py
```

결과는 [`outputs/reflectance_validation/`](outputs/reflectance_validation/)에 생성된다. `tiff_metadata.json`은 원본 태그·타입·NoData·scale·offset을, `band_statistics.csv`는 밴드별 유효 픽셀 통계를, `candidate_index_comparison.csv`와 `index_comparison_statistics.csv`는 90개 점의 지수 비교를, `index_scale_sensitivity.png`는 비교 그래프를 담는다.

### 확인한 근거와 한계

- 실제 12개 TIFF 모두 두 밴드 `uint16`이다. 두 번째 밴드는 0/1 dataMask다. TIFF에 NoData 값, 반사도 단위, 유효한 반사도 scale/offset 태그는 없다. GDAL에서 읽히는 scale=1·offset=0은 항등값이며 반사도 변환 근거가 아니다.
- 2025년 각 밴드는 총 909,364픽셀 중 dataMask=1이 909,212개, 2026년은 총 482,052개가 유효하다. 일부 2025년 밴드에는 유효 픽셀값 65535가 존재한다. 전체 밴드 통계는 CSV에 있다.
- [Copernicus Browser 관련 공식 커뮤니티 답변](https://forum.dataspace.copernicus.eu/t/issues-with-download-image-analytical/3566)은 Visualise 탭 16-bit TIFF 값이 0–65535 범위로 스케일된다고 설명한다. [Sentinel Hub Evalscript 문서](https://docs.sentinel-hub.com/api/latest/evalscript/v3/)에는 반사도 ×65535의 UINT16 출력 예와 DN을 그대로 UINT16으로 출력하는 예가 모두 있다. [CDSE L2A 예시](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Process/Examples/S2L2A.html)도 두 출력 방식을 구분한다.
- 따라서 **반사도 ×65535, 추가 배율·오프셋 없음**으로 저장했다는 내보내기 설정 또는 evalscript가 확인될 때에만 `/65535`를 역변환으로 사용할 수 있다. 현재 TIFF의 태그와 파일명만으로는 이 조건이나 밴드별 추가 정규화 여부를 입증할 수 없다. `/10000`은 SAFE DN에 흔히 쓰이는 별도 변환을 가정한 **민감도 비교용**이며, 이 Browser TIFF에 검증된 변환이 아니다. 원본 SAFE의 보정 오프셋도 이 비교에는 적용하지 않았다.

### 90개 점 비교 결과

Week03에 저장된 지수를 원본 TIFF의 `/65535` 가정으로 재계산하면 최대 절대 오차가 NDVI 4.1×10⁻⁸, EVI 4.5×10⁻⁸, NDMI 4.6×10⁻⁸이다. 이는 Week03 계산의 재현성을 확인하지만 변환의 물리적 정확성을 증명하지는 않는다.

`/10000`을 **가정**하면 NDVI·NDMI는 공통 배율이 비율에서 상쇄되어 90개 점에서 변하지 않는다. EVI는 분모의 `+1` 항 때문에 평균이 0.511에서 1.335로 바뀐다(평균 차이 +0.824). 오프셋이나 밴드별 배율이 다르면 NDVI·NDMI도 변할 수 있다.

**판단:** 현재 증거만으로 지수를 새 변환값으로 교체하지 않는다. `/65535`는 잠정값으로 표시하고, Browser 내보내기 설정·evalscript 또는 동일 위치의 검증된 FLOAT32 반사도를 확보한 뒤 변환을 확정해야 한다. 그 결과가 기존 가정과 다르면 EVI는 반드시 다시 계산하고, 오프셋·밴드별 배율 차이가 있으면 NDVI·NDMI도 다시 계산한다.

## Random Forest AGBD 베이스라인

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe train_baseline.py
```

`train_baseline.py`는 최종 90행 CSV의 `ndvi`, `evi`, `ndmi`를 입력으로, `agbd`(Mg/ha)를 목표로 사용한다. `RandomForestRegressor(n_estimators=100, random_state=42)`를 두 분할에서 각각 새로 학습한다. 학습 전 전체 데이터로 표준화·결측 대체·특성 선택·튜닝을 하지 않는다. 검증 지표는 각 모델이 학습에 사용하지 않은 점에서만 계산한다.

### 데이터 점검

- 필수 입력·목표·좌표 결측치 0개, `point_id`·GEDI shot·좌표 중복 0개. 식생지수와 AGBD의 범위·평균·사분위수 등은 `data_diagnostics.json`에 기록한다.
- IQR 1.5배 기준으로 AGBD 1개(상한 약 290.55 Mg/ha), NDVI 7개, EVI 6개, NDMI 11개가 이상치 **후보**다. 자동 삭제하지 않았다. 이 기준은 분포 표시용이며 측정 오류 판정이 아니다.
- 90개 점은 모두 하나의 GEDI 파일에서 왔다. 최근접 점 거리는 중앙값 약 56.6m이고, UTM 52N 동서 약 5.72km·남북 약 2.97km 범위에서 여러 선형 묶음으로 분포한다. 금정산 전체의 임의 표본으로 볼 수 없다.

### 검증 설계와 결과

무작위 검증은 `train_test_split(test_size=0.2, random_state=42)`을 사용한다. 공간 검증은 관측점의 긴 동서 축을 따라 폭 1km의 easting 블록을 만들고, 가장 동쪽 블록 20개를 검증에 사용한다. 다른 블록의 점 중 검증점에서 500m 이내인 15개는 학습에서 제외했다. 블록과 버퍼는 **좌표만으로** 정하며 AGBD에 따라 바꾸지 않는다. 학습·검증 ID 교집합을 검사한다.

| 검증 방법 | 학습 | 검증 | 최소 학습·검증 거리 | MAE (Mg/ha) | RMSE (Mg/ha) | R² |
|---|---:|---:|---:|---:|---:|---:|
| 무작위 80/20 | 72 | 18 | 55.2m | 42.87 | 52.34 | -0.139 |
| 동쪽 1km 블록 + 500m 버퍼 | 55 | 20 | 500.6m | 57.71 | 70.61 | -2.165 |

공간 검증의 MAE는 무작위 검증보다 14.83 Mg/ha, RMSE는 18.27 Mg/ha 높고 R²는 2.026 낮다. 검증 점과 학습 점이 무작위 분할에서는 매우 가까우므로 공간 자기상관에 따른 낙관적 평가 가능성이 있다. 다만 두 검증은 표본 구성과 학습 크기도 다르므로 성능 차이 전체를 데이터 누수 하나로 설명할 수 없다. 두 R² 모두 음수로, 각 검증 집합의 평균 AGBD를 예측하는 기준보다 성능이 낮다. 공간 홀드아웃도 단일 동쪽 블록 한 번의 평가이며 다른 지역·방향의 성능을 대표하지 않는다. 표의 수치는 Week04 가상환경의 scikit-learn 1.9.1에서 실행한 결과다.

결과는 [`outputs/baseline_model/`](outputs/baseline_model/)의 `model_metrics.csv`, `actual_vs_predicted.png`, `feature_importance.png`, `validation_comparison.png`에 저장한다. 추가로 `data_diagnostics.json`, `validation_predictions.csv`, `spatial_split.png`를 남겨 자료 상태와 분할을 확인할 수 있게 했다. 변수 중요도는 각 **학습 집합**에서 얻은 impurity 기반 값이며 인과 효과나 안정적인 변수 순위를 뜻하지 않는다.

이 결과는 GEDI L4A의 **기존 모델 추정 AGBD**에 대한 탐색적 예측 점수다. 독립적인 현장 실측 검증 성능으로 해석하지 않는다. 반사도 `/65535` 변환도 아직 확정되지 않았으므로 특히 EVI와 이를 이용한 모델 점수에 불확실성이 있다. 90개 점의 점수로 금정산 전체 예측 정확도를 일반화할 수 없다.

## 학습 평균 Dummy와 동일 분할 비교

```powershell
.\.venv\Scripts\python.exe compare_baselines.py
```

`compare_baselines.py`는 `train_baseline.py`의 무작위 분할 규칙과 공간 분할 함수를 그대로 호출한다. 각 분할의 동일한 학습·검증 인덱스로 RF와 `DummyRegressor(strategy="mean")`을 각각 학습·평가한다. 공간 버퍼로 제외된 15개는 두 모델에 모두 사용하지 않는다. 기존 `baseline_model/`의 지표, 검증점 ID·순서·RF 예측값 및 공간 분할 ID를 대조한 결과 **기존 RF 성능이 1e-8 허용오차 안에서 재현**되었다. 기존 결과 파일은 변경하지 않는다.

| 검증 | 모델 | MAE (Mg/ha) | RMSE (Mg/ha) | R² |
|---|---|---:|---:|---:|
| 무작위 80/20 | RF | 42.87 | 52.34 | -0.139 |
| 무작위 80/20 | 학습 평균 Dummy | 43.34 | 50.94 | -0.080 |
| 공간 블록·500m 버퍼 | RF | 57.71 | 70.61 | -2.165 |
| 공간 블록·500m 버퍼 | 학습 평균 Dummy | 48.07 | 57.16 | -1.074 |

무작위 분할에서 RF는 MAE만 0.47 Mg/ha 낮고 RMSE는 1.39 Mg/ha 높다. 공간 분리에서는 RF의 MAE가 9.64 Mg/ha, RMSE가 13.45 Mg/ha 더 높다. **현재 한 번씩의 분할에서는 RF가 단순 기준보다 일관되게 낫다는 근거가 없다.** 분할별 학습·검증 AGBD 평균 차이, 한 궤도에 몰린 90개 점, 세 식생지수만 사용한 입력, 잠정적 반사도 변환이 가능한 한계다. 이 비교만으로 통계적 유의성이나 원인을 확정할 수 없다.

[`outputs/baseline_comparison/`](outputs/baseline_comparison/)의 `comparison_metrics.csv`에 네 평가를, `comparison_chart.png`에 **동일 검증점의 실제값과 두 모델 예측값을 같은 축**으로 표시한 산점도 및 지표 그래프를, `comparison_summary.md`에 재현 확인과 해석을 저장한다. GEDI L4A AGBD는 기존 모델 추정값으로 현장 실측값이 아니다.
