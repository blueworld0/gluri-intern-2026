# Week04: AGBD 예측 준비

현재 단계는 Week03의 SCL 품질 검증을 통과한 식생 후보를 읽고 데이터 상태를 확인한다. 예측 모델은 아직 포함하지 않는다.

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

스크립트는 전체 행 수, 컬럼 목록, AGBD·식생지수 기초 통계, 각 컬럼의 결측치 수를 출력한다. `outputs/`는 향후 Week04 결과용 폴더이며 현재 실행에서는 파일을 생성하지 않는다.

SCL 클래스 4는 식생 분류이며 산림 여부를 보증하지 않는다. Week03의 2025년 지수 계산에 사용된 반사도 계수는 잠정적이다.
