# `/65535` 변환 타당성 판정

## 판정: **C — 근거 부족으로 판단 불가**

| 구분 | 필요한 조건 | 현재 상태 |
|---|---|---|
| **A: 타당** | 이 2025 TIFF가 추가 보정·offset·비선형 처리 없이 **밴드 반사율×65535**로 기록됐다는 실제 요청·evalscript 또는 동일 픽셀의 독립 FLOAT32 반사율 일치 | 미확보 |
| **B: 잘못됨** | 실제 출력식·제품 메타데이터·독립 반사율이 `/65535`와 모순됨을 동일 장면·픽셀에서 입증 | 미확보 |
| **C: 판단 불가** | 저장된 지수 계산은 재현되지만 실제 출력식은 미확인 | **해당** |

2025 TIFF 메타데이터는 `scale=1`, `offset=0`, `nodata=null`이며 단위·evalscript 태그가 없다(`week05/step2/audit_metadata.json`). 2025 밴드 B08에는 32767을 넘는 유효 값 32,535개가 있어 이 TIFF를 원본 15-bit DN과 동일시하기 어렵다(`week05/step2/band_quality.csv`). [공식 L2A 원제품 명세](https://sentinels.copernicus.eu/documents/d/sentinel/s2-pdgs-cs-di-psd-v15-0)의 `(DN + BOA_ADD_OFFSET)/QUANTIFICATION_VALUE`를 Browser 출력 TIFF에 바로 적용할 근거도 없다. `/10000`은 대체 정답으로 판정하지 않는다.

[공식 Evalscript 문서](https://docs.sentinel-hub.com/api/latest/evalscript/v3/)는 UINT16이 출력 자료형이며 반환값을 정수 범위로 반올림·제한한다고 설명한다. [공식 L2A 예제](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Process/Examples/S2L2A.html)는 `REFLECTANCE ×65535` 출력을 **명시적으로 작성한 경우**를 보여준다. 현재 TIFF의 생성식이 그 예제와 같다는 증거는 없다. TIFF에 보이는 65535 값은 클리핑 가능성을 시사하지만 원인은 판정할 수 없다.

## 모델에 미칠 수 있는 영향 — 추정과 확인 사실 분리

| 대상 | 확인 사실 | 가능한 영향·한계 |
|---|---|---|
| 밴드값 | 현재 모델은 원시 밴드를 직접 사용하지 않고 NDVI·EVI·NDMI를 사용한다(`week04/train_baseline.py:27-32`). | 진짜 출력 변환에 밴드별 배율·offset·비선형 처리·클리핑이 있으면 지수의 입력값이 잘못됐을 수 있다. 실제 변환식이 없어 크기는 계산 불가. |
| NDVI·NDMI | 모든 관련 밴드에 동일한 **양의 상수 배율**, offset 0만 적용한다면 비율식에서 배율이 소거된다(`week03/src/vegetation_indices.py:15-27`). | 밴드별 offset 또는 서로 다른 배율은 이 불변성을 깨뜨린다. 비선형 출력·클리핑·NoData 누락도 별도 문제다. |
| EVI | 분모에 `+1`이 있어 공통 배율에도 값이 달라진다. Week04의 `/65535`와 `/10000` **가정 비교**에서 90점 EVI 평균은 0.511 대 1.335였다(`week04/outputs/reflectance_validation/index_comparison_statistics.csv`). | 이 차이는 민감도를 보여줄 뿐 어느 쪽이 참인지 말하지 않는다. 픽셀별 순서도 바뀔 수 있어 RF 분할에 영향을 줄 수 있다. |
| RF 학습·AGBD 예측 | Week04 RF는 잠정 지수로 학습됐고 Random/Spatial R²는 모두 음수였다(`week05/step1/outputs/model_comparison.csv`). | 실제 보정값이 다르면 학습 경계·예측·성능이 달라질 수 있다. 수치적 효과·성능 개선 방향은 검증 반사율과 재학습 없이는 알 수 없다. 이번에는 재학습하지 않았다. |

단순 상수 배율을 밴드값에 적용하는 문제와 **밴드별 offset/비선형 출력/NoData 처리**는 구분해야 한다. 같은 배율만 달랐다면 NDVI·NDMI는 유지될 수 있으나, 그 사실만으로 전체 입력이 검증된 것은 아니다. 기존 90점의 좌표·마스크·SCL 점검은 통과했으나 방사 값의 물리적 정확성을 입증하지는 않는다(`week05/step2/report.md`).
