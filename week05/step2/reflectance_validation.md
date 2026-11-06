# 반사율 변환 검증

## 판정

**`/65535`는 이 파일에서 사용한 계산을 재현하지만, 올바른 반사율 변환으로 검증되지 않았다.** `uint16`은 출력 자료형의 범위일 뿐 저장값이 물리 반사율×65535임을 보장하지 않는다. TIFF의 scale=1, offset=0, NoData 미설정, 단위·밴드 설명 태그 없음은 영상값에 적용할 반사율 계수를 제공하지 않는다. 근거: `audit_metadata.json`, `week04/outputs/reflectance_validation/tiff_metadata.json`.

## 제품 문서와 이 TIFF의 차이

- [Copernicus Sentinel-2 L2A 제품 명세](https://sentinels.copernicus.eu/documents/d/sentinel/s2-pdgs-cs-di-psd-v15-0)와 [L2A 품질 보고서](https://sentinels.copernicus.eu/documents/d/sentinel/ompc-cs-dqr-002-03-2023-i60r0-msi-l2a-dqr-april-2023)는 **원제품 DN**의 BOA 반사율을 `(DN + BOA_ADD_OFFSET) / QUANTIFICATION_VALUE`로 산출한다고 명시한다. 기준 04.00 이후 offset 항이 있으며 DN=0은 NoData이다. 이 규칙을 Browser가 다시 내보낸 TIFF에 바로 적용할 수는 없다.
- [Sentinel Hub L2A 자료 문서](https://docs.sentinel-hub.com/api/latest/data/sentinel-2-l2a/)는 입력 단위 `REFLECTANCE`와 `DN`, 기본 harmonizeValues 설정을 구분한다. harmonized DN과 원제품 DN의 의미도 다르다.
- [Evalscript 문서](https://docs.sentinel-hub.com/api/latest/evalscript/v3/)에서 UINT16은 **출력 형식**이다. evalscript 반환값은 정수로 반올림·범위 제한되며 65535로 자동 정규화되지 않는다. [공식 Process API L2A 예제](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Process/Examples/S2L2A.html)는 반사율에 65535를 곱해 UINT16로 쓰는 **명시적** 사례를 보여준다. [CDSE 포럼의 Browser 담당자 설명](https://forum.dataspace.copernicus.eu/t/issues-with-download-image-analytical/3566)은 Visualise의 16-bit Analytical 출력이 0–65535 범위로 조정되고 원본 JP2 DN과 절대값을 비교하면 안 된다고 설명하지만, 현재 여섯 TIFF에 적용된 정확한 evalscript·배율·offset은 제시하지 않는다.

따라서 `/10000`도 현재 TIFF의 확인된 대체 계수가 아니다. 원제품 DN 공식에 필요한 제품 XML의 `BOA_ADD_OFFSET`·`QUANTIFICATION_VALUE`가 없고, TIFF가 원제품 DN인지도 확인되지 않았다. 2025 B08 유효 픽셀 중 **32,535개가 32767 초과**, 6개가 65535인 사실은 이 TIFF를 원제품의 15-bit DN으로 취급해서는 안 된다는 추가 근거다(`band_quality.csv`).

## 데이터 자체에서 확인되는 것

| 밴드 | 원시 최솟값 | 최댓값 | 평균 | 중앙값 | 유효 65535 | `/65535` 가정 평균 |
|---|---:|---:|---:|---:|---:|---:|
| B02 | 0 | 65535 | 3159.23 | 1848 | 8 | 0.0482 |
| B03 | 0 | 65535 | 5546.82 | 4653 | 8 | 0.0846 |
| B04 | 26 | 65535 | 4332.80 | 2549 | 11 | 0.0661 |
| B08 | 13 | 65535 | 21687.95 | 21830 | 6 | 0.3309 |
| B11 | 944 | 54387 | 12809.04 | 12412 | 0 | 0.1955 |
| B12 | 596 | 58477 | 7393.22 | 6141 | 0 | 0.1128 |

유효 픽셀은 밴드마다 909,212개, dataMask=0은 152개다. 위 반사율 숫자는 **가정 결과**이다. 값이 주로 0–1에 든다는 것만으로 물리 단위가 입증되지 않는다. 공식 Sentinel Hub L2A 문서는 반사율이 1을 넘을 수도 있다고 명시한다. 65535 값은 출력 상한에 걸렸을 가능성이 있으나 그 원인(밝은 물체, evalscript 증폭, 원래 포화 등)은 식별되지 않는다. 최종 90점에는 65535가 없고 원시 0이 한 점의 한 밴드에 있다. dataMask가 유효하므로 원시 0을 무조건 NoData라고 볼 수 없으며, 해당 픽셀을 별도 확인할 필요가 있다.

## 식생지수에 대한 영향

동일한 **상수 배율만** 바꾸고 offset이 0이면 NDVI와 NDMI는 대수적으로 불변이다. EVI는 분모의 `+1` 때문에 불변이 아니다. Week04의 90점 민감도 비교에서 EVI 평균은 `/65535` 가정 **0.511**, `/10000` 가정 **1.335**였다(`week04/outputs/reflectance_validation/index_comparison_statistics.csv`). 이는 `/10000`이 정답이라는 뜻이 아니다.

**밴드별 offset**이 존재하거나 밴드마다 배율이 다르면 NDVI·NDMI도 달라진다. **NoData·dataMask 오류**는 배율과 별개로 잘못된 픽셀을 유효값으로 넣는다. 이번 감사에서 90점의 마스크 오류는 발견되지 않았지만, 원제품 offset과 Browser 출력 과정을 확인할 자료가 없어 반사율 정확도는 **검증 불가**다. Week04의 원시 TIFF 재계산과 저장 지수 일치(`week04/outputs/reflectance_validation/candidate_index_comparison.csv`)는 계산 재현성만 입증한다.
