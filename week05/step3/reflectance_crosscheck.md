# 반사율 교차 검증 가능성 검토

## 우선순위별 조사

| 우선 | 경로 | 확보된 증거 | 직접 수치 비교 가능 여부 |
|---:|---|---|---|
| 1 | 기존 API 요청·evalscript 복원 | 현재 저장소에 2025 요청 코드·본문·evalscript가 없다. 2026 Browser 사용 설명과 해시 설정은 2025 출력 설정을 특정하지 못한다. | **불가**. `units`, `sampleType`, `harmonizeValues`, 반환식, 모자이크·리샘플링·출력 격자 필요. |
| 2 | 원본 TIFF 메타데이터 | 여섯 2025 TIFF는 동일 EPSG:32652 격자, 첫 밴드 `uint16`, 둘째 0/1 마스크다. `scale=1`, `offset=0`, `nodata=null`, 단위·밴드 태그·제품 ID·정확한 촬영 시각 없음. | **불가**. 1/0은 반사율 변환식을 뜻하지 않는다. `week05/step2/audit_metadata.json` 참조. |
| 3 | 동일 장면의 공식 제품 메타데이터 | SCL JP2 파일명과 사용자 확인 제품 정보가 있지만 원본 SAFE/제품 XML(`MTD_MSIL2A.xml`)과 datastrip 메타데이터가 현재 저장소에 없다. TIFF의 장면 ID도 없다. | **불가**. 제품 XML을 확보해도 Browser TIFF와 동일 장면·출력 과정을 먼저 연결해야 한다. |
| 4 | 검증된 FLOAT32 밴드 반사율 | 저장소의 `float32` GeoTIFF는 `/65535`를 사용해 만든 NDVI·EVI·NDMI **지수 결과**다. 독립 밴드 반사율이 아니다. 2026 TIFF는 다른 시점·범위여서 검증 기준이 될 수 없다. | **불가**. 같은 장면·밴드·위치·격자 또는 정확히 알려진 재표본화가 필요. |

## 교차 비교를 하지 않은 이유

Week04는 2025 원시 TIFF에 `/65535`를 적용하면 저장된 지수가 거의 완벽하게 재현됨을 보였다(`week04/outputs/reflectance_validation/candidate_index_comparison.csv`). 이는 **동일한 가정으로 계산을 반복한 내부 재현성 검사**다. 독립 반사율과의 교차 비교가 아니다. `/10000` 시나리오도 검증 자료가 아니므로 Step3 수치 비교 CSV·그래프를 추가하지 않았다.

SCL JP2의 촬영 시각을 근거로 2025 TIFF와 다른 영상을 곧바로 비교하면 시점·장면 차이가 변환 차이와 섞인다. 현재 2026 밴드와 2025 밴드도 동일 시점이 아니다. 이번에는 데이터 API·브라우저 계정 요청을 실행하지 않았다.

## 검증 가능한 비교 절차

1. 2025 원본 TIFF의 **실제 생성 요청**을 확보한다: 장면/product ID, AOI, timeRange, `input.data.type`, `units`, `harmonizeValues`, mosaicking, resampling, 출력 CRS·크기·변환, `sampleType`, evalscript `evaluatePixel`, `dataMask` 처리. 설정값을 문서의 기본값으로 추정해 채우지 않는다.
2. 요청의 장면 ID를 SCL JP2 및 공식 SAFE/제품 XML과 대조한다. 원제품의 `BOA_ADD_OFFSET`·`QUANTIFICATION_VALUE`를 읽더라도, Browser 출력이 원제품 DN인지 먼저 확인한다. [Copernicus 제품 명세](https://sentinels.copernicus.eu/documents/d/sentinel/s2-pdgs-cs-di-psd-v15-0)의 DN 공식은 출력 TIFF에 자동 적용되지 않는다.
3. 동일 장면에서 **검증된 `REFLECTANCE` FLOAT32 밴드**를 확보한다. 가능하면 동일 2025 TIFF 격자와 `dataMask`로 요청한다. 다르면 재표본화 규칙을 문서화하고 픽셀 중심·밴드별 정합을 검증한다.
4. B02/B03/B04/B08/B11/B12의 공통 유효 픽셀과 최종 GEDI 90점에서 `TIFF/65535 − REFLECTANCE`를 밴드별로 계산한다. MAE, 최대 절대차, 기울기·절편, 값 0·65535 주변의 편향과 마스크 불일치를 보고한다. 같은 촬영 시점·제품·픽셀임이 확인되지 않으면 **직접 수치 검증으로 분류하지 않는다**.

## API 실행 조건

현재 저장소에는 실행 가능한 요청 설정이나 인증 정보가 확인되지 않았다. 임의 인증 정보를 만들지 않았다. [CDSE Sentinel Hub 처리 단위 문서](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Overview/ProcessingUnit.html)에 따르면 성공한 Process API 요청은 처리 단위(PU)를 차감하며 FLOAT32 TIFF에는 출력 형식 계수가 적용된다. 따라서 데이터 API 호출은 하지 않았고, 위 절차만 기록했다. 향후 호출 전 계정 권한과 예상 PU 사용량을 확인해야 한다.
