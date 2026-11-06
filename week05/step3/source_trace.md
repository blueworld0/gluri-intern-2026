# Week05 Step3 — 원본 출력 설정 추적

## 탐색 범위와 확인 결과

저장소의 Week03–Week05 코드·문서·설정·산출물과 현재 Git 이력(4개 커밋)을 확인했다. 검색 대상에는 `evalscript`, `sampleType`, `harmonizeValues`, Sentinel Hub Process API 엔드포인트와 TIFF 생성 관련 문자열을 포함했다. **아래 '미확인'은 현재 저장소에서 찾지 못했다는 뜻**이며, 사용자의 Browser 기록이나 외부 계정에도 없다는 뜻은 아니다.

| 항목 | 현재 저장소에서 확인된 내용 | 근거 |
|---|---|---|
| TIFF 다운로드/생성 코드 | 2025·2026 **원본** `_(Raw).tiff`를 생성하거나 다운로드하는 코드 미확인. `week03/main.py`와 `week03/analyze_sentinel_2025.py`는 이미 존재하는 TIFF를 읽고 **지수 결과 TIFF**를 저장한다. | `week03/main.py:112-153`, `week03/analyze_sentinel_2025.py:269-315`; 저장소 코드 검색 |
| API 요청 코드·evalscript | 현재 저장소에서 미확인. 문서에 evalscript **일반 설명 링크**는 있지만 이 파일을 만든 요청 본문은 없다. | `week04/validate_reflectance.py:27-31`; 저장소 코드·설정 검색 |
| `sampleType`, `harmonizeValues`, `output`, 입력 단위 | 실제 다운로드 요청값 미확인. TIFF 첫 밴드가 `uint16`인 것은 관찰 사실이며 `sampleType: UINT16` 요청의 기록은 아니다. `dataMask` 형태의 둘째 밴드가 있으나 반환식은 미확인. | `week05/step2/audit_metadata.json`; `week04/outputs/reflectance_validation/tiff_metadata.json` |
| 출처 | **2026 TIFF**의 Copernicus Browser `Visualise → Download Image → Analytical → TIFF 16-bit → Raw Bands, dataMask ON` 사용자 확인 기록과 파일 SHA-256 설정이 있다. **2025 TIFF**는 유사한 파일명·구조에 근거해 같은 유형이라고 잠정 취급했다. 2025 다운로드 로그나 별도 설정 확인 기록은 현재 저장소에 없다. | `week03/reflectance_config.json:2-17`, `week03/README_SENTINEL-2.md:9`, `week03/README_SENTINEL_2025.md:13-16` |
| 데이터 컬렉션·제품 | 2025 TIFF 파일명에는 `Sentinel-2_L2A`와 `2025-04-28 00:00–23:59` 요청 구간이 있다. SCL JP2 이름에는 `T52SDE_20250428T020711`이 있고 기존 기록에는 사용자 확인 `S2C_MSIL2A`, 처리 기준 `N05.11`이 있다. TIFF 자체의 제품 ID·타일·처리 기준 태그는 없다. | `week03/Sentinel-2/2025-04-28/`, `week03/outputs/sentinel_2025/scl_validation/scl_quality_report.json`의 `product`, `tiff_association` |
| TIFF 생성 날짜 | TIFF 내부 생성 시각 태그와 요청 로그를 확인하지 못했다. 파일시스템 생성·수정 시각은 다운로드/복사 때 바뀔 수 있어 실제 출력 생성일의 증거로 사용하지 않았다. | `week05/step2/audit_metadata.json`의 `dataset_tags`, `band_tags` |
| 실제 촬영 날짜 | SCL **파일명**은 2025-04-28 02:07:11을 나타낸다. TIFF 이름은 그날의 **요청 범위**만 표시한다. TIFF와 SCL이 동일 장면이라는 독립 증거가 없어 TIFF의 정확한 촬영 시각은 미확인이다. | `week03/README_SCL_2025.md:11`, `week03/README_SENTINEL_2025.md:13`, `scl_quality_report.json`의 `tiff_association` |

2025 밴드는 B02, B03, B04, B08, B11, B12이고 모델에 들어간 NDVI·EVI·NDMI는 이 가운데 B02·B04·B08·B11로 계산됐다(`week03/src/vegetation_indices.py:9-39`, `week04/train_baseline.py:27-32`). 2025 코드의 `SCALE = 1/65535`는 명시적으로 잠정 가정이다(`week03/analyze_sentinel_2025.py:3-4,36`).

## 문서의 적용 범위

[Sentinel Hub Evalscript 공식 문서](https://docs.sentinel-hub.com/api/latest/evalscript/v3/)는 `setup`의 입력 단위·출력 `sampleType` 및 `evaluatePixel` 반환식이 실제 TIFF 값을 결정한다고 설명한다. [Sentinel-2 L2A 문서](https://docs.sentinel-hub.com/api/latest/data/sentinel-2-l2a/)는 `REFLECTANCE`·`DN`, `harmonizeValues` 선택에 따라 값의 의미가 달라짐을 설명한다. 두 문서는 **가능한 출력 방식**을 정의하지만, 현재 2025 TIFF에 어떤 설정을 썼는지는 입증하지 않는다.
