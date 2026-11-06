# Sentinel-2 입력 처리 과정 추적

## 자료 계보

| 단계 | 확인 내용 | 근거 |
|---|---|---|
| 내보내기 출처 | **2026 TIFF**에 대해 사용자가 Copernicus Browser **Visualise → Download Image → Analytical → TIFF 16-bit → Raw Bands**, `dataMask` 포함을 확인했다. **2025 TIFF**는 파일명·동일 형식을 근거로 같은 유형의 출력으로 취급했지만 해당 요청 로그·evalscript나 별도 사용자 확인 기록은 없다. | `week03/reflectance_config.json:2`, `week03/README_SENTINEL-2.md:9`, `week03/README_SENTINEL_2025.md:14` |
| 제품 | TIFF 파일명은 `Sentinel-2_L2A`이고 2025-04-28 요청 구간을 표시한다. SCL JP2 파일명은 `T52SDE_20250428T020711_SCL_20m.jp2`이다. `S2C_MSIL2A`, 타일 T52SDE, 처리 기준 N05.11은 사용자 확인 기록이며 TIFF 내부에서 제품 ID·기준을 독립 확인할 수 없다. | `week03/Sentinel-2/2025-04-28/`, `week03/outputs/sentinel_2025/scl_validation/scl_quality_report.json`의 `product`, `tiff_association` |
| 원본 밴드 | B02, B03, B04, B08, B11, B12. 각 TIFF는 첫 밴드 `uint16` 값, 둘째 밴드 0/1 `dataMask`. 2025 TIFF 6개가 같은 EPSG:32652, 1244×731 격자를 갖는다. | `week05/step2/audit_metadata.json`, `week05/step2/band_quality.csv` |
| 반사율 가정 | 2025 처리에서 `SCALE = 1/65535`, offset 0을 잠정값으로 지정했다. 실제 계산은 유효 픽셀 마스킹·필요한 격자 정렬 후 `original * scale + offset`이다. 2026 설정 파일의 동일 형식 가정이 이유이며 2025 TIFF의 검증된 계수는 아니다. | `week03/analyze_sentinel_2025.py:3-4,36,276`, `week03/src/preprocess.py:45-113`, `week03/README_SENTINEL_2025.md:13-16` |
| 식생지수 | NDVI=(B08−B04)/(B08+B04), EVI=2.5(B08−B04)/(B08+6B04−7.5B02+1), NDMI=(B08−B11)/(B08+B11). 2025 지수는 `float32` GeoTIFF, NoData −9999로 저장했다. B03·B12는 이 세 지수에 사용하지 않았다. | `week03/src/vegetation_indices.py:9-39`, `week03/main.py:22-34`, `week03/analyze_sentinel_2025.py:310-315` |
| GEDI 결합 | GEDI 위경도(EPSG:4326)를 영상 CRS로 변환해 해당 지수 raster의 행·열을 찾아 단일 픽셀값을 추출했다. `point_id`, 중복 픽셀, 날짜 차이, 유효값을 기록했다. 표적 AGBD는 GEDI L4A 추정값이다. | `week03/merge_sentinel_gedi.py:96-173`, `week03/analyze_sentinel_2025.py:286-303` |
| SCL 선별 | 20m SCL을 최근접 방식으로 10m 지수 격자에 맞췄다. GEDI 좌표의 원래 SCL과 정렬된 지수 픽셀의 SCL이 모두 4(식생)이고 결합에 성공한 90점을 선택했다. | `week03/validate_scl_2025.py:91-115,184-205`, `week03/outputs/sentinel_2025/scl_validation/scl_quality_report.json` |
| Week04 입력 | 최종 90점 CSV의 `ndvi`, `evi`, `ndmi`가 RF 입력이다. | `week04/load_data.py:8-34`, `week04/train_baseline.py:27-32` |

## 출처와 시간의 한계

2025 TIFF 파일명에는 요청 날짜 범위만 있다. TIFF 메타데이터에는 정확한 촬영 시각·타일·제품 ID·처리 기준·evalscript가 없다(`audit_metadata.json`). SCL JP2 파일명과 사용자 확인 출처는 2025-04-28 촬영과 일치하지만, TIFF 여섯 개가 **그 SCL과 동일 제품의 출력인지**는 독립적으로 입증되지 않는다. GEDI 관측일은 2025-04-25 UTC이고, 저장된 **요청 날짜 기준** 차이는 3일이다. SCL 촬영 시각이 실제 TIFF에도 해당한다면 약 56.94시간이다. `week03/README_SENTINEL_2025.md:13`, `week03/README_SCL_2025.md:11`, `week03/outputs/sentinel_2025/scl_validation/scl_quality_report.json` 참고.

코드의 `project_and_sample(..., nominal_date)`는 요청 날짜를 `sentinel_acquisition_date` 컬럼에 기록한다(`week03/analyze_sentinel_2025.py:291`, `week03/merge_sentinel_gedi.py:139-143`). 값 자체가 틀렸다는 증거는 없지만, 컬럼명은 검증된 촬영 날짜처럼 보이므로 출처 상태를 함께 읽어야 한다.
