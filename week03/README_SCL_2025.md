# 2025 Sentinel-2 SCL 품질 검증

기존 2025 Sentinel-2 지수와 GEDI 97개 점을 유지한 채 `T52SDE_20250428T020711_SCL_20m.jp2`를 검사한다.

```powershell
.\.venv\Scripts\python.exe validate_scl_2025.py
```

기존 결과가 있으면 덮어쓰지 않고 중단한다. 출력은 `outputs/sentinel_2025/scl_validation/`에 저장한다. 기존 `requirements.txt`의 패키지만 사용한다.

JP2의 파일명에 있는 타일 `T52SDE`와 촬영 시각 `2025-04-28 02:07:11 UTC`, 20m 격자 및 EPSG:32652를 확인한다. TIFF 여섯 개의 요청 날짜·L2A 표기·CRS·공간 범위를 대조한다. TIFF에는 타일·정확한 촬영 시각·제품 ID 태그가 없으므로 같은 제품이라는 사실을 파일만으로 독립 검증할 수 없다. 사용자께서 제공한 제품 출처와 모순이 없다는 근거로 진행하며, 모순이 발견되면 중단한다.

SCL은 최근접 이웃으로 10m 지수 격자에 정렬한다. 각 GEDI 좌표에서 원래 20m SCL을 읽고, 대응하는 10m 지수 픽셀에서도 정렬된 SCL을 읽는다. 두 값이 모두 `4`(식생)인 점만 보수적인 분석 후보로 남긴다. 두 값이 다르면 `scl_boundary_disagreement`로 표시한다. `3`, `8`, `9`, `10`은 구름·그림자로 분류하며 다른 비식생·어두운 영역·미분류 클래스는 별도로 보고한다. 클래스 의미는 [Sentinel-2 L2A 자료 설명](https://docs.sentinel-hub.com/api/latest/data/sentinel-2-l2a/)을 따른다.

출력에는 모든 GEDI 점과 원래/정렬 SCL값(`gedi_scl_points.csv`), 최종 후보(`gedi_scl_valid_points.csv`), 클래스별 개수, 전후 상관계수, 지도 두 장과 JSON 품질 보고서가 있다. SCL=4는 식생 분류일 뿐 산림 경계를 증명하지 않으며 SCL 자체도 구름을 완벽히 판별하지 못한다.

SCL 분류값에는 반사도 변환을 적용하지 않는다. 기존 지수 래스터에 사용한 Analytical TIFF의 `÷65535` 계수는 여전히 잠정적이며, 원본 제품 DN의 변환 방식과 동일시하지 않는다. [Copernicus Browser Analytical 출력 설명](https://forum.dataspace.copernicus.eu/t/issues-with-download-image-analytical/3566)
