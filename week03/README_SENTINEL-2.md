# GLURI R&D Week03: Sentinel-2 입력 데이터 준비

`Sentinel-2/`의 GeoTIFF를 검사·전처리하고 NDVI, EVI, NDMI를 계산합니다. 원본 TIFF는 변경하지 않습니다. 이번 단계에서는 GEDI, DEM, 임상도, AGB 예측을 사용하지 않습니다.

## 현재 데이터와 반사도 변환

확인된 TIFF 6개는 B02, B03, B04, B08, B11, B12입니다. 누락된 예상 밴드는 없습니다. 모두 578 × 834 픽셀, EPSG:32652, 약 9.9975 × 9.9941 m이며 같은 공간 격자입니다. 각 파일의 1번 밴드는 영상값이고, 2번 밴드는 0/1의 `dataMask` 형태입니다. 현재 마스크 픽셀은 모두 1입니다.

사용자가 여섯 파일을 같은 설정으로 Copernicus Browser의 **Visualise → Download Image → Analytical → TIFF (16-bit) → Raw Bands, dataMask ON**에서 받았고 촬영일은 2026-09-12라고 확인했습니다. 사용자 제공 변환 계수는 `reflectance = pixel / 65535`(offset 0)입니다. 이 설정과 파일별 SHA-256 해시를 [reflectance_config.json](reflectance_config.json)에 기록했습니다. `python main.py`는 현재 파일이 설정과 일치할 때 이 계수를 적용합니다. 새 파일이나 변경된 파일이 있으면 계수를 자동 적용하지 않습니다. GeoTIFF 자체에는 물리 단위나 반사도 계수가 없으므로 이 변환의 근거는 확인된 내보내기 설정입니다. [Copernicus Browser 포럼 설명](https://forum.dataspace.copernicus.eu/t/issues-with-download-image-analytical/3566)은 Visualise 탭의 16비트 다운로드가 값 범위를 0–65535로 조정한다고 설명합니다.

## 설치 및 실행

Python 3.10 이상을 사용합니다. 현재 프로젝트의 `.venv`가 있으면 재사용합니다.

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py
```

가상환경이 없다면:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py
```

새 데이터에 대해 계수를 별도로 확인했다면 `--reflectance-scale`과 `--reflectance-offset`으로 지정할 수 있습니다. 이 옵션은 모든 입력 밴드에 동일한 계수를 적용합니다. 촬영일을 따로 확인했다면 `--acquisition-date YYYY-MM-DD`를 사용합니다. 파일명과 다른 날짜는 거부합니다.

```powershell
.\.venv\Scripts\python.exe main.py --reflectance-scale 0.000015259021896696 --reflectance-offset 0 --acquisition-date 2026-09-12
```

## 결과물

- `outputs/rasters/ndvi.tif`, `evi.tif`, `ndmi.tif`: FLOAT32, EPSG:32652, 원본 격자, NoData=-9999.
- `outputs/figures/rgb.png`, `ndvi.png`, `evi.png`, `ndmi.png`, `ndvi_histogram.png`: 촬영일과 좌표가 포함된 시각화. RGB는 채널별 원시값 2–98 백분위 스트레칭 미리보기입니다.
- `outputs/statistics.csv`: 지수별 유효 픽셀 수, 최소·최대, 평균, 표준편차, 중앙값.

지수 지도 색상은 유효값의 2–98 백분위 구간으로 표시하지만, 제목은 실제 데이터의 최소·최대값을 표기합니다. 극단값도 GeoTIFF와 통계에 그대로 보존됩니다. RGB는 채널별 스트레칭이므로 단일 수치 Colorbar 대신 채널과 표시 범위를 적었습니다.

## 품질 확인과 해석 주의

반사도 변환 후 모든 값은 0–1 범위입니다. B02와 B03에는 같은 위치(row 499, col 341)의 `65535` 포화 픽셀이 각각 1개 있습니다. NDVI와 NDMI는 모두 [-1, 1] 범위입니다. EVI는 50개 픽셀이 [-1, 1] 밖이고, 그중 12개는 절댓값 2를 초과합니다. EVI 최대값 296.513은 row 780, col 536에서 분모가 약 0.00203으로 작아져 발생했습니다. 최소값 -173.422는 row 403, col 249에서 분모가 약 -0.00237입니다. 두 픽셀 모두 포화 픽셀과 다른 위치입니다. 원식을 따른 결과이므로 값을 임의로 자르거나 삭제하지 않았습니다. **향후 AGB 입력 특성으로 EVI를 사용할 때는 비식생·밝은 표면 및 작은 분모 영역에 대한 품질 필터를 결정해야 합니다.**

NaN, NoData, GDAL 마스크, 명시적 `dataMask`를 합쳐 유효 픽셀을 판별합니다. 다른 격자는 `rasterio.warp.reproject`로 B04 격자에 정렬합니다. 날짜 표기가 서로 다르면 결합을 중단합니다. 이후 GEDI와 결합할 때는 실제 원본 제품의 촬영 시각 및 구름/그림자 품질 정보를 추가로 확보하는 것이 좋습니다.
