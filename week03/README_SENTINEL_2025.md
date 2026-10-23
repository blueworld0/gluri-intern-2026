# Sentinel-2 2025 영상과 GEDI 비교

`analyze_sentinel_2025.py`는 `Sentinel-2/2025-04-28/`의 여섯 밴드를 검사하고, 기존 `src/` 전처리·식생지수·시각화 함수와 `merge_sentinel_gedi.py`의 좌표 샘플링 함수를 재사용한다. 원본 TIFF, 기존 2026년 결과, GEDI 품질 필터는 수정하지 않는다.

```powershell
.\.venv\Scripts\python.exe analyze_sentinel_2025.py
```

기본 결과는 `outputs/sentinel_2025/`에 저장된다. 같은 폴더가 이미 있으면 덮어쓰지 않고 중단한다. 새 실행 결과가 필요하면 `--output-dir outputs/sentinel_2025_run2`처럼 새 경로를 지정한다. 기존 `requirements.txt` 패키지만 사용한다.

## 해석 조건

- 파일명은 `2025-04-28 00:00–23:59` 요청 구간을 표시한다. TIFF 내부에 실제 촬영 시각이 없어 **정확한 촬영일은 검증되지 않았다**. GEDI 관측일 `2025-04-25 UTC`와의 3일 차이는 파일명 기준의 잠정 시차다.
- TIFF의 scale/offset 메타데이터는 1/0이며 반사도 의미의 계수는 없다. 2026년 16비트 Copernicus Browser Analytical Raw 다운로드와 동일한 형식으로 보고 `reflectance = pixel / 65535`를 **잠정 적용**했다. 내보내기 설정 또는 원본 제품 메타데이터로 확인하기 전에는 이 값을 검증된 반사도로 취급하지 않는다. [Copernicus 포럼 설명](https://forum.dataspace.copernicus.eu/t/issues-with-download-image-analytical/3566)
- 두 번째 TIFF 밴드는 설명 태그가 없으나 0/1 값만 갖고 여섯 파일에 동일하게 들어 있어 dataMask로 판별했다. 0 픽셀은 계산에서 제외한다.
- SCL, CLD, CLM이 현재 제공되지 않아 구름과 구름 그림자는 평가 불가다. `quality_review_points.csv`는 구름 판정이 필요한 점도 모두 포함한다. 이 목록은 구름 영향이 확인된 점들의 목록이 아니다. SCL을 같은 입력 폴더에 추가하면 스크립트가 픽셀별 SCL 분류를 기록한다. [SCL 클래스 정의](https://docs.sentinel-hub.com/api/latest/data/sentinel-2-l2a/)
- EVI는 분모가 0에 가까우면 큰 값을 낸다. 전체 영상의 극단값과 97개 샘플값을 구분해 확인한다. 계산식 결과를 임의로 자르거나 클리핑하지 않는다.

## 출력

- `rasters/ndvi.tif`, `evi.tif`, `ndmi.tif`: float32, EPSG:32652, 원본 2025 격자, NoData=-9999. 요청 날짜와 계수의 검증 상태를 태그로 저장한다.
- `figures/`: RGB, 지수 지도·NDVI 히스토그램, GEDI 점을 표시한 RGB 품질 검토 지도, 같은 97개 점의 2025/2026 비교 그래프, AGBD 산점도 및 상관계수 비교.
- `gedi_sentinel_2025_merged.csv`: 2025년 샘플값과 위치·dataMask·포화·SCL 상태.
- `paired_2025_2026.csv`: 같은 GEDI 파일·BEAM·shot_number의 두 시기 값과 차이.
- `index_comparison.csv`: 두 시기에 모두 값이 있는 같은 점을 사용한 지수 분포 및 AGBD와의 Pearson·Spearman 상관계수.
- `quality_review_points.csv`: NoData, 포화 또는 미평가/의심 SCL 상태가 있는 점.
- `raster_statistics.csv`, `quality_report.json`: 전체 영상 통계, 공간 격자·품질·시차·변환 가정.

GEDI AGBD는 모델 추정치이며, 이 상관계수는 탐색적 통계다. 두 시기 중 어느 영상이 더 정확한 AGB 모델을 보장한다는 뜻은 아니다.
