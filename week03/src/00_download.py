"""1단계: 공개 자료 내려받기 (금정산 주변만 잘라서).

- ESA CCI Biomass v4: 지상부 바이오매스(AGB, t/ha), 100 m, 2017~2020
- ESA WorldCover: 토지피복 10 m (나무 픽셀 고르기용), 2020·2021
- Copernicus DEM: 표고 30 m (임시 경계 만들 때 사용)
"""
from config import DATA, BBOX, CCI_YEARS, CCI_URL, WORLDCOVER, DEM_URL
from common import clip_remote

print("[ESA CCI Biomass v4]")
for y in CCI_YEARS:
    clip_remote(CCI_URL.format(y=y, sd=""), DATA / "cci_biomass_v4" / f"agb_{y}.tif", BBOX)
    clip_remote(CCI_URL.format(y=y, sd="_SD"), DATA / "cci_biomass_v4" / f"agb_sd_{y}.tif", BBOX)

print("[ESA WorldCover]")
for y, url in WORLDCOVER.items():
    clip_remote(url, DATA / "worldcover" / f"worldcover_{y}.tif", BBOX)

print("[Copernicus DEM]")
clip_remote(DEM_URL, DATA / "dem" / "cop_dem30.tif", BBOX)
print("완료")
