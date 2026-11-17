"""1단계: 공개 자료 내려받기 (금정산 주변만 잘라서).

- ESA CCI Biomass v4 (2017–2020, NASA MAAP) / v6 (2015–2022, CEDA): 지상부 바이오매스 AGB, t/ha, 100 m
- ESA WorldCover: 토지피복 10 m (나무 픽셀 고르기용), 2020·2021
- Copernicus DEM: 표고 30 m (임시 경계 만들 때 사용)
이미 받은 파일은 건너뛰므로, 중간에 끊겨도 다시 실행하면 남은 것만 받습니다.
CEDA(영국) 서버는 느려서 v6 는 연도당 몇 분씩 걸릴 수 있습니다.
"""
from config import DATA, BBOX, CCI, WORLDCOVER, DEM_URL
from common import clip_remote, find_ceda_tifs

print("[ESA CCI Biomass v4 — NASA MAAP]")
v4 = CCI["v4"]
for y in v4["years"]:
    clip_remote(v4["url"].format(y=y, sd=""), v4["dir"] / f"agb_{y}.tif", BBOX)
    clip_remote(v4["url"].format(y=y, sd="_SD"), v4["dir"] / f"agb_sd_{y}.tif", BBOX)

print("[ESA CCI Biomass v6 — CEDA]")
v6 = CCI["v6"]
for y in v6["years"]:
    if (v6["dir"] / f"agb_{y}.tif").exists() and (v6["dir"] / f"agb_sd_{y}.tif").exists():
        print(f"  {y} 이미 있음")
        continue
    agb_url, sd_url = find_ceda_tifs(v6["url"].format(y=y))
    if not agb_url:
        print(f"  {y}: CEDA 에서 파일을 못 찾음 → 건너뜀")
        continue
    clip_remote(agb_url, v6["dir"] / f"agb_{y}.tif", BBOX)
    if sd_url:
        clip_remote(sd_url, v6["dir"] / f"agb_sd_{y}.tif", BBOX)

print("[ESA WorldCover]")
for y, url in WORLDCOVER.items():
    clip_remote(url, DATA / "worldcover" / f"worldcover_{y}.tif", BBOX)

print("[Copernicus DEM]")
clip_remote(DEM_URL, DATA / "dem" / "cop_dem30.tif", BBOX)
print("완료")
