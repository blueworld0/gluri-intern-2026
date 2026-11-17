"""week04 · ESA CCI Biomass v6.0 (2015–2022)를 CEDA에서 금정산 주변만 잘라 받기.

실행:  week04 폴더에서  python download_cci_v6.py
결과:  data/cci_biomass_v6/agb_{연도}.tif, agb_sd_{연도}.tif  (100 m, t/ha)

- CEDA는 로그인 없이 받을 수 있고, 라이선스상 이용 목적 제한이 없음(Permitted use: Any).
- 타일 전체(수백 MB)를 받지 않고 필요한 부분만 인터넷으로 읽어서 자름.
- 원격 읽기가 안 되면 타일 전체를 임시로 받아서 자르고 지움(느림).
"""
import os
import re
import sys
import urllib.request
from pathlib import Path

import rasterio
from rasterio.windows import from_bounds

BBOX = (128.98, 35.12, 129.14, 35.36)          # week03 config.py 와 같은 범위
YEARS = list(range(2015, 2023))
BASE = "https://dap.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v6.0/geotiff"
TILE = "N40E120"                                # 30–40°N, 120–130°E 타일 (금정산 포함)
OUT = Path(__file__).resolve().parent / "data" / "cci_biomass_v6"

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif,.TIF")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "5")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")


def list_tifs(url):
    """CEDA 폴더 페이지에서 금정산 타일의 tif 파일 주소 찾기."""
    try:
        html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "ignore")
    except Exception as e:
        print("   폴더 열기 실패:", url, "→", e)
        return []
    names = sorted(set(re.findall(rf'({TILE}[^"<>/]*?\.tif)', html)))
    return [url.rstrip("/") + "/" + n for n in names]


def clip(url, out_path):
    if out_path.exists():
        print("   이미 있음:", out_path.name)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        src = rasterio.open("/vsicurl/" + url)
    except Exception as e:
        print("   원격 읽기 실패 → 타일 전체를 임시로 받습니다(느림):", e)
        tmp = out_path.with_suffix(".full.tif")
        urllib.request.urlretrieve(url, tmp)
        src = rasterio.open(tmp)
    with src:
        win = from_bounds(*BBOX, transform=src.transform).round_offsets().round_lengths()
        arr = src.read(1, window=win)
        prof = src.profile.copy()
        prof.update(driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                    transform=src.window_transform(win), compress="deflate", count=1)
        for k in ("blockxsize", "blockysize", "tiled", "photometric", "interleave"):
            prof.pop(k, None)
        with rasterio.open(out_path, "w", **prof) as dst:
            dst.write(arr, 1)
    tmp = out_path.with_suffix(".full.tif")
    if tmp.exists():
        tmp.unlink()
    print("   저장:", out_path.name, arr.shape, "최대값", float(arr.max()))


found_any = False
for y in YEARS:
    print(f"[{y}]")
    tifs = list_tifs(f"{BASE}/{y}/") or [u for u in list_tifs(f"{BASE}/") if f"-{y}-" in u]
    agb = [u for u in tifs if "AGB-MERGED" in u or ("AGB" in u and "_SD" not in u and "SD-" not in u)]
    sd = [u for u in tifs if "AGB_SD" in u or "SD-MERGED" in u]
    if not agb:
        print("   이 연도 파일을 못 찾음. 찾은 목록:", tifs or "없음")
        continue
    found_any = True
    print("   AGB:", agb[0].split("/")[-1])
    clip(agb[0], OUT / f"agb_{y}.tif")
    if sd:
        print("   SD :", sd[0].split("/")[-1])
        clip(sd[0], OUT / f"agb_sd_{y}.tif")

if not found_any:
    print("\n파일을 하나도 못 찾았어요. 아래 주소를 브라우저로 열어서 폴더 구조를 캡처해 주세요:")
    print(" ", BASE + "/")
    sys.exit(1)
print("\n끝! 저장 위치:", OUT)
