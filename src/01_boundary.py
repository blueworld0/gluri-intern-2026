"""2단계: 임시 공원 경계 만들기.

공식 경계(국립공원공단/환경부 고시 도면)를 아직 못 받았기 때문에,
'금정산~백양산 능선에 이어진 산림' 영역을 자동으로 그려서 임시로 씁니다.

방법 (쉽게 말하면):
  1) WorldCover 2021 에서 '나무' 픽셀만 고른다.
  2) 그중 표고(DEM)가 h 미터 이상인 곳만 남긴다 (도심 공원·가로수 제외).
  3) 작은 틈은 메우고, 고당봉·백양산 정상과 이어진 덩어리만 남긴다.
  4) h 를 바꿔가며 면적이 공식 면적(66.859 km²)에 가장 가까운 것을 고른다.

※ 이건 어디까지나 '대용(proxy)' 경계입니다. 공식 경계를 받으면 config.py 의
  OFFICIAL_BOUNDARY 에 경로를 넣으면 모든 계산이 공식 경계로 바뀝니다.
"""
import json
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.features import shapes
from scipy import ndimage
from shapely.geometry import shape, mapping, Point
from shapely.ops import unary_union, transform as shp_transform
from pyproj import Transformer

from config import DATA, PROXY_BOUNDARY, OFFICIAL_AREA_KM2, TREE_CLASS, CRS_M

wc = rasterio.open(DATA / "worldcover" / "worldcover_2021.tif")
tree = wc.read(1) == TREE_CLASS
tr = wc.transform

# DEM 을 WorldCover 격자(약 10 m)에 맞춰 다시 그리기
dem = np.zeros(tree.shape, "float32")
with rasterio.open(DATA / "dem" / "cop_dem30.tif") as d:
    reproject(d.read(1), dem, src_transform=d.transform, src_crs=d.crs,
              dst_transform=tr, dst_crs=wc.crs, resampling=Resampling.bilinear)

# 두 봉우리 위치: 북쪽의 최고점 = 고당봉, 남쪽(35.20°N 이남)의 최고점 = 백양산
rows_lat = tr.f + tr.e * (np.arange(dem.shape[0]) + 0.5)
north = rows_lat > 35.25
south = rows_lat < 35.20
def peak(rowmask):
    d = np.where(rowmask[:, None], dem, -1)
    r, c = np.unravel_index(np.argmax(d), d.shape)
    return r, c, float(dem[r, c])
peaks = {"고당봉(북쪽 최고점)": peak(north), "백양산(남쪽 최고점)": peak(south)}
for k, (r, c, z) in peaks.items():
    lon, lat = tr * (c + 0.5, r + 0.5)
    print(f"{k}: {lat:.4f}N {lon:.4f}E, DEM {z:.0f} m")

# 픽셀 하나의 실제 면적(m²) — 위도에 따라 조금씩 다름
lat_c = np.deg2rad(rows_lat)
px_area = (abs(tr.a) * 111320 * np.cos(lat_c)) * (abs(tr.e) * 110574)

def region(h):
    m = tree & (dem >= h)
    m = ndimage.binary_closing(m, iterations=3)          # 작은 틈 메우기
    lab, _ = ndimage.label(m)
    keep = {lab[r, c] for r, c, _ in peaks.values()} - {0}
    m = np.isin(lab, list(keep))
    m = ndimage.binary_fill_holes(m)                     # 안쪽 구멍 메우기
    m = ndimage.binary_opening(m, iterations=2)          # 가는 꼬리 자르기
    area = (m * px_area[:, None]).sum() / 1e6
    return m, area

best = None
for h in range(40, 260, 10):
    m, a = region(h)
    print(f"  표고 ≥ {h:3d} m → {a:6.1f} km²")
    if best is None or abs(a - OFFICIAL_AREA_KM2) < abs(best[2] - OFFICIAL_AREA_KM2):
        best = (h, m, a)
h, m, a = best
print(f"선택: 표고 ≥ {h} m, 면적 {a:.1f} km² (공식 {OFFICIAL_AREA_KM2} km²)")

# 래스터 → 폴리곤, 약간 단순화
polys = [shape(g) for g, v in shapes(m.astype("uint8"), mask=m, transform=tr) if v == 1]
geom = unary_union(polys)
to_m = Transformer.from_crs("EPSG:4326", CRS_M, always_xy=True).transform
to_ll = Transformer.from_crs(CRS_M, "EPSG:4326", always_xy=True).transform
gm = shp_transform(to_m, geom).buffer(30).buffer(-30).simplify(15)
area_km2 = gm.area / 1e6
geom = shp_transform(to_ll, gm)

PROXY_BOUNDARY.parent.mkdir(parents=True, exist_ok=True)
fc = {"type": "FeatureCollection", "features": [{"type": "Feature", "geometry": mapping(geom),
      "properties": {"name": "금정산 임시경계(proxy)", "method": f"WorldCover2021 나무 & DEM>={h}m, 고당봉·백양산 연결 산림",
                     "area_km2": round(area_km2, 2), "official_area_km2": OFFICIAL_AREA_KM2}}]}
PROXY_BOUNDARY.write_text(json.dumps(fc, ensure_ascii=False), encoding="utf-8")
print("저장:", PROXY_BOUNDARY, f"{area_km2:.2f} km²")
