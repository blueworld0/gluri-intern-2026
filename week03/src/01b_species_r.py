"""2-2단계 (week03 개선): 임상도로 '금정산 맞춤' 뿌리비율 R 지도 만들기.

왜 하나?
  IPCC 기본값 R = 0.43 은 전 세계 공통값이라 금정산 숲과 맞지 않을 수 있습니다.
  뿌리비율은 수종마다 다르므로(소나무 0.25 ~ 신갈나무 0.39), 금정산에 실제로
  어떤 나무가 얼마나 있는지(임상도)를 보고 칸마다 그 수종의 R 을 넣습니다.

순서
  1) 임상도(부산·양산, 2019년판)에서 경계 안 부분만 잘라낸다.
  2) 수종 코드(KOFTR_GROU)마다 R 을 붙인다 (config.py 의 SPECIES_R).
  3) 수종별 면적표 + 면적가중 평균 R 을 outputs 에 저장한다.
  4) WorldCover 와 같은 10 m 격자에 'R 지도'를 그려 data/imsang 에 저장한다.
"""
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

from config import DATA, OUT, OFFICIAL_BOUNDARY, PROXY_BOUNDARY, IMSANG, SPECIES_R, R_MAP

OUT.mkdir(exist_ok=True)
bpath = OFFICIAL_BOUNDARY or PROXY_BOUNDARY
b = gpd.read_file(bpath).to_crs(5179)               # 임상도 좌표계(EPSG:5179, 미터)

# 1) 경계 안만 잘라 읽기 (부산 + 양산)
parts = []
for src in IMSANG:
    g = gpd.read_file(src, encoding="cp949", bbox=tuple(b.total_bounds))
    parts.append(gpd.clip(g[["FRTP_CD", "KOFTR_GROU", "geometry"]], b))
g = pd.concat(parts, ignore_index=True)
g["ha"] = g.area / 1e4

# 2) 수종별 R 붙이기 (표에 없는 코드 = 비산림·죽림 → R 없음)
g["수종"] = g.KOFTR_GROU.map(lambda c: SPECIES_R.get(c, ("(비산림/죽림 등)", np.nan))[0])
g["R"] = g.KOFTR_GROU.map(lambda c: SPECIES_R.get(c, (None, np.nan))[1])

# 3) 수종별 면적표 + 가중평균
t = g.groupby(["KOFTR_GROU", "수종"], dropna=False).agg(면적_ha=("ha", "sum"), R=("R", "first")).reset_index()
t["비율_pct"] = 100 * t.면적_ha / t.면적_ha.sum()
t = t.sort_values("면적_ha", ascending=False)
forest = t.dropna(subset=["R"])
r_mean = float((forest.R * forest.면적_ha).sum() / forest.면적_ha.sum())
t.round(3).to_csv(OUT / "species_R_table.csv", index=False, encoding="utf-8-sig")
frtp = g.groupby("FRTP_CD").ha.sum()
frtp_pct = {({"0": "비산림", "1": "침엽수림", "2": "활엽수림", "3": "혼효림", "4": "죽림"}.get(k, k)): round(float(100 * v / frtp.sum()), 1)
            for k, v in frtp.items()}
json.dump(dict(임상도="2019년판(부산·양산)", 경계=bpath.name, 산림_면적가중평균R=round(r_mean, 4),
               임상별_면적비율_pct=frtp_pct), open(OUT / "species_R_summary.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print(t.head(10).round(2).to_string(index=False))
print("임상별 비율(%):", frtp_pct)
print(f"금정산 산림 면적가중 평균 R = {r_mean:.3f}  (IPCC 기본 0.43)")

# 4) R 지도 (WorldCover 격자와 똑같이)
with rasterio.open(DATA / "worldcover" / "worldcover_2020.tif") as wc:
    prof, tr, shp, crs = wc.profile.copy(), wc.transform, (wc.height, wc.width), wc.crs
gf = g.dropna(subset=["R"]).to_crs(crs)
rmap = rasterize(zip(gf.geometry, gf.R), out_shape=shp, transform=tr, fill=0, dtype="float32")
prof.update(dtype="float32", nodata=0, compress="deflate")
with rasterio.open(R_MAP, "w", **prof) as dst:
    dst.write(rmap, 1)
print("저장:", R_MAP)
