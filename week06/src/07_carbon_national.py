"""6주차 2단계: 새 R 지도(국가 고유계수)로 위성 D 저장량·흡수량 다시 계산.

week04 의 02_carbon.py 와 같은 계산이고, 칸별 R 만 06 단계에서 만든 새 지도로 바꿈.
  탄소(tC) = Σ [ 나무 칸 AGB(t/ha) × 칸 면적(ha) × (1 + R칸) ] × 0.5
  연간 흡수량 = 2015–2022 탄소의 직선 추세 기울기
"""
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import reproject, Resampling
from shapely.geometry import shape

from config import *

summ = json.load(open(OUT / "national_coef_summary.json", encoding="utf-8"))
geom = shape(json.load(open(PROXY_BOUNDARY, encoding="utf-8"))["features"][0]["geometry"])
wc = rasterio.open(WORLDCOVER)
tr, shp = wc.transform, (wc.height, wc.width)
tree = (wc.read(1) == TREE_CLASS) & geometry_mask([geom], out_shape=shp, transform=tr, invert=True)
lat = tr.f + tr.e * (np.arange(shp[0]) + 0.5)
px_ha = np.broadcast_to(((abs(tr.a) * 111320 * np.cos(np.deg2rad(lat))) * (abs(tr.e) * 110574) / 1e4)[:, None], shp)
area_tree_ha = px_ha[tree].sum()

with rasterio.open(R_MAP_NEW) as s_:
    rmap = s_.read(1)
r_px = np.where(rmap > 0, rmap, summ["면적가중_R"]).astype("float32")   # 임상도상 비산림 칸 = 평균 R

def on_grid(path):
    out = np.zeros(shp, "float32")
    with rasterio.open(path) as s:
        reproject(s.read(1).astype("float32"), out, src_transform=s.transform, src_crs=s.crs,
                  dst_transform=tr, dst_crs=wc.crs, resampling=Resampling.nearest)
    return out

old = pd.read_csv(W4OUT / "stock_by_year_v6.csv")
rows = []
for y in YEARS:
    agb = on_grid(CCI_DIR / f"agb_{y}.tif")
    agb_t = (agb * px_ha)[tree].sum()
    c = (agb * px_ha * (1 + r_px))[tree].sum() * CF
    c_old = old.loc[old.연도 == y, "탄소_tC"].item()
    rows.append(dict(연도=y, 산림면적_ha=area_tree_ha, 평균AGB_t_ha=agb_t / area_tree_ha, AGB_총량_t=agb_t,
                     탄소_tC=c, 탄소밀도_tC_ha=c / area_tree_ha, CO2환산_tCO2=c * CO2_PER_C,
                     week04_탄소_tC=c_old, 차이_pct=100 * (c - c_old) / c_old))
st = pd.DataFrame(rows)

def trend(df):
    return np.polyfit(df.연도, df.탄소_tC, 1)[0]
ch = []
for lab, sub in [("2015–2022 직선추세 (대표)", st), ("2015–2021 직선추세 (2022 제외)", st[st.연도 < 2022])]:
    b = trend(sub); b_old = np.polyfit(sub.연도, sub.week04_탄소_tC, 1)[0]
    ch.append(dict(구분=lab, 흡수량_tC_년=b, 흡수량_tCO2_년=b * CO2_PER_C, tC_ha_년=b / area_tree_ha,
                   week04_tC_년=b_old, 차이_pct=100 * (b - b_old) / b_old))
c0, c1 = st.탄소_tC.iloc[0], st.탄소_tC.iloc[-1]
ch.append(dict(구분="2015→2022 양 끝점 평균", 흡수량_tC_년=(c1 - c0) / 7, 흡수량_tCO2_년=(c1 - c0) / 7 * CO2_PER_C,
               tC_ha_년=(c1 - c0) / 7 / area_tree_ha))
ch = pd.DataFrame(ch)

st.round(2).to_csv(OUT / "stock_by_year_national.csv", index=False, encoding="utf-8-sig")
ch.round(2).to_csv(OUT / "sink_national.csv", index=False, encoding="utf-8-sig")
json.dump(dict(나무칸평균R=float(r_px[tree].mean()), 산림면적_ha=float(area_tree_ha)),
          open(OUT / "run_info_week06.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
pd.set_option("display.width", 220)
print(f"나무 칸 평균 R = {r_px[tree].mean():.4f}")
print(st.round(1).to_string(index=False))
print(ch.round(2).to_string(index=False))
