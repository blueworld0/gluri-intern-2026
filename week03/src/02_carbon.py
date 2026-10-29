"""3단계: 탄소 저장량과 연도별 변화 계산 (설악산 논문 방식).

설악산 논문(전청옥 외 2024, 국립공원연구지)과 같은 순서입니다.
  1) 공원 경계 안에서
  2) WorldCover 2020 '나무' 픽셀만 골라 (10 m)
  3) 그 자리의 ESA CCI 바이오매스(AGB, t/ha, 100 m)를 읽어
  4) AGB 총량(t) = Σ (AGB t/ha × 픽셀면적 ha)
  5) 탄소(tC) = AGB × (1 + R) × CF       R=뿌리비율, CF=탄소비율 0.5
     - R_MODE="ipcc"    : R=0.43 하나 (설악산 논문 방식)
     - R_MODE="species" : 임상도 수종별 R 을 칸마다 적용 (01b_species_r.py, week03 개선)
  6) CO2(tCO2) = 탄소 × 44/12
  7) 연도별 차이로 변화량, 연평균 흡수량 계산 (축적차이법)
"""
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import reproject, Resampling
from shapely.geometry import shape

from config import (DATA, OUT, OFFICIAL_BOUNDARY, PROXY_BOUNDARY, CCI_YEARS, TREE_CLASS,
                    R_ROOT, CF, CO2_PER_C, R_SENSITIVITY, R_MODE, R_MAP)

OUT.mkdir(exist_ok=True)
bpath = OFFICIAL_BOUNDARY or PROXY_BOUNDARY
boundary_kind = "공식" if OFFICIAL_BOUNDARY else "임시(proxy)"
geom = shape(json.loads(open(bpath, encoding="utf-8").read())["features"][0]["geometry"])
print(f"경계: {boundary_kind} — {bpath}")

# 기준 격자 = WorldCover 2020 (약 10 m, 경위도)
wc = rasterio.open(DATA / "worldcover" / "worldcover_2020.tif")
tr, shp = wc.transform, (wc.height, wc.width)
inside = geometry_mask([geom], out_shape=shp, transform=tr, invert=True)   # 경계 안 = True
tree = (wc.read(1) == TREE_CLASS) & inside

# 픽셀 면적(ha): 경위도 격자라 위도마다 조금씩 다름
lat = tr.f + tr.e * (np.arange(shp[0]) + 0.5)
px_ha = (abs(tr.a) * 111320 * np.cos(np.deg2rad(lat))) * (abs(tr.e) * 110574) / 1e4
px_ha = np.broadcast_to(px_ha[:, None], shp)

area_park_ha = px_ha[inside].sum()
area_tree_ha = px_ha[tree].sum()
print(f"경계 면적 {area_park_ha/100:.2f} km², 그중 나무 {area_tree_ha/100:.2f} km² ({area_tree_ha/area_park_ha:.1%})")

def on_grid(path):
    """100 m CCI 값을 10 m 격자에 그대로(최근린) 옮기기."""
    out = np.zeros(shp, "float32")
    with rasterio.open(path) as s:
        reproject(s.read(1).astype("float32"), out, src_transform=s.transform, src_crs=s.crs,
                  dst_transform=tr, dst_crs=wc.crs, resampling=Resampling.nearest)
    return out

# 칸별 뿌리비율 R 지도
if R_MODE == "species":
    with rasterio.open(R_MAP) as s_:
        rmap = s_.read(1)
    r_fill = json.load(open(OUT / "species_R_summary.json", encoding="utf-8"))["산림_면적가중평균R"]
    r_px = np.where(rmap > 0, rmap, r_fill).astype("float32")   # 임상도상 비산림·죽림인 나무 칸은 평균값
    print(f"R 방식: 수종별 (임상도 2019) — 나무 칸 평균 R = {r_px[tree].mean():.3f}")
else:
    r_px = np.full(shp, R_ROOT, "float32")
    print(f"R 방식: IPCC 기본값 {R_ROOT}")

rows, comp = [], []
for y in CCI_YEARS:
    agb = on_grid(DATA / "cci_biomass_v4" / f"agb_{y}.tif")
    sd = on_grid(DATA / "cci_biomass_v4" / f"agb_sd_{y}.tif")
    agb_t = (agb * px_ha)[tree].sum()               # 나무 픽셀만
    agb_all_t = (agb * px_ha)[inside].sum()         # 참고: 경계 안 전체 픽셀
    sd_t = (sd * px_ha)[tree].sum()                 # 오차가 모두 같은 방향이라 가정한 '최대' 불확도
    c = (agb * px_ha * (1 + r_px))[tree].sum() * CF          # 칸마다 (1+R) 을 곱해서 더함
    c_ipcc = agb_t * (1 + R_ROOT) * CF                          # 비교용: R=0.43
    comp.append(dict(연도=y, 탄소_IPCC043_tC=c_ipcc, 탄소_수종별R_tC=c, 차이_pct=100 * (c - c_ipcc) / c_ipcc))
    rows.append(dict(연도=y, 산림면적_ha=round(area_tree_ha, 1),
                     평균AGB_t_ha=agb_t / area_tree_ha, AGB_총량_t=agb_t,
                     탄소_tC=c, 탄소밀도_tC_ha=c / area_tree_ha, CO2환산_tCO2=c * CO2_PER_C,
                     AGB_불확도상한_pct=100 * sd_t / agb_t,
                     참고_경계전체AGB_t=agb_all_t))
df = pd.DataFrame(rows)

# 변화량 (축적차이법): 기준연도 대비
base, last = CCI_YEARS[0], CCI_YEARS[-1]
c0, c1 = df.loc[df.연도 == base, "탄소_tC"].item(), df.loc[df.연도 == last, "탄소_tC"].item()
slope = np.polyfit(df.연도, df.탄소_tC, 1)[0]       # 4개 연도 직선 추세 (tC/년)
change = pd.DataFrame([
    dict(구분=f"{base}→{last} 변화량", 탄소_tC=c1 - c0, CO2_tCO2=(c1 - c0) * CO2_PER_C),
    dict(구분=f"{base}→{last} 연평균(양 끝점)", 탄소_tC=(c1 - c0) / (last - base), CO2_tCO2=(c1 - c0) / (last - base) * CO2_PER_C),
    dict(구분=f"{base}–{last} 연평균(직선추세)", 탄소_tC=slope, CO2_tCO2=slope * CO2_PER_C),
    dict(구분="2018→2020 연평균(참고)", 탄소_tC=(c1 - df.loc[df.연도 == 2018, "탄소_tC"].item()) / 2,
         CO2_tCO2=(c1 - df.loc[df.연도 == 2018, "탄소_tC"].item()) / 2 * CO2_PER_C),
])

# 민감도: 뿌리비율 R 을 한국 수종별 값으로 바꾸면?
agb_last = df.loc[df.연도 == last, "AGB_총량_t"].item()
sens = pd.DataFrame([dict(계수=k, R=r, CF=CF, 탄소_tC=agb_last * (1 + r) * CF,
                          탄소밀도_tC_ha=agb_last * (1 + r) * CF / area_tree_ha) for k, r in R_SENSITIVITY.items()])

df.round(2).to_csv(OUT / "stock_by_year.csv", index=False, encoding="utf-8-sig")
if R_MODE == "species":
    pd.DataFrame(comp).round(2).to_csv(OUT / "r_comparison.csv", index=False, encoding="utf-8-sig")
change.round(1).to_csv(OUT / "change.csv", index=False, encoding="utf-8-sig")
sens.round(2).to_csv(OUT / "sensitivity_R.csv", index=False, encoding="utf-8-sig")
json.dump(dict(boundary=boundary_kind, boundary_file=str(bpath.name), park_area_km2=area_park_ha / 100,
               tree_area_km2=area_tree_ha / 100, R_mode=R_MODE,
               R=float(r_px[tree].mean()) if R_MODE == "species" else R_ROOT, CF=CF), open(OUT / "run_info.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)

pd.set_option("display.width", 200)
print(df.round(1).to_string(index=False)); print(); print(change.round(0).to_string(index=False))
print(); print(sens.round(2).to_string(index=False))
if R_MODE == "species":
    print(); print(pd.DataFrame(comp).round(1).to_string(index=False))
