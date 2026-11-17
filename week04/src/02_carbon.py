"""3단계: 탄소 저장량과 연도별 변화 계산 (설악산 논문 방식, week04: CCI v4·v6 비교).

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
  8) week04: 같은 계산을 CCI v4(2017–2020)와 v6(2015–2022)로 각각 해서 비교
"""
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import reproject, Resampling
from shapely.geometry import shape

from config import (DATA, OUT, OFFICIAL_BOUNDARY, PROXY_BOUNDARY, CCI, CCI_MAIN, TREE_CLASS,
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



def stock_table(ver):
    """한 판(v4 또는 v6)의 연도별 저장량 표. 파일이 있는 연도만 계산."""
    info = CCI[ver]
    years = [y for y in info["years"] if (info["dir"] / f"agb_{y}.tif").exists()]
    rows, comp = [], []
    for y in years:
        agb = on_grid(info["dir"] / f"agb_{y}.tif")
        sd_path = info["dir"] / f"agb_sd_{y}.tif"
        sd = on_grid(sd_path) if sd_path.exists() else np.full(shp, np.nan, "float32")
        agb_t = (agb * px_ha)[tree].sum()               # 나무 픽셀만
        sd_t = (sd * px_ha)[tree].sum()                 # 오차가 모두 같은 방향이라 가정한 '최대' 불확도
        c = (agb * px_ha * (1 + r_px))[tree].sum() * CF          # 칸마다 (1+R) 을 곱해서 더함
        c_ipcc = agb_t * (1 + R_ROOT) * CF                          # 비교용: R=0.43
        comp.append(dict(연도=y, 탄소_IPCC043_tC=c_ipcc, 탄소_수종별R_tC=c, 차이_pct=100 * (c - c_ipcc) / c_ipcc))
        rows.append(dict(판=ver, 연도=y, 산림면적_ha=round(area_tree_ha, 1),
                         평균AGB_t_ha=agb_t / area_tree_ha, AGB_총량_t=agb_t,
                         탄소_tC=c, 탄소밀도_tC_ha=c / area_tree_ha, CO2환산_tCO2=c * CO2_PER_C,
                         AGB_불확도상한_pct=100 * sd_t / agb_t))
    return pd.DataFrame(rows), pd.DataFrame(comp)


def change_table(df):
    """축적차이법 변화량: 양 끝점, 직선추세, 해마다 차이의 들쭉날쭉한 정도."""
    ver = df.판.iloc[0]
    base, last = int(df.연도.min()), int(df.연도.max())
    c0, c1 = df.loc[df.연도 == base, "탄소_tC"].item(), df.loc[df.연도 == last, "탄소_tC"].item()
    slope = np.polyfit(df.연도, df.탄소_tC, 1)[0]
    yoy = df.sort_values("연도").탄소_tC.diff().dropna()          # 해마다 증감 (tC)
    out = [
        dict(판=ver, 구분=f"{base}→{last} 변화량", 탄소_tC=c1 - c0),
        dict(판=ver, 구분=f"{base}→{last} 연평균(양 끝점)", 탄소_tC=(c1 - c0) / (last - base)),
        dict(판=ver, 구분=f"{base}–{last} 연평균(직선추세)", 탄소_tC=slope),
        dict(판=ver, 구분="해마다 증감의 표준편차(들쭉날쭉한 정도)", 탄소_tC=yoy.std()),
        dict(판=ver, 구분="감소한 해의 수 / 전체", 탄소_tC=np.nan, 비고=f"{int((yoy < 0).sum())} / {len(yoy)}"),
    ]
    t = pd.DataFrame(out)
    t["CO2_tCO2"] = t.탄소_tC * CO2_PER_C
    t["tC_per_ha"] = t.탄소_tC / area_tree_ha
    return t


stocks, changes = {}, {}
for ver in CCI:
    df_v, comp_v = stock_table(ver)
    if df_v.empty:
        print(f"[{ver}] 파일 없음 → 건너뜀")
        continue
    stocks[ver], changes[ver] = df_v, change_table(df_v)
    df_v.round(2).to_csv(OUT / f"stock_by_year_{ver}.csv", index=False, encoding="utf-8-sig")
    changes[ver].round(2).to_csv(OUT / f"change_{ver}.csv", index=False, encoding="utf-8-sig")
    if R_MODE == "species":
        comp_v.round(2).to_csv(OUT / f"r_comparison_{ver}.csv", index=False, encoding="utf-8-sig")

# v4 와 v6 를 같은 연도끼리 나란히 (겹치는 2017–2020)
if len(stocks) == 2:
    a = stocks["v4"][["연도", "평균AGB_t_ha", "탄소_tC"]].rename(columns={"평균AGB_t_ha": "AGB_v4", "탄소_tC": "탄소_v4_tC"})
    b = stocks["v6"][["연도", "평균AGB_t_ha", "탄소_tC"]].rename(columns={"평균AGB_t_ha": "AGB_v6", "탄소_tC": "탄소_v6_tC"})
    cmp = a.merge(b, on="연도", how="outer").sort_values("연도")
    cmp["차이_v6-v4_pct"] = 100 * (cmp.탄소_v6_tC - cmp.탄소_v4_tC) / cmp.탄소_v4_tC
    for col in ("탄소_v4_tC", "탄소_v6_tC"):
        cmp[col.replace("탄소", "전년대비")] = cmp[col].diff()
    cmp.round(1).to_csv(OUT / "compare_v4_v6.csv", index=False, encoding="utf-8-sig")

# 대표 판(CCI_MAIN) 기준 민감도와 실행 정보
main = CCI_MAIN if CCI_MAIN in stocks else list(stocks)[-1]
df = stocks[main]
last = int(df.연도.max())
agb_last = df.loc[df.연도 == last, "AGB_총량_t"].item()
sens = pd.DataFrame([dict(계수=k, R=r, CF=CF, 탄소_tC=agb_last * (1 + r) * CF,
                          탄소밀도_tC_ha=agb_last * (1 + r) * CF / area_tree_ha) for k, r in R_SENSITIVITY.items()])
sens.round(2).to_csv(OUT / "sensitivity_R.csv", index=False, encoding="utf-8-sig")
json.dump(dict(boundary=boundary_kind, boundary_file=str(bpath.name), park_area_km2=area_park_ha / 100,
               tree_area_km2=area_tree_ha / 100, R_mode=R_MODE,
               R=float(r_px[tree].mean()) if R_MODE == "species" else R_ROOT, CF=CF,
               cci_main=main, cci_years={v: [int(x) for x in stocks[v].연도] for v in stocks}),
          open(OUT / "run_info.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)

pd.set_option("display.width", 200)
for ver in stocks:
    print(f"\n── CCI {ver} ──")
    print(stocks[ver].drop(columns="판").round(1).to_string(index=False))
    print(changes[ver].drop(columns="판").round(1).to_string(index=False))
if len(stocks) == 2:
    print("\n── v4 vs v6 ──"); print(cmp.round(1).to_string(index=False))
