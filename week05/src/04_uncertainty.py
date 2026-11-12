"""5주차 ②: 저장량·연간 흡수량의 불확도(범위)를 몬테카를로로 계산.

몬테카를로 = "모르는 값들을 그럴듯한 범위 안에서 무작위로 뽑아 계산을 수만 번 반복" →
결과 수만 개를 줄 세워 가운데 90%를 '범위'로 보고.

흔들어 보는 것 (오차의 원인 4가지)
  ① 지도 오차   : CCI가 칸마다 주는 표준편차(SD). 이웃 칸끼리 오차가 같이 움직이는 범위(블록)를 가정
  ② 연도 흔들림 : 8개 연도 점이 직선에서 벗어난 정도 → 추세 기울기(흡수량)의 표준오차
  ③ 뿌리비율 R  : 0.25 ~ 0.43 (대표 0.32)
  ④ 탄소비율 CF : 0.47 ~ 0.50
계산
  저장량(2022) = AGB총량 × (1 + 지도오차) × (1+R)/(1+0.32) × CF/0.5 … week04 값에 비율을 곱하는 방식
  흡수량       = 추세 기울기(+연도 흔들림) × (1 + 지도오차) × (1+R)/(1+0.32) × CF/0.5
                 ※ 지도가 전체적으로 x% 높게 잡으면 증가량도 x% 높게 잡힌다고 가정
"""
import json
import warnings
warnings.filterwarnings("ignore", message="Glyph")
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import reproject, Resampling
from shapely.geometry import shape

from config import *

OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(SEED)

# ── 0. week04 연도별 저장량 불러오기 (수종별 R, CF 0.5로 계산된 값) ──
st = pd.read_csv(W4OUT / "stock_by_year_v6.csv")
area_ha = st.산림면적_ha.iloc[0]

# ── 1. 지도 오차: 칸별 SD → 공원 전체 AGB의 상대오차 ──────────────
geom = shape(json.load(open(PROXY_BOUNDARY, encoding="utf-8"))["features"][0]["geometry"])
wc = rasterio.open(WORLDCOVER)
tr, shp = wc.transform, (wc.height, wc.width)
tree = (wc.read(1) == TREE_CLASS) & geometry_mask([geom], out_shape=shp, transform=tr, invert=True)
lat = tr.f + tr.e * (np.arange(shp[0]) + 0.5)
dx_m = abs(tr.a) * 111320 * np.cos(np.deg2rad(lat))          # 칸 가로 길이 (m)
dy_m = abs(tr.e) * 110574                                    # 칸 세로 길이 (m)
px_ha = np.broadcast_to((dx_m * dy_m / 1e4)[:, None], shp)
rows, cols = np.nonzero(tree)
y_m = rows * dy_m                                            # 격자 위의 대략적인 위치 (m)
x_m = cols * dx_m.mean()

def on_grid(path):
    out = np.zeros(shp, "float32")
    with rasterio.open(path) as s:
        reproject(s.read(1).astype("float32"), out, src_transform=s.transform, src_crs=s.crs,
                  dst_transform=tr, dst_crs=wc.crs, resampling=Resampling.nearest)
    return out

def rel_map_error(year, block_km):
    """공원 전체 AGB 총량의 상대오차(1σ).
    블록 안 칸들은 오차가 같은 방향(더함), 블록끼리는 서로 독립(제곱합의 제곱근)."""
    agb = on_grid(CCI_DIR / f"agb_{year}.tif")[tree] * px_ha[tree]
    sd = on_grid(CCI_DIR / f"agb_sd_{year}.tif")[tree] * px_ha[tree]
    if block_km is None:
        return sd.sum() / agb.sum()
    b = block_km * 1000
    key = (y_m // b).astype(np.int64) * 100000 + (x_m // b).astype(np.int64)
    block_sd = pd.Series(sd).groupby(key).sum().to_numpy()
    return np.sqrt((block_sd ** 2).sum()) / agb.sum()

last = max(YEARS)
blk = []
for bk in BLOCK_KM_TEST:
    e = [rel_map_error(y, bk) for y in YEARS]
    blk.append(dict(블록=("공원 전체(한 방향)" if bk is None else f"{bk:g} km"),
                    **{str(y): 100 * v for y, v in zip(YEARS, e)}, 평균_pct=100 * np.mean(e)))
blk = pd.DataFrame(blk)
blk.round(2).to_csv(OUT / "map_error_by_block.csv", index=False, encoding="utf-8-sig")
map_err = np.mean([rel_map_error(y, BLOCK_KM_MAIN) for y in YEARS])   # 대표값: 1 km, 8년 평균
print("지도 오차(공원 전체 AGB의 1σ, %)\n", blk[["블록", "평균_pct"]].round(1).to_string(index=False))

# ── 2. 연도 흔들림: 추세 기울기와 표준오차 ───────────────────────
def slope_se(df):
    x, y = df.연도.to_numpy(float), df.탄소_tC.to_numpy(float)
    b, a = np.polyfit(x, y, 1)
    res = y - (a + b * x)
    se = np.sqrt((res ** 2).sum() / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
    return b, se, len(x) - 2

# ── 3. 몬테카를로 ───────────────────────────────────────────────
def draws(n, use=("map", "year", "R", "CF")):
    z_map = rng.standard_normal(n) if "map" in use else np.zeros(n)
    t_year = rng.standard_t(df_, n) if "year" in use else np.zeros(n)
    R = rng.triangular(R_MIN, R_MODE, R_MAX, n) if "R" in use else np.full(n, R_BASE)
    CF = rng.uniform(CF_MIN, CF_MAX, n) if "CF" in use else np.full(n, CF_BASE)
    k = (1 + z_map * map_err) * (1 + R) / (1 + R_BASE) * CF / CF_BASE     # 공통 배율
    return k, t_year

def q(a):
    return np.percentile(a, [5, 50, 95])

stock_last = st.loc[st.연도 == last, "탄소_tC"].item()
rows_out = []
for label, sub in [("2015–2022 (대표)", st), ("2015–2021 (2022 제외)", st[st.연도 < 2022])]:
    b, se, df_ = slope_se(sub)
    k, t = draws(N_SIM)
    sink = (b + t * se) * k
    lo, md, hi = q(sink)
    rows_out.append(dict(항목=f"연간 흡수량 {label}", 단위="tC/년", 기본값=b, 중앙값=md, 하한_5pct=lo, 상한_95pct=hi,
                         **{"tCO2_기본값": b * CO2_PER_C, "tCO2_하한": lo * CO2_PER_C, "tCO2_상한": hi * CO2_PER_C}))
    if label.startswith("2015–2022"):
        # 원인별 기여: 하나만 흔들고 나머지는 고정했을 때 90% 범위의 폭
        contrib = []
        for src, nm in [("map", "① 지도 오차"), ("year", "② 연도 흔들림"), ("R", "③ 뿌리비율 R"), ("CF", "④ 탄소비율 CF")]:
            k1, t1 = draws(N_SIM, use=(src,))
            lo1, _, hi1 = q((b + t1 * se) * k1)
            contrib.append(dict(원인=nm, 흡수량_90pct범위_tC=f"{lo1:,.0f} ~ {hi1:,.0f}", 폭_tC=hi1 - lo1,
                                폭_기본값대비_pct=100 * (hi1 - lo1) / b))
        sink_main = sink
    # 저장량
k, _ = draws(N_SIM, use=("map", "R", "CF"))
stock = stock_last * k
lo, md, hi = q(stock)
rows_out.insert(0, dict(항목=f"저장량 {last}", 단위="tC", 기본값=stock_last, 중앙값=md, 하한_5pct=lo, 상한_95pct=hi,
                        tCO2_기본값=stock_last * CO2_PER_C, tCO2_하한=lo * CO2_PER_C, tCO2_상한=hi * CO2_PER_C))
summ = pd.DataFrame(rows_out)
summ["범위_기본값대비_pct"] = [f"{100*(l/b0-1):+.0f}% ~ {100*(h/b0-1):+.0f}%" for l, h, b0 in
                            zip(summ.하한_5pct, summ.상한_95pct, summ.기본값)]
summ.round(0).to_csv(OUT / "uncertainty_summary.csv", index=False, encoding="utf-8-sig")
contrib = pd.DataFrame(contrib).sort_values("폭_tC", ascending=False)
contrib.round(1).to_csv(OUT / "uncertainty_contribution.csv", index=False, encoding="utf-8-sig")
json.dump(dict(N_SIM=N_SIM, SEED=SEED, block_km=BLOCK_KM_MAIN, map_err_1sigma_pct=100 * map_err,
               R=[R_MIN, R_MODE, R_MAX], CF=[CF_MIN, CF_MAX], area_ha=area_ha),
          open(OUT / "uncertainty_settings.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# ── 4. 그림: 흡수량 분포 ─────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in ["Malgun Gothic", "Noto Sans CJK KR", "NanumGothic", "Noto Sans KR", "WenQuanYi Zen Hei"]:
    if any(f == x.name for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams["axes.unicode_minus"] = False
v = sink_main * CO2_PER_C / 1e4
lo, md, hi = q(v)
fig, ax = plt.subplots(figsize=(7, 3.8))
ax.hist(v, bins=80, color="#4C8C6A", alpha=.85)
for x_, ls in [(lo, "--"), (md, "-"), (hi, "--")]:
    ax.axvline(x_, color="k", ls=ls, lw=1)
ax.set_title(f"금정산 연간 흡수량 몬테카를로 {N_SIM:,}회 (2015–2022 추세)")
ax.set_xlabel("연간 흡수량 (만 tCO₂/년)")
ax.set_ylabel("횟수")
ax.text(0.98, 0.95, f"중앙값 {md:.1f}만\n90% 범위 {lo:.1f} ~ {hi:.1f}만", transform=ax.transAxes,
        ha="right", va="top")
fig.tight_layout()
fig.savefig(OUT / "sink_uncertainty_hist.png", dpi=150)

pd.set_option("display.width", 200)
print(f"\n대표 지도 오차(블록 {BLOCK_KM_MAIN:g} km, 8년 평균) = ±{100*map_err:.1f}% (1σ)")
print(summ[["항목", "단위", "기본값", "하한_5pct", "중앙값", "상한_95pct", "범위_기본값대비_pct"]].round(0).to_string(index=False))
print("\n원인별 기여 (하나만 흔들었을 때)\n", contrib.round(1).to_string(index=False))
