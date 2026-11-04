"""4단계: 보고용 그림 만들기 (지도 + 연도별 막대그래프)."""
import json
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

from config import DATA, OUT, OFFICIAL_BOUNDARY, PROXY_BOUNDARY, R_ROOT, CF, R_MODE, CCI

# 한글 글꼴 (Windows: 맑은 고딕, 리눅스: Noto Sans CJK / 나눔)
for f in ["Malgun Gothic", "Noto Sans CJK KR", "NanumGothic", "Noto Sans KR"]:
    if any(f == x.name for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams["axes.unicode_minus"] = False

bpath = OFFICIAL_BOUNDARY or PROXY_BOUNDARY
b = gpd.read_file(bpath)
label = "공식 경계" if OFFICIAL_BOUNDARY else "임시 경계(proxy)"
info = json.load(open(OUT / "run_info.json", encoding="utf-8"))
main = info["cci_main"]                                   # 대표 판 (보통 v6)
df = pd.read_csv(OUT / f"stock_by_year_{main}.csv")
last = int(df.연도.max())

# (1) 탄소밀도 지도: AGB × (1+R) × CF  (tC/ha)
r_used = info["R"]   # 수종별이면 평균 R
r_label = f"수종별 R(평균 {r_used:.2f})" if R_MODE == "species" else f"R {R_ROOT}"
with rasterio.open(CCI[main]["dir"] / f"agb_{last}.tif") as s:
    c = s.read(1).astype(float) * (1 + r_used) * CF
    ext = [s.bounds.left, s.bounds.right, s.bounds.bottom, s.bounds.top]
fig, ax = plt.subplots(figsize=(6, 8.5))
im = ax.imshow(np.ma.masked_equal(c, 0), extent=ext, cmap="Greens", vmin=0, vmax=140)
b.boundary.plot(ax=ax, color="#222", lw=1.4)
ax.set_title(f"금정산 산림 탄소밀도 {last} (ESA CCI Biomass {main}, 100 m)\n검은 선 = {label}", fontsize=11)
ax.set_xlabel("경도"); ax.set_ylabel("위도")
cb = fig.colorbar(im, ax=ax, shrink=0.6); cb.set_label(f"탄소밀도 (tC/ha) = AGB×(1+R)×0.5, {r_label}")
fig.tight_layout(); fig.savefig(OUT / f"map_carbon_density_{last}.png", dpi=150); plt.close(fig)

# (2) 연도별 탄소 저장량: v4 와 v6 를 나란히 (같은 해끼리 비교)
tabs = {v: pd.read_csv(OUT / f"stock_by_year_{v}.csv") for v in info["cci_years"]}
years = sorted(set().union(*[set(t.연도) for t in tabs.values()]))
colors = {"v4": "#9db8a6", "v6": "#2f7d4f"}
w = 0.38 if len(tabs) == 2 else 0.6
fig, ax = plt.subplots(figsize=(7.5, 4))
for k, (v, t) in enumerate(tabs.items()):
    off = (k - (len(tabs) - 1) / 2) * w
    xs = [years.index(y) + off for y in t.연도]
    ys = t.탄소_tC / 1e4
    ax.bar(xs, ys, width=w * 0.95, color=colors.get(v, "#555"), label=f"CCI {v}")
    for xi, yi in zip(xs, ys):
        ax.text(xi, yi + 0.6, f"{yi:.1f}", ha="center", fontsize=8, color="#333")
ax.set_xticks(range(len(years))); ax.set_xticklabels([str(y) for y in years])
ax.set_ylabel("탄소 저장량 (만 tC)")
ax.set_title(f"연도별 탄소 저장량: CCI v4 vs v6 (지상부+지하부, {label}, {r_label})", fontsize=10)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.grid(axis="y", color="#ddd", lw=0.6); ax.set_axisbelow(True)
ymax = max(t.탄소_tC.max() for t in tabs.values()) / 1e4
ax.set_ylim(0, ymax * 1.18); ax.legend(frameon=False, loc="upper left", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "bar_carbon_v4_v6.png", dpi=150); plt.close(fig)
print("저장:", OUT / f"map_carbon_density_{last}.png", OUT / "bar_carbon_v4_v6.png")
