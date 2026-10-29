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

from config import DATA, OUT, OFFICIAL_BOUNDARY, PROXY_BOUNDARY, R_ROOT, CF, R_MODE

# 한글 글꼴 (Windows: 맑은 고딕, 리눅스: Noto Sans CJK / 나눔)
for f in ["Malgun Gothic", "Noto Sans CJK KR", "NanumGothic", "Noto Sans KR"]:
    if any(f == x.name for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams["axes.unicode_minus"] = False

bpath = OFFICIAL_BOUNDARY or PROXY_BOUNDARY
b = gpd.read_file(bpath)
label = "공식 경계" if OFFICIAL_BOUNDARY else "임시 경계(proxy)"
df = pd.read_csv(OUT / "stock_by_year.csv")
last = int(df.연도.max())

# (1) 탄소밀도 지도: AGB × (1+R) × CF  (tC/ha)
r_used = json.load(open(OUT / "run_info.json", encoding="utf-8"))["R"]   # 수종별이면 평균 R
r_label = f"수종별 R(평균 {r_used:.2f})" if R_MODE == "species" else f"R {R_ROOT}"
with rasterio.open(DATA / "cci_biomass_v4" / f"agb_{last}.tif") as s:
    c = s.read(1).astype(float) * (1 + r_used) * CF
    ext = [s.bounds.left, s.bounds.right, s.bounds.bottom, s.bounds.top]
fig, ax = plt.subplots(figsize=(6, 8.5))
im = ax.imshow(np.ma.masked_equal(c, 0), extent=ext, cmap="Greens", vmin=0, vmax=140)
b.boundary.plot(ax=ax, color="#222", lw=1.4)
ax.set_title(f"금정산 산림 탄소밀도 {last} (ESA CCI Biomass v4, 100 m)\n검은 선 = {label}", fontsize=11)
ax.set_xlabel("경도"); ax.set_ylabel("위도")
cb = fig.colorbar(im, ax=ax, shrink=0.6); cb.set_label(f"탄소밀도 (tC/ha) = AGB×(1+R)×0.5, {r_label}")
fig.tight_layout(); fig.savefig(OUT / f"map_carbon_density_{last}.png", dpi=150); plt.close(fig)

# (2) 연도별 탄소 저장량 막대그래프
fig, ax = plt.subplots(figsize=(6, 3.8))
x = df.연도.astype(str)
y = df.탄소_tC / 1e4
ax.bar(x, y, color="#2f7d4f", width=0.55)
for xi, yi in zip(x, y):
    ax.text(xi, yi + 0.6, f"{yi:.1f}", ha="center", fontsize=9, color="#333")
ax.set_ylabel("탄소 저장량 (만 tC)")
ax.set_title(f"연도별 탄소 저장량 (지상부+지하부, {label}, {r_label})", fontsize=10)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.grid(axis="y", color="#ddd", lw=0.6); ax.set_axisbelow(True)
ax.set_ylim(0, y.max() * 1.15)
fig.tight_layout(); fig.savefig(OUT / "bar_carbon_by_year.png", dpi=150); plt.close(fig)
print("저장:", OUT / f"map_carbon_density_{last}.png", OUT / "bar_carbon_by_year.png")
