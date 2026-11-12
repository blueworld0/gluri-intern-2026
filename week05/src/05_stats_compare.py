"""5주차 ①: 국가 통계(ha당 임목축적)로 위성(D) 결과 검증.

비교 두 가지
  (1) 흡수 속도: 통계 축적의 증가 속도(m³/ha/년) → 탄소(tC/ha/년)  vs  위성 D 2.07 (90% 1.24–2.91)
  (2) 수준: 같은 해(2015, 2020) 통계 축적 → 지상부 바이오매스(t/ha)  vs  위성 D 평균 AGB
환산 (IPCC, 국가 고유계수)
  지상부 바이오매스(t/ha) = 축적 × WD × BEF
  탄소(tC/ha)             = 축적 × WD × BEF × (1+R) × CF
  계수는 침엽/활엽을 금정산 임상 비율로 가중평균. '기본/낮게/높게' 세 경우로 계산
"""
import warnings
warnings.filterwarnings("ignore", message="Glyph")
import numpy as np
import pandas as pd
from config import *

OUT.mkdir(exist_ok=True)
df = pd.read_csv(STATS_CSV, encoding="utf-8-sig").dropna(subset=["ha당임목축적_m3ha"])
st = pd.read_csv(W4OUT / "stock_by_year_v6.csv")        # 위성 D 연도별 (금정산)
w = {k: v / sum(SHARE.values()) for k, v in SHARE.items()}

def factors(c):
    agb = sum(w[k] * c["WD"][k] * c["BEF"][k] for k in w)                 # m³ → t 바이오매스
    carbon = sum(w[k] * c["WD"][k] * c["BEF"][k] * (1 + c["R"][k]) for k in w) * CF_STAT  # m³ → tC
    return agb, carbon

# (1) 흡수 속도
rate_rows = []
for reg, g in df.groupby("지역"):
    g = g.sort_values("연도")
    periods = [(a, b) for a, b in zip(g.연도[:-1], g.연도[1:])] + ([(g.연도.min(), g.연도.max())] if len(g) > 2 else [])
    for y0, y1 in periods:
        v0 = g.loc[g.연도 == y0, "ha당임목축적_m3ha"].item(); v1 = g.loc[g.연도 == y1, "ha당임목축적_m3ha"].item()
        dv = (v1 - v0) / (y1 - y0)
        row = dict(지역=reg, 기간=f"{y0}→{y1}", 축적증가_m3_ha_년=dv)
        for nm, c in COEF.items():
            row[f"통계흡수_tC_ha_년_{nm}"] = dv * factors(c)[1]
        rate_rows.append(row)
rate = pd.DataFrame(rate_rows)
# 위성 D: 전체 추세 + 통계와 같은 기간(있으면)
d_same = {}
for y0, y1 in [(2015, 2020)]:
    if {y0, y1} <= set(st.연도):
        d_same[f"{y0}→{y1}"] = (st.loc[st.연도 == y1, "탄소밀도_tC_ha"].item() - st.loc[st.연도 == y0, "탄소밀도_tC_ha"].item()) / (y1 - y0)
rate["위성D_같은기간_tC_ha_년"] = rate.기간.map(d_same)
rate["위성D_2015_2022추세"] = D_SINK["대표"]
rate["기본값이_위성90%범위안"] = np.where((rate.통계흡수_tC_ha_년_기본 >= D_SINK["하한"]) & (rate.통계흡수_tC_ha_년_기본 <= D_SINK["상한"]), "예", "아니오")

# (2) 수준 (같은 해 AGB)
lvl_rows = []
for _, r in df.iterrows():
    if r.연도 in set(st.연도):
        row = dict(지역=r.지역, 연도=int(r.연도), 축적_m3_ha=r.ha당임목축적_m3ha,
                   위성D_AGB_t_ha=st.loc[st.연도 == r.연도, "평균AGB_t_ha"].item())
        for nm, c in COEF.items():
            row[f"통계AGB_t_ha_{nm}"] = r.ha당임목축적_m3ha * factors(c)[0]
        row["위성나누기통계_기본"] = row["위성D_AGB_t_ha"] / row["통계AGB_t_ha_기본"]
        lvl_rows.append(row)
lvl = pd.DataFrame(lvl_rows)

fac = pd.DataFrame([dict(경우=nm, AGB환산_t_per_m3=factors(c)[0], 탄소환산_tC_per_m3=factors(c)[1]) for nm, c in COEF.items()])
rate.round(2).to_csv(OUT / "stats_vs_satellite_rate.csv", index=False, encoding="utf-8-sig")
lvl.round(1).to_csv(OUT / "stats_vs_satellite_level.csv", index=False, encoding="utf-8-sig")
fac.round(3).to_csv(OUT / "stats_conversion_factors.csv", index=False, encoding="utf-8-sig")

# 그림
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in ["Malgun Gothic", "Noto Sans CJK KR", "NanumGothic", "Noto Sans KR", "WenQuanYi Zen Hei"]:
    if any(f == x.name for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams["axes.unicode_minus"] = False
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
# 왼쪽: 지상부 바이오매스 수준
a1.plot(st.연도, st.평균AGB_t_ha, "o-", color="#4C8C6A", label="위성 D (금정산, CCI v6)")
for reg, g in df.groupby("지역"):
    g = g.sort_values("연도")
    _, cb = factors(COEF["기본"])
    base = g.ha당임목축적_m3ha * factors(COEF["기본"])[0]
    lo = g.ha당임목축적_m3ha * factors(COEF["낮게"])[0]
    hi = g.ha당임목축적_m3ha * factors(COEF["높게"])[0]
    a1.plot(g.연도, base, "s--", color="#555", label=f"국가 통계 ({reg}, 기본 계수)")
    a1.fill_between(g.연도, lo, hi, color="#999", alpha=.25, label="계수 낮게~높게")
a1.set_ylabel("지상부 바이오매스 (t/ha)"); a1.set_title("① 수준 비교: ha당 지상부 바이오매스"); a1.legend(fontsize=8)
# 오른쪽: 흡수 속도
r = rate[rate.기간 == "2015→2020"] if (rate.기간 == "2015→2020").any() else rate.iloc[[-1]]
labels = ["위성 D\n2015–2022 추세", "위성 D\n2015→2020"] + [f"통계 {x}\n{p}" for x, p in zip(r.지역, r.기간)]
vals = [D_SINK["대표"], d_same.get("2015→2020", np.nan)] + list(r.통계흡수_tC_ha_년_기본)
a2.bar(labels, vals, color=["#4C8C6A", "#7FB89A"] + ["#8C8C8C"] * len(r))
a2.errorbar([0], [D_SINK["대표"]], yerr=[[D_SINK["대표"] - D_SINK["하한"]], [D_SINK["상한"] - D_SINK["대표"]]], color="k", capsize=6)
for i, (_, rr) in enumerate(r.iterrows()):
    a2.errorbar([2 + i], [rr.통계흡수_tC_ha_년_기본], yerr=[[rr.통계흡수_tC_ha_년_기본 - rr.통계흡수_tC_ha_년_낮게], [rr.통계흡수_tC_ha_년_높게 - rr.통계흡수_tC_ha_년_기본]], color="k", capsize=6)
a2.set_ylabel("흡수 속도 (tC/ha/년)"); a2.set_title("② 흡수 속도 비교 (막대: 위성 90% 범위 / 계수 낮게~높게)")
fig.tight_layout(); fig.savefig(OUT / "stats_vs_satellite.png", dpi=150)

pd.set_option("display.width", 220)
print("환산계수\n", fac.round(3).to_string(index=False))
print("\n흡수 속도\n", rate.round(2).to_string(index=False))
print("\n수준 (지상부 바이오매스)\n", lvl.round(1).to_string(index=False))
