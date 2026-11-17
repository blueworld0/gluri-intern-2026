"""6주차 4단계: 국가 통계(ha당 임목축적) vs 위성 D, 국가 고유계수로 다시 비교.

5주차 05_stats_compare.py 와 같은 비교이고, 환산계수만 06 단계 결과로 바꿈.
  지상부 바이오매스(t/ha) = 축적(m³/ha) × [WD × BEF]의 금정산 가중평균
  탄소(tC/ha)             = 축적 × [WD × BEF × (1+R)]의 가중평균 × 0.5
  경우: 기본(임상도 수종별) / 낮게·높게(침엽·활엽 각각 표의 최저·최고 수종) / 5주차 임시값(참고)
통계 입력: week06/data_stats/forest_stock.csv (전국 2010·2015·2020 + 부산광역시·금정구 2015·2020 + 북구 2020). 다른 지역·연도는 같은 형식으로 줄 추가.
"""
import warnings
warnings.filterwarnings("ignore", message="Glyph")
import numpy as np
import pandas as pd
from config import *

fac = pd.read_csv(OUT / "conversion_factors_national.csv")
CASES = {"기본": fac.iloc[0], "낮게": fac.iloc[1], "높게": fac.iloc[2], "5주차임시": fac.iloc[3]}
df = pd.read_csv(STATS_CSV, encoding="utf-8-sig").dropna(subset=["ha당임목축적_m3ha"])
st = pd.read_csv(OUT / "stock_by_year_national.csv")
unc = pd.read_csv(OUT / "uncertainty_summary.csv")
area = st.산림면적_ha.iloc[0]
u = unc[unc.항목.str.contains("2015–2022")].iloc[0]
D = dict(대표=u.기본값 / area, 하한=u.하한_5pct / area, 상한=u.상한_95pct / area)

# (1) 흡수 속도
rate_rows = []
for reg, g in df.groupby("지역"):
    g = g.sort_values("연도")
    yrs = list(g.연도)
    periods = list(zip(yrs[:-1], yrs[1:])) + ([(yrs[0], yrs[-1])] if len(yrs) > 2 else [])
    for y0, y1 in periods:
        v0 = g.loc[g.연도 == y0, "ha당임목축적_m3ha"].item(); v1 = g.loc[g.연도 == y1, "ha당임목축적_m3ha"].item()
        dv = (v1 - v0) / (y1 - y0)
        row = dict(지역=reg, 기간=f"{y0}→{y1}", 축적증가_m3_ha_년=dv)
        for nm, c in CASES.items():
            row[f"통계흡수_tC_ha_년_{nm}"] = dv * c.탄소환산_tC_per_m3
        rate_rows.append(row)
rate = pd.DataFrame(rate_rows)
d_same = {}
for y0, y1 in [(2015, 2020)]:
    d_same[f"{y0}→{y1}"] = (st.loc[st.연도 == y1, "탄소밀도_tC_ha"].item() - st.loc[st.연도 == y0, "탄소밀도_tC_ha"].item()) / (y1 - y0)
rate["위성D_같은기간_tC_ha_년"] = rate.기간.map(d_same)
rate["위성D_2015_2022추세"] = D["대표"]
rate["차이_기본vs위성추세_pct"] = 100 * (rate.통계흡수_tC_ha_년_기본 - D["대표"]) / D["대표"]
rate["기본값이_위성90%범위안"] = np.where(rate.통계흡수_tC_ha_년_기본.between(D["하한"], D["상한"]), "예", "아니오")

# (2) 수준 (같은 해 AGB)
lvl_rows = []
for _, r in df.iterrows():
    if r.연도 in set(st.연도):
        row = dict(지역=r.지역, 연도=int(r.연도), 축적_m3_ha=r.ha당임목축적_m3ha,
                   위성D_AGB_t_ha=st.loc[st.연도 == r.연도, "평균AGB_t_ha"].item())
        for nm, c in CASES.items():
            row[f"통계AGB_t_ha_{nm}"] = r.ha당임목축적_m3ha * c.AGB환산_t_per_m3
        row["차이_위성vs통계기본_pct"] = 100 * (row["위성D_AGB_t_ha"] / row["통계AGB_t_ha_기본"] - 1)
        lvl_rows.append(row)
lvl = pd.DataFrame(lvl_rows)
rate.round(3).to_csv(OUT / "stats_vs_satellite_rate.csv", index=False, encoding="utf-8-sig")
lvl.round(2).to_csv(OUT / "stats_vs_satellite_level.csv", index=False, encoding="utf-8-sig")

pd.set_option("display.width", 240)
print(f"위성 D 추세 {D['대표']:.2f} tC/ha/년 (90% {D['하한']:.2f}–{D['상한']:.2f})")
print(rate.round(2).to_string(index=False))
print(lvl.round(1).to_string(index=False))

# ── (3) 임상별 보정: 침엽·활엽·혼효의 ha당 축적을 따로 환산 ──────────────
# 왜? 구 통계 평균(177 m³/ha)에는 '부피는 크지만 가벼운' 침엽수가 많이 섞여 있음(금정구 산림의 46%).
#     금정산 경계 안은 활엽수가 더 많음(임상도 2019). 같은 평균 부피라도 임상 구성에 따라 무게가 달라짐.
# 방법: KOSIS 「행정구역별 임상별 면적 및 축적」(2020)으로 임상별 ha당 축적을 구하고,
#       ① 금정산 임상 비율(임상도)로 가중, ② 임상별 환산계수(금정산 수종 가중)를 따로 곱함
TYPE_CSV = ROOT / "data_stats" / "busan_forest_type_2020_kosis.csv"
if TYPE_CSV.exists():
    import json
    raw = pd.read_csv(TYPE_CSV, encoding="utf-8-sig", header=None, skiprows=3).replace("-", 0)
    raw.columns = ["시도", "구", "면적계", "면적침엽", "면적활엽", "면적혼효", "면적죽림", "면적무립목",
                   "축적계", "축적침엽", "축적활엽", "축적혼효"]
    num = raw.columns[2:]
    raw[num] = raw[num].apply(pd.to_numeric, errors="coerce").fillna(0)
    sp = pd.read_csv(OUT / "species_coef_table.csv").dropna(subset=["WD"])
    sp["임상"] = np.select([sp.침엽몫 == 1, sp.침엽몫 == 0], ["침엽", "활엽"], "혼효")
    tf = {k: (g.면적_ha * g.WDxBEF).sum() / g.면적_ha.sum() for k, g in sp.groupby("임상")}          # m³ → t
    tc = {k: (g.면적_ha * g.WDxBEFx1R).sum() / g.면적_ha.sum() * CF for k, g in sp.groupby("임상")}   # m³ → tC
    share = {k: g.면적_ha.sum() for k, g in sp.groupby("임상")}
    share = {k: v / sum(share.values()) for k, v in share.items()}                                   # 금정산 임상 비율
    agb_sat = st.loc[st.연도 == 2020, "평균AGB_t_ha"].item()
    rows3 = []
    for _, r in raw.iterrows():
        vol = {k: float(r[f"축적{k}"]) / float(r[f"면적{k}"]) for k in ("침엽", "활엽", "혼효")}
        own = {k: float(r[f"면적{k}"]) for k in ("침엽", "활엽", "혼효")}
        own = {k: v / sum(own.values()) for k, v in own.items()}
        agb_mix = sum(share[k] * vol[k] * tf[k] for k in vol)            # 금정산 임상 비율 + 임상별 계수
        rows3.append(dict(지역=f"{r.시도} {r.구}", 침엽_m3_ha=vol["침엽"], 활엽_m3_ha=vol["활엽"], 혼효_m3_ha=vol["혼효"],
                          통계_침엽비율_pct=100 * own["침엽"], 금정산_침엽비율_pct=100 * share["침엽"],
                          통계AGB_평균축적방식_t_ha=float(r.축적계) / float(r.면적계) * CASES["기본"].AGB환산_t_per_m3,
                          통계AGB_임상보정_t_ha=agb_mix, 위성D_AGB_2020_t_ha=agb_sat,
                          차이_위성vs임상보정_pct=100 * (agb_sat / agb_mix - 1)))
    t3 = pd.DataFrame(rows3)
    t3.round(2).to_csv(OUT / "stats_vs_satellite_by_forest_type.csv", index=False, encoding="utf-8-sig")
    json.dump(dict(임상별_AGB환산_t_per_m3=tf, 임상별_탄소환산_tC_per_m3=tc, 금정산_임상비율=share),
              open(OUT / "forest_type_factors.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n(3) 임상별 보정\n", t3.round(1).to_string(index=False))
    print({k: round(v, 3) for k, v in tf.items()}, {k: round(v, 3) for k, v in share.items()})

    # (3-2) 흡수 속도도 임상 보정: 2015 표가 있으면 2015→2020 탄소밀도 차이
    TYPE15 = ROOT / "data_stats" / "busan_forest_type_2015_kosis.csv"
    if TYPE15.exists():
        def read_type(path):
            d = pd.read_csv(path, encoding="utf-8-sig", header=None, skiprows=3).replace("-", 0)
            d.columns = raw.columns[:len(d.columns)]
            d[d.columns[2:]] = d[d.columns[2:]].apply(pd.to_numeric, errors="coerce").fillna(0)
            d["구"] = d.구.replace({"소계": "계"})
            return d.set_index("구")
        a15, a20 = read_type(TYPE15), raw.assign(구=raw.구.replace({"소계": "계"})).set_index("구")
        rows4 = []
        for gu in [g for g in a15.index if g in a20.index]:
            def cden(r):   # 금정산 임상 비율 × 임상별 ha당 축적 × 임상별 탄소계수 → tC/ha
                return sum(share[k] * float(r[f"축적{k}"]) / float(r[f"면적{k}"]) * tc[k] for k in ("침엽", "활엽", "혼효"))
            c15, c20 = cden(a15.loc[gu]), cden(a20.loc[gu])
            avg = (float(a20.loc[gu].축적계) / float(a20.loc[gu].면적계) - float(a15.loc[gu].축적계) / float(a15.loc[gu].면적계)) / 5
            rows4.append(dict(지역=f"부산광역시 {gu}", 탄소밀도_2015_임상보정=c15, 탄소밀도_2020_임상보정=c20,
                              흡수_임상보정_tC_ha_년=(c20 - c15) / 5,
                              흡수_평균축적방식_tC_ha_년=avg * CASES["기본"].탄소환산_tC_per_m3,
                              위성D_2015_2022추세=D["대표"], 위성D_2015_2020=d_same["2015→2020"]))
        t4 = pd.DataFrame(rows4)
        t4.round(3).to_csv(OUT / "stats_vs_satellite_rate_by_forest_type.csv", index=False, encoding="utf-8-sig")
        print("\n(3-2) 흡수 속도 임상 보정\n", t4.round(2).to_string(index=False))


# ── 그림 ─────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in ["Malgun Gothic", "NanumGothic", "Noto Sans CJK KR", "Noto Sans KR"]:
    if any(f == x.name for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams["axes.unicode_minus"] = False
COL = {"전국": "#555555", "부산광역시": "#3A6FC0", "부산 금정구": "#C0703A", "부산 북구": "#8A4FB0"}
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
a1.plot(st.연도, st.평균AGB_t_ha, "o-", color="#4C8C6A", lw=2, label="위성 D (금정산, CCI v6)")
for k, (reg, g) in enumerate(df.groupby("지역", sort=False)):
    g = g.sort_values("연도"); v = g.ha당임목축적_m3ha
    base, lo, hi = (v * CASES[c].AGB환산_t_per_m3 for c in ("기본", "낮게", "높게"))
    x = g.연도 + 0.15 * k
    a1.errorbar(x, base, yerr=[base - lo, hi - base], fmt="s--" if len(g) > 1 else "s", color=COL.get(reg, "k"),
                capsize=3, label=f"통계 {reg}")
if TYPE_CSV.exists():
    for reg, nm in [("부산광역시 계", "부산광역시"), ("부산광역시 금정구", "부산 금정구")]:
        if reg in set(t3.지역):
            a1.plot(2020 + 0.6, t3.loc[t3.지역 == reg, "통계AGB_임상보정_t_ha"].item(), "*", ms=13, color=COL[nm],
                    label=f"통계 {nm} (임상 보정)")
a1.set_ylabel("지상부 바이오매스 (t/ha)"); a1.set_title("① 수준: ha당 지상부 바이오매스 (막대: 계수 낮게~높게)")
a1.legend(fontsize=7.5, loc="upper left")
labels, vals, cols, errs = ["위성 D\n2015–22 추세", "위성 D\n2015→20"], [D["대표"], d_same["2015→2020"]], ["#4C8C6A", "#7FB89A"], [None, None]
for _, rr in rate[rate.기간 == "2015→2020"].iterrows():
    labels.append(f"통계\n{rr.지역.replace('부산 ', '')}"); vals.append(rr.통계흡수_tC_ha_년_기본); cols.append(COL.get(rr.지역, "#888"))
    errs.append((rr.통계흡수_tC_ha_년_기본 - rr.통계흡수_tC_ha_년_낮게, rr.통계흡수_tC_ha_년_높게 - rr.통계흡수_tC_ha_년_기본))
xs = np.arange(len(vals))
a2.bar(xs, vals, color=cols, alpha=.85)
a2.errorbar([0], [D["대표"]], yerr=[[D["대표"] - D["하한"]], [D["상한"] - D["대표"]]], color="k", capsize=6)
for x_, e in zip(xs, errs):
    if e: a2.errorbar([x_], [vals[x_]], yerr=[[e[0]], [e[1]]], color="k", capsize=4)
if TYPE_CSV.exists() and TYPE15.exists():
    for x_, lab in zip(xs, labels):
        key = {"통계\n금정구": "부산광역시 금정구", "통계\n부산광역시": "부산광역시 계"}.get(lab)
        if key and key in set(t4.지역):
            a2.plot(x_, t4.loc[t4.지역 == key, "흡수_임상보정_tC_ha_년"].item(), "*", ms=13, color="k",
                    label="통계 (임상 보정)" if "임상" not in str(a2.get_legend_handles_labels()[1]) else None)
a2.axhspan(D["하한"], D["상한"], color="#4C8C6A", alpha=.08)
a2.set_xticks(xs); a2.set_xticklabels(labels, fontsize=8.5)
a2.set_ylabel("흡수 속도 (tC/ha/년)"); a2.set_title("② 흡수 속도 (연두 띠: 위성 90% 범위, ★: 임상 보정)")
a2.legend(fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(OUT / "stats_vs_satellite.png", dpi=150)
