"""6주차 1단계: 국가 고유계수(WD·BEF·R)를 금정산 수종 구성으로 가중평균.

왜 하나?
  5주차 통계 비교에서 쓴 계수(WD·BEF·R)는 '침엽/활엽 평균'의 임시값이었다.
  이번에는 출처가 있는 국가 고유계수 표(이선정 외 2015, Table 1)를 쓰고,
  임상도에서 금정산에 실제로 있는 수종 면적만큼 가중평균한다.

순서
  1) 임상도(2019년판, 부산·양산)를 경계로 잘라 수종 코드별 면적을 구한다
  2) 수종 코드마다 국가계수(WD, BEF, R)를 붙인다 (config.SPECIES_MAP)
  3) 면적가중 평균 → 부피(m³)를 바이오매스·탄소로 바꾸는 환산계수
  4) 위성 계산용 칸별 R 지도를 새 R 값으로 다시 그린다
"""
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

from config import *

OUT.mkdir(parents=True, exist_ok=True)

# ── 1. 수종 코드 → 계수 (혼효림은 섞어서) ─────────────────────────
def coef_of(code):
    name, keys = SPECIES_MAP[code]
    m = np.array(MIX_WEIGHT.get(code, [1 / len(keys)] * len(keys)))
    wd = np.array([NATIONAL[k][0] for k in keys]); bef = np.array([NATIONAL[k][1] for k in keys])
    r = np.array([NATIONAL[k][2] for k in keys])
    con = sum(mi for mi, k in zip(m, keys) if NATIONAL[k][3] == "침엽")   # 침엽 몫 (0~1)
    return dict(수종=name, 적용계수="+".join(dict.fromkeys(keys)) if len(keys) > 1 else keys[0],
                WD=(m * wd).sum(), BEF=(m * bef).sum(), R=(m * r).sum(),
                WDxBEF=(m * wd * bef).sum(), WDxBEFx1R=(m * wd * bef * (1 + r)).sum(), 침엽몫=con)

# ── 2. 임상도 자르기 ───────────────────────────────────────────
b = gpd.read_file(PROXY_BOUNDARY).to_crs(5179)
parts = []
for src in IMSANG:
    g = gpd.read_file(src, encoding="cp949", bbox=tuple(b.total_bounds))
    parts.append(gpd.clip(g[["FRTP_CD", "KOFTR_GROU", "geometry"]], b))
g = pd.concat(parts, ignore_index=True)
g["ha"] = g.area / 1e4
g["KOFTR_GROU"] = g.KOFTR_GROU.astype(str)

area = g.groupby("KOFTR_GROU").ha.sum()
rows = []
for code, ha in area.sort_values(ascending=False).items():
    if code in SPECIES_MAP:
        rows.append(dict(KOFTR_GROU=code, 면적_ha=ha, **coef_of(code)))
    else:
        rows.append(dict(KOFTR_GROU=code, 면적_ha=ha, 수종="(비산림/죽림 등)"))
t = pd.DataFrame(rows)
t["비율_pct"] = 100 * t.면적_ha / t.면적_ha.sum()
forest = t.dropna(subset=["WD"]).copy()
w = forest.면적_ha / forest.면적_ha.sum()

# ── 3. 면적가중 평균과 환산계수 ──────────────────────────────────
agb_f = float((w * forest.WDxBEF).sum())               # 1 m³ → 지상부 바이오매스 (t)
c_f = float((w * forest.WDxBEFx1R).sum() * CF)         # 1 m³ → 탄소 (tC)
r_mean = float((w * forest.R).sum())

# 낮게/높게: 수종 구성이 임상도와 조금 다를 때.
#   침엽 몫은 침엽수 중 가장 낮은(높은) 종, 활엽 몫은 활엽수 중 가장 낮은(높은) 종으로 바꿔 계산
nat = pd.DataFrame(NATIONAL, index=["WD", "BEF", "R", "임상"]).T
nat[["WD", "BEF", "R"]] = nat[["WD", "BEF", "R"]].astype(float)
nat["WDxBEF"] = nat.WD * nat.BEF
nat["WDxBEFx1R"] = nat.WDxBEF * (1 + nat.R)
con_share = float((w * forest.침엽몫).sum())
def bound(col, fn):
    return con_share * fn(nat.loc[nat.임상 == "침엽", col]) + (1 - con_share) * fn(nat.loc[nat.임상 == "활엽", col])
fac = pd.DataFrame([
    dict(경우="기본 (임상도 수종별 가중)", AGB환산_t_per_m3=agb_f, 탄소환산_tC_per_m3=c_f, 평균R=r_mean),
    dict(경우="낮게 (침엽·활엽 각각 표의 최저 수종)", AGB환산_t_per_m3=bound("WDxBEF", min),
         탄소환산_tC_per_m3=bound("WDxBEFx1R", min) * CF, 평균R=np.nan),
    dict(경우="높게 (침엽·활엽 각각 표의 최고 수종)", AGB환산_t_per_m3=bound("WDxBEF", max),
         탄소환산_tC_per_m3=bound("WDxBEFx1R", max) * CF, 평균R=np.nan),
])
# 5주차 임시값으로 같은 계산 (비교)
o = OLD_WEEK05; s = sum(o["SHARE"].values())
old_agb = sum(o["SHARE"][k] / s * o["WD"][k] * o["BEF"][k] for k in o["SHARE"])
old_c = sum(o["SHARE"][k] / s * o["WD"][k] * o["BEF"][k] * (1 + o["R"][k]) for k in o["SHARE"]) * CF
fac.loc[len(fac)] = dict(경우="(참고) 5주차 임시값", AGB환산_t_per_m3=old_agb, 탄소환산_tC_per_m3=old_c, 평균R=np.nan)

t.round(4).to_csv(OUT / "species_coef_table.csv", index=False, encoding="utf-8-sig")
fac.round(4).to_csv(OUT / "conversion_factors_national.csv", index=False, encoding="utf-8-sig")
nat.round(3).to_csv(OUT / "national_coef_source_table.csv", encoding="utf-8-sig", index_label="수종")
json.dump(dict(출처="이선정 외 2015 한국기후변화학회지 6(4):303-310 Table 1 (국가 고유계수)",
               임상도="2019년판(부산·양산)", 산림면적_ha=float(forest.면적_ha.sum()),
               침엽몫=con_share, 면적가중_WD=float((w * forest.WD).sum()), 면적가중_BEF=float((w * forest.BEF).sum()),
               면적가중_R=r_mean, 이전_R=OLD_R_MEAN, AGB환산_t_per_m3=agb_f, 탄소환산_tC_per_m3=c_f),
          open(OUT / "national_coef_summary.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# ── 4. 새 R 지도 (WorldCover 격자) ────────────────────────────────
r_by_code = {c: coef_of(c)["R"] for c in SPECIES_MAP}
gf = g[g.KOFTR_GROU.isin(r_by_code)].copy()
gf["R"] = gf.KOFTR_GROU.map(r_by_code)
with rasterio.open(WORLDCOVER) as wc:
    prof, tr, shp, crs = wc.profile.copy(), wc.transform, (wc.height, wc.width), wc.crs
gf = gf.to_crs(crs)
rmap = rasterize(zip(gf.geometry, gf.R), out_shape=shp, transform=tr, fill=0, dtype="float32")
prof.update(dtype="float32", nodata=0, compress="deflate")
with rasterio.open(R_MAP_NEW, "w", **prof) as dst:
    dst.write(rmap, 1)

pd.set_option("display.width", 220)
print(t.head(12).round(3).to_string(index=False))
print(f"\n침엽 몫 {con_share:.3f} / 면적가중 R {r_mean:.4f} (이전 {OLD_R_MEAN})")
print(fac.round(3).to_string(index=False))
print("저장:", R_MAP_NEW)
