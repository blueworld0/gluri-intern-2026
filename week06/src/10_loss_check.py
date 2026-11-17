"""6주차 5단계: 2022년 저장량 감소가 '실제로 숲이 사라져서'인지 확인 (Hansen 산림 손실 지도).

Hansen Global Forest Change v1.12 (30 m, 2001–2024 손실 연도) 를 경계 안에서 연도별로 집계.
  - 실제 벌채·고사·산불 등으로 나무가 사라진 면적이 2022년에 크게 늘었다면 → 실제 감소
  - 손실 면적이 평소와 비슷하면 → 위성 바이오매스 지도의 연도별 흔들림일 가능성이 큼
손실 면적 × 그해 평균 AGB 로 '사라진 바이오매스(최대치)'를 어림해 CCI 의 감소량과 비교.
※ 이 자료 범위는 경도 129.00° 부터라 경계의 서쪽 끝 약 4.5% 는 빠짐 (결과는 그만큼 조금 작게 나옴)
"""
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import shape

from config import *

LOSS = ROOT / "data_extra" / "lossyear.tif"      # 프로젝트 data/hansen_gfc_v1.12/lossyear.tif 를 복사
geom = shape(json.load(open(PROXY_BOUNDARY, encoding="utf-8"))["features"][0]["geometry"])
with rasterio.open(LOSS) as s:
    ly = s.read(1); tr = s.transform
inside = geometry_mask([geom], out_shape=ly.shape, transform=tr, invert=True)
lat = tr.f + tr.e * (np.arange(ly.shape[0]) + 0.5)
px_ha = np.broadcast_to(((abs(tr.a) * 111320 * np.cos(np.deg2rad(lat))) * (abs(tr.e) * 110574) / 1e4)[:, None], ly.shape)

st = pd.read_csv(OUT / "stock_by_year_national.csv")
rows = []
for y in range(2010, 2025):
    ha = px_ha[inside & (ly == y - 2000)].sum()
    row = dict(연도=y, 손실면적_ha=ha)
    if y in set(st.연도):
        agb = st.loc[st.연도 == y, "평균AGB_t_ha"].item()
        row["손실로_사라진_AGB_최대_t"] = ha * agb
        if y - 1 in set(st.연도):
            row["CCI_AGB_전년대비_t"] = st.loc[st.연도 == y, "AGB_총량_t"].item() - st.loc[st.연도 == y - 1, "AGB_총량_t"].item()
    rows.append(row)
t = pd.DataFrame(rows)
t.round(1).to_csv(OUT / "forest_loss_by_year.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 200)
print(t.round(1).to_string(index=False))
m = t[(t.연도 >= 2015) & (t.연도 <= 2021)].손실면적_ha.mean()
print(f"\n2015–2021 평균 손실 {m:.1f} ha/년, 2022 손실 {t.loc[t.연도 == 2022, '손실면적_ha'].item():.1f} ha")
