"""Independent checks of generated CV predictions, masks, and carbon accounting."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from sklearn.metrics import r2_score,mean_absolute_error

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'data/processed'
cv=pd.read_csv(P/'model_cv.csv')
assert len(cv)==12, 'Expected 3 thresholds x 2 models x 2 CV schemes'
source=pd.read_csv(P/'gedi_l4a_2021_filtered.csv',dtype={'shot_number':'uint64'})
ids=set(source.shot_number.tolist())
for row in cv.itertuples():
    path=P/f'{row.model}_{row.split}_{row.threshold_Mg_ha}_oof.csv'
    d=pd.read_csv(path,dtype={'shot_number':'uint64'})
    assert len(d)==row.n and d.shot_number.is_unique and set(d.shot_number)<=ids
    assert d['fold'].nunique()==5 and np.isfinite(d.prediction).all()
    if row.split=='spatial':
        assert d.groupby('block_id')['fold'].nunique().max()==1
    error=d.prediction-d.agbd
    assert np.isclose(row.R2,r2_score(d.agbd,d.prediction),atol=1e-9)
    assert np.isclose(row.RMSE,np.sqrt(np.mean(error**2)),atol=1e-9)
    assert np.isclose(row.Bias,error.mean(),atol=1e-9)
    assert np.isclose(row.MAE,mean_absolute_error(d.agbd,d.prediction),atol=1e-9)
with rasterio.open(P/'features_10m.tif') as features:
    assert features.count==31 and all(features.descriptions)
    assert features.crs.to_epsg()==32652 and features.res==(10.,10.)
    assert set(features.dtypes)=={'float32'}
carbon=pd.read_csv(P/'carbon_totals.csv')
checks=[]
for row in carbon.itertuples():
    count,total=0,0.
    wc_path=ROOT/'data/raw/worldcover/ESA_WorldCover_10m_2021_v200_N36E126_Map.tif'
    with rasterio.open(P/f'agb_map_{row.model}_10m.tif') as ds,rasterio.open(wc_path) as wc:
        with WarpedVRT(wc,crs=ds.crs,transform=ds.transform,width=ds.width,height=ds.height,
                       resampling=Resampling.nearest,warp_mem_limit=16) as aligned:
            for _,win in ds.block_windows(1):
                arr=ds.read(1,window=win,masked=True)
                valid=~np.ma.getmaskarray(arr)
                assert np.isfinite(arr.data[valid]).all() and (arr.data[valid]>=0).all()
                assert (aligned.read(1,window=win)[valid]==10).all()
                count+=int(valid.sum())
                total+=float(arr.data[valid].sum(dtype='float64'))*.01
    assert np.isclose(count*.01,row.mapped_tree_area_ha,atol=1e-8)
    assert np.isclose(total,row.AGB_total_Mg,rtol=1e-9)
    assert np.isclose(total*1.25*.47,row.carbon_tC,rtol=1e-9)
    assert np.isclose(row.carbon_tC*44/12,row.CO2_equivalent_t,rtol=1e-9)
    assert np.isclose(row.tree_cover_area_ha-row.mapped_tree_area_ha,row.unmapped_tree_area_ha)
    checks.append({'model':row.model,'checked_map_pixels':count,'AGB_sum_Mg':total})
selected=cv.loc[(cv.threshold_Mg_ha==600)&(cv.split=='spatial')].sort_values(['RMSE','model']).iloc[0].model
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while block:=f.read(1024*1024):
            h.update(block)
    return h.hexdigest()
assert sha(P/'agb_map_10m.tif')==sha(P/f'agb_map_{selected}_10m.tif')
for filename in ['obs_vs_pred.png','feature_importance.png','agb_map.png']:
    assert (ROOT/'figures'/filename).stat().st_size>10000
assert (ROOT/'step2_summary.md').is_file()
(P/'step2_verification.json').write_text(json.dumps({'status':'passed','cv_rows':len(cv),
    'selected_model':selected,'map_checks':checks},ensure_ascii=False,indent=2),encoding='utf-8')
print('PASS: 12 CV settings, unique uint64 IDs, spatial block isolation, 31-band grid, tree masks, raster sums, carbon arithmetic, selected map and figures')
