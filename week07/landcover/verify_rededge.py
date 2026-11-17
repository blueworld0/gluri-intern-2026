"""Independently reproduce bilinear extra-band samples from native 20m DN."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'forest_evaluation'


def main():
    config=json.loads((OUT/'rededge_config.json').read_text(encoding='utf-8'))
    pixels=pd.read_csv(OUT/'rededge_pixel_values.csv')
    data=pd.read_csv(OUT/'rededge_forest184_input.csv').set_index('point_id')
    granule=next((ROOT/'week06/step5/source').glob('*.SAFE/GRANULE/*'))
    with rasterio.open(ROOT/'week06/step5/outputs/rasters/reflectance_B04.tif') as ref:
        x=ref.transform.c+(pixels.raster_col.to_numpy()+.5)*ref.transform.a
        y=ref.transform.f+(pixels.raster_row.to_numpy()+.5)*ref.transform.e
    checks=[]
    for b in ['B05','B06','B07','B8A']:
        raw_path=next((granule/'IMG_DATA/R20m').glob(f'*_{b}_20m.jp2'))
        with rasterio.open(raw_path) as ds:
            t=ds.transform
            fractional_col=(x-t.c)/t.a-.5
            fractional_row=(y-t.f)/t.e-.5
            c=np.floor(fractional_col).astype(int)
            r=np.floor(fractional_row).astype(int)
            fx,fy=fractional_col-c,fractional_row-r
            left,top=c.min(),r.min()
            win=Window(int(left),int(top),int(c.max()-left+2),int(r.max()-top+2))
            raw=ds.read(1,window=win)
            rr,cc=r-top,c-left
        with rasterio.open(granule/f'QI_DATA/MSK_QUALIT_{b}.jp2') as ds:
            assert ds.transform==t
            qi=ds.read(list(range(1,9)),window=win)
        with rasterio.open(granule/f'QI_DATA/MSK_DETFOO_{b}.jp2') as ds:
            assert ds.transform==t
            detector=ds.read(1,window=win)
        manual=np.zeros(len(pixels))
        for dr,dc,w in [(0,0,(1-fx)*(1-fy)),(0,1,fx*(1-fy)),(1,0,(1-fx)*fy),(1,1,fx*fy)]:
            dn=raw[rr+dr,cc+dc]
            assert (~np.isin(dn,[0,65535])).all()
            assert (qi[:,rr+dr,cc+dc]==0).all()
            assert np.isin(detector[rr+dr,cc+dc],np.arange(1,13)).all()
            rho=(dn.astype('float32')+config['new_band_offsets'][b])/config['quantification']
            manual+=rho.astype(float)*w
        error=float(np.max(abs(manual-pixels[b].to_numpy())))
        assert error<1e-7
        means=(pd.Series(manual)*pixels.intersection_area_m2).groupby(pixels.point_id).sum()/pixels.intersection_area_m2.groupby(pixels.point_id).sum()
        mean_error=float(np.max(abs(means.sort_index()-data['reflectance_'+b].sort_index())))
        assert mean_error<1e-7
        checks.append(dict(band=b,aligned_pixel_n=len(pixels),native_contributors_all_valid=True,
            maximum_manual_bilinear_difference=error,maximum_manual_circle_mean_difference=mean_error))
    pred=pd.read_csv(OUT/'rededge_predictions.csv')
    metrics=pd.read_csv(OUT/'rededge_overall_metrics.csv').set_index('model')
    assignments=pd.read_csv(OUT/'forest184_split_assignments.csv')
    for name,g in pred.groupby('model'):
        assert len(g)==184 and g.point_id.is_unique
        residual=g.agbd-g.prediction
        assert np.isclose(np.sqrt(np.mean(residual**2)),metrics.loc[name,'rmse'],atol=1e-12,rtol=1e-12)
        for fold,h in g.groupby('fold'):
            roles=assignments[assignments.split_id.eq(f'new441_spatial5_fold{fold}')]
            assert set(h.point_id)==set(roles.loc[roles.role.eq('test'),'point_id'])
    result=dict(status='passed',method='Explicit four-native-pixel bilinear arithmetic, no GDAL warp or reused processing helper',
        band_checks=checks,metric_and_original_test_ids_verified=True,reflectance_tolerance=1e-7)
    (OUT/'rededge_independent_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
