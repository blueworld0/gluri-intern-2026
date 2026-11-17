"""Verify exported means and approximate circle weights using quadrature."""
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import rasterio
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'forest_evaluation'


def exact_circle_rectangle(left,bottom,right,top,r=12.5):
    lo,hi=max(left,-r),min(right,r)
    if lo>=hi:
        return 0.
    points=[0.]
    for y in [bottom,top]:
        if abs(y)<r:
            edge=math.sqrt(r*r-y*y)
            points.extend([-edge,edge])
    def height(x):
        edge=math.sqrt(max(0.,r*r-x*x))
        return max(0.,min(top,edge)-max(bottom,-edge))
    return quad(height,lo,hi,points=sorted(set(x for x in points if lo<x<hi)),epsabs=1e-8,epsrel=1e-10)[0]


def main():
    assert np.isclose(exact_circle_rectangle(-20,-20,20,20),math.pi*12.5**2,atol=1e-7)
    assert np.isclose(exact_circle_rectangle(0,0,20,20),math.pi*12.5**2/4,atol=1e-7)
    assert exact_circle_rectangle(20,20,30,30)==0
    weights=pd.read_csv(OUT/'footprint_pixel_weights.csv')
    averaged=pd.read_csv(OUT/'footprint_forest184_input.csv').set_index('point_id')
    assert weights.common_valid.all() and len(averaged)==184
    for b in ['B02','B03','B04','B08','B11','B12']:
        mean=(weights[b]*weights.intersection_area_m2).groupby(weights.point_id).sum()/weights.intersection_area_m2.groupby(weights.point_id).sum()
        assert np.allclose(mean.sort_index(),averaged['reflectance_'+b].sort_index(),atol=1e-12,rtol=1e-12)
    blue,red,nir,swir=[averaged['reflectance_'+b] for b in ['B02','B04','B08','B11']]
    assert np.allclose(averaged.ndvi,(nir-red)/(nir+red),atol=1e-12,rtol=1e-12)
    assert np.allclose(averaged.evi,2.5*(nir-red)/(nir+6*red-7.5*blue+1),atol=1e-12,rtol=1e-12)
    assert np.allclose(averaged.ndmi,(nir-swir)/(nir+swir),atol=1e-12,rtol=1e-12)
    # Deterministic, outcome-independent subset; quadrature uses the true
    # circle equation, without Shapely or the polygon intersection algorithm.
    selected=sorted(averaged.index)[::9]
    errors=[]
    with rasterio.open(ROOT/'week06/step5/outputs/rasters/reflectance_B04.tif') as ds:
        t=ds.transform
        for row in weights[weights.point_id.isin(selected)].itertuples(index=False):
            p=averaged.loc[row.point_id]
            left,top=t*(row.raster_col,row.raster_row)
            right,bottom=t*(row.raster_col+1,row.raster_row+1)
            area=exact_circle_rectangle(left-p.easting_m,bottom-p.northing_m,right-p.easting_m,top-p.northing_m)
            errors.append(abs(area-row.intersection_area_m2))
    # The 256-edge inscribed polygon omits <0.05 m2 from the full circle.
    assert max(errors)<.051
    result=dict(status='passed',all184_exported_reflectance_means_recomputed=True,
        all184_indices_recomputed=True,independent_geometry_method='True circle/rectangle numerical integration; no Shapely',
        sampled_circle_n=len(selected),sampled_pixel_n=len(errors),
        max_true_circle_vs_polygon_pixel_area_difference_m2=max(errors),area_tolerance_m2=.051,
        analytic_whole_circle_quadrant_and_disjoint_checks=True)
    (OUT/'footprint_independent_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
