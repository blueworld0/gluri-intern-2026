"""Physical intersections, rather than center sampling or unweighted zonal means."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import geopandas as gpd
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box
from kfcs.volume.raster import RasterProvider

def fixture(tmp_path, values, crs='EPSG:5179', transform=None, unit='m3_per_ha'):
    values=np.array(values,dtype='float64'); path=tmp_path/'volume.tif'
    with rasterio.open(path,'w',driver='GTiff',height=values.shape[0],width=values.shape[1],count=1,dtype=values.dtype,crs=crs,
                       transform=transform or from_origin(900000,1900010,10,10),nodata=-9999) as dst: dst.write(values,1)
    return {'volume':{'raster_path':str(path),'raster_unit':unit,'edition':'synthetic','band':1,'edge_segments':64},
            'units':{'m2_per_ha':10000},'spatial':{'crs':'EPSG:5179'}}

def frame(geoms,crs='EPSG:5179'):
    g=gpd.GeoDataFrame({'eligible':[True]*len(geoms)},geometry=geoms,crs=crs)
    g['area_ha']=g.to_crs('EPSG:5179').area/10000
    return g

def test_fractional_pixel_weights_and_order(tmp_path):
    c=fixture(tmp_path,[[100,300]])
    g=frame([box(900005,1900000,900020,1900010),box(900000,1900000,900010,1900010)])
    g.index=[9,2]
    r=RasterProvider(c).provide(g)
    assert r.index.tolist()==[9,2]
    assert r.m3_per_ha.tolist()==pytest.approx([700/3,100])

def test_nodata_partial_and_no_valid_pixels(tmp_path):
    c=fixture(tmp_path,[[100,-9999]])
    g=frame([box(900005,1900000,900020,1900010),box(900010,1900000,900020,1900010),box(910000,1900000,910010,1900010)])
    p=RasterProvider(c);r=p.provide(g)
    assert r.iloc[0].m3_per_ha==100
    assert r.iloc[0].raster_valid_area_ha==pytest.approx(.005)
    assert r.iloc[1:].m3_per_ha.isna().all()
    assert r.volume_flag.tolist()==['직접','결측','결측']
    assert p.audit['missing_area_ha']==pytest.approx(.02)

def test_m3_per_pixel(tmp_path):
    c=fixture(tmp_path,[[1,3]],unit='m3_per_pixel')
    r=RasterProvider(c).provide(frame([box(900005,1900000,900020,1900010)]))
    assert r.m3_per_ha.iloc[0]==pytest.approx(700/3)

def test_polygon_crs_differs(tmp_path):
    c=fixture(tmp_path,[[100,300]])
    g=frame([box(900005,1900000,900020,1900010)]).to_crs('EPSG:4326')
    r=RasterProvider(c).provide(g)
    assert r.m3_per_ha.iloc[0]==pytest.approx(700/3,rel=1e-6)

def test_geographic_raster_physical_pixel_area(tmp_path):
    c=fixture(tmp_path,[[1,1]],crs='EPSG:4326',transform=from_origin(129,35.01,.01,.01),unit='m3_per_pixel')
    g=frame([box(129,35,129.01,35.01)],crs='EPSG:4326')
    r=RasterProvider(c).provide(g)
    # One cubic metre over a whole geographic pixel; no degree²→ha shortcut.
    expected=1/g.area_ha.iloc[0]
    assert r.m3_per_ha.iloc[0]==pytest.approx(expected,rel=1e-4)

def test_excluded_forest_type_is_zero(tmp_path):
    c=fixture(tmp_path,[[-9999]])
    g=frame([box(900000,1900000,900010,1900010)]);g['eligible']=False
    result=RasterProvider(c).provide(g)
    assert result.m3_per_ha.iloc[0]==0
    assert result.volume_source.iloc[0]=='raster:synthetic'
