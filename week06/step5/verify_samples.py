"""Independent source-pixel and manual bilinear checks of emitted CSV values."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

HERE=Path(__file__).resolve().parent


def run():
    p=pd.read_csv(HERE/"outputs/gedi_sentinel_all_points.csv",dtype={"shot_number":"string"})
    safe=next((HERE/"source").glob("*.SAFE"))
    granule=next((safe/"GRANULE").iterdir())
    errors={}
    for band in ["B02","B03","B04","B08"]:
        source=next((granule/"IMG_DATA/R10m").glob(f"*_{band}_10m.jp2"))
        with rasterio.open(source) as ds:
            raw=np.array([v[0] for v in ds.sample(zip(p.easting_m,p.northing_m))])
        np.testing.assert_array_equal(raw,p[f"{band}_native_dn"])
        expected=(raw.astype("float64")-1000)/10000
        err=float(np.max(abs(expected-p[f"reflectance_{band}"])))
        assert err<1e-6
        errors[band]=err
    with rasterio.open(HERE/"outputs/rasters/ndvi.tif") as ds:transform=ds.transform
    sample=p.sample(30,random_state=42)
    for band in ["B11","B12"]:
        source=next((granule/"IMG_DATA/R20m").glob(f"*_{band}_20m.jp2"))
        diffs=[]
        with rasterio.open(source) as ds:
            for row in sample.itertuples():
                x,y=transform*(row.raster_col+.5,row.raster_row+.5)
                col,r=(~ds.transform)*(x,y)
                col-=.5;r-=.5
                c0,r0=int(np.floor(col)),int(np.floor(r))
                a=ds.read(1,window=Window(c0,r0,2,2)).astype("float64")
                assert ((a!=0)&(a!=65535)).all()
                a=(a-1000)/10000
                u,v=col-c0,r-r0
                expected=(1-v)*((1-u)*a[0,0]+u*a[0,1])+v*((1-u)*a[1,0]+u*a[1,1])
                diffs.append(abs(expected-getattr(row,f"reflectance_{band}")))
        errors[band]=float(max(diffs));assert errors[band]<1e-6
    final=pd.read_csv(HERE/"outputs/vegetation_candidates.csv")
    assert len(final)==441 and final.scl_native_20m.eq(4).all() and final.scl_aligned_10m.eq(4).all()
    assert not final.duplicate_pixel.any()
    assert final.forest_status.eq("not_assessed").all()
    report=(HERE/"outputs/report.md").read_text(encoding="utf-8")
    assert "\ufffd" not in report and "???" not in report
    validation=json.loads((HERE/"outputs/validation.json").read_text())
    validation.update(native_10m_source_checks_points_per_band=1055,manual_20m_bilinear_points_per_band=30,
                      independent_reflectance_max_absolute_errors=errors,mask_regression_tests_passed=2,
                      utf8_report_verified=True)
    (HERE/"outputs/validation.json").write_text(json.dumps(validation,indent=2),encoding="utf-8")
    print(json.dumps(validation,indent=2))


if __name__=="__main__":run()
