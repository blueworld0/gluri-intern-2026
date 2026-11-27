"""Compare bands and bands plus DSM elevation/slope on the fixed February holdout."""

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_origin, rowcol
from rasterio.warp import reproject
from sklearn.ensemble import RandomForestRegressor

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
APRIL = ROOT / "week07" / "landcover" / "forest_evaluation"
HOLDOUT = ROOT / "week08" / "outputs" / "february_holdout_predictions.csv"
DEM_TILES = [
    HERE / "source" / f"Copernicus_DSM_COG_10_N35_00_E{longitude}_00_DEM.tif"
    for longitude in (128, 129)
]
SOUTH_DEM = HERE / "source" / "Copernicus_DSM_COG_10_N34_00_E128_00_DEM.tif"
OUT = HERE / "outputs"
sys.path.insert(0, str(ROOT / "week08"))
from evaluate_february_holdout import digest, scores

RESOLUTION = 30


def sample_dsm(points):
    east, north = points.easting_m.to_numpy(), points.northing_m.to_numpy()
    xmin = np.floor(east.min() / RESOLUTION) * RESOLUTION - 3 * RESOLUTION
    ymax = np.ceil(north.max() / RESOLUTION) * RESOLUTION + 3 * RESOLUTION
    width = int(np.ceil((east.max() - xmin) / RESOLUTION)) + 3
    height = int(np.ceil((ymax - north.min()) / RESOLUTION)) + 3
    grid = from_origin(xmin, ymax, RESOLUTION, RESOLUTION)
    elevations = np.full((height, width), np.nan, dtype="float32")
    tiles = [*DEM_TILES]
    if SOUTH_DEM.is_file():
        tiles.append(SOUTH_DEM)
    for path in tiles:
        with rasterio.open(path) as source:
            if source.crs.to_epsg() != 4326 or source.count != 1:
                raise ValueError(f"Unexpected Copernicus DSM grid in {path}")
            reproject(
                rasterio.band(source, 1), elevations, dst_transform=grid, dst_crs="EPSG:32652",
                dst_nodata=np.nan, resampling=Resampling.bilinear, init_dest_nodata=False,
            )
    rr, cc = np.asarray(rowcol(grid, east, north))
    dz_dx = (elevations[rr, cc + 1] - elevations[rr, cc - 1]) / (2 * RESOLUTION)
    dz_dy = (elevations[rr - 1, cc] - elevations[rr + 1, cc]) / (2 * RESOLUTION)
    degrees = np.rad2deg(np.arctan(np.hypot(dz_dx, dz_dy)))
    result = pd.DataFrame({
        "point_id": points.point_id,
        "dsm_elevation_m": elevations[rr, cc],
        "dsm_slope_deg": degrees,
    })
    if not np.isfinite(result[["dsm_elevation_m", "dsm_slope_deg"]]).all().all():
        raise ValueError("DSM elevation or slope is missing for part of the fixed cohort")
    return result


def main():
    train_path = APRIL / "rededge_forest184_input.csv"
    config_path = APRIL / "rededge_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    train, test = pd.read_csv(train_path), pd.read_csv(HOLDOUT)
    if len(train) != 184 or len(test) != 35 or train.point_id.isin(test.point_id).any():
        raise ValueError("Week07 training or week08 evaluation points do not match expected counts")
    coords = ["point_id", "easting_m", "northing_m"]
    dsm = sample_dsm(pd.concat([train[coords], test[coords]], ignore_index=True))
    train = train.merge(dsm, on="point_id", validate="one_to_one")
    test = test.merge(dsm, on="point_id", validate="one_to_one")
    features = config["feature_sets"]["bands10"]
    topo_columns = ["dsm_elevation_m", "dsm_slope_deg"]
    test["dummy_mean_reproduced"] = train.agbd.mean()
    for name, columns in (("bands10", features), ("bands10_with_dsm", features + topo_columns)):
        model = RandomForestRegressor(**config["rf_params"])
        model.fit(train[columns], train.agbd)
        test[name] = model.predict(test[columns])
    np.testing.assert_allclose(test.bands10, test.rf_bands10, rtol=0, atol=1e-10)
    np.testing.assert_allclose(test.dummy_mean_reproduced, test.dummy_mean, rtol=0, atol=1e-10)
    metrics = {
        scope: {
            name: scores(frame.agbd_Mg_ha, frame[name])
            for name in ("dummy_mean", "bands10", "bands10_with_dsm")
        }
        for scope, frame in (("all", test), ("agbd_ge100", test.loc[test.high_agbd_ge100]))
    }
    split = {name: df[topo_columns].agg(["min", "max"]).to_dict()
             for name, df in (("april_training", train), ("february_holdout", test))}
    result = {
        "analysis": "exploratory DSM comparison on the existing February holdout",
        "train_n": len(train), "test_n": len(test),
        "terrain_range": split, "metrics": metrics,
        "sha256": {path.name: digest(path) for path in (*DEM_TILES, HOLDOUT, config_path, train_path)},
    }
    (OUT / "topography_comparison.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    dsm.to_csv(OUT / "topography_samples.csv", index=False, encoding="utf-8-sig")
    test[["point_id", "agbd_Mg_ha", "high_agbd_ge100", *topo_columns, "dummy_mean",
          "bands10", "bands10_with_dsm"]].to_csv(OUT / "topography_predictions.csv",
                                                  index=False, encoding="utf-8-sig")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
