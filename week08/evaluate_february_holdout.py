"""Apply week07 footprint quality checks and evaluate the February date/orbit holdout."""

import argparse
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from rasterio.enums import Resampling
from rasterio.transform import array_bounds, rowcol
from rasterio.warp import transform
from shapely.geometry import Point, box
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "week06" / "step5"))
from process_candidate_a import crop, digest, read_grid, warp

APRIL = ROOT / "week07" / "landcover" / "forest_evaluation"
PRODUCT = "S2C_MSIL2A_20250128T021001_N0511_R103_T52SDE_20250128T060411"
ITEM = "S2C_52SDE_20250128_0_L2A"


def reflectance_calibration(path, bands, product):
    root = ET.parse(path).getroot()
    if root.findtext(".//PRODUCT_URI") != product:
        raise ValueError("Unexpected Sentinel-2 product metadata")
    scale = float(root.findtext(".//BOA_QUANTIFICATION_VALUE"))
    ids = {e.attrib["bandId"]: e.attrib["physicalBand"] for e in root.iter("Spectral_Information")}
    offsets = {ids[e.attrib["band_id"]]: float(e.text) for e in root.iter("BOA_ADD_OFFSET")}
    special = {
        e.findtext("SPECIAL_VALUE_TEXT"): int(e.findtext("SPECIAL_VALUE_INDEX"))
        for e in root.iter("Special_Values")
    }
    return scale, {b: offsets[b if b == "B8A" else "B" + str(int(b[1:]))] for b in bands}, special


def prepare_rasters(source, bounds, bands, product):
    scale, offsets, special = reflectance_calibration(source / "product_metadata.xml", bands, product)
    _, target = crop(source / "B04.jp2", bounds)
    if target["crs"].to_epsg() != 32652 or target["transform"].a != 10:
        raise ValueError("Unexpected 10 m reference grid")
    footprint_bounds = array_bounds(target["height"], target["width"], target["transform"])
    scl = read_grid(source / "SCL.jp2", target)
    classi = read_grid(source / "CLASSI_B00.jp2", target, [1, 2, 3])
    if not np.isin(classi, [0, 1]).all():
        raise ValueError("Unexpected CLASSI flags")
    arrays = {}
    for band in bands:
        raw, native = crop(source / f"{band}.jp2", footprint_bounds)
        quality = read_grid(source / f"QUALIT_{band}.jp2", native, list(range(1, 9)))
        detector = read_grid(source / f"DETFOO_{band}.jp2", native)
        if not np.isin(quality, [0, 1]).all() or not np.isin(detector, np.arange(13)).all():
            raise ValueError(f"Unexpected quality or detector flags for {band}")
        valid = (
            (raw != special["NODATA"]) & (raw != special["SATURATED"])
            & (quality == 0).all(axis=0) & (detector > 0)
        )
        reflectance = (raw.astype("float32") + offsets[band]) / scale
        reflectance[~valid] = np.nan
        aligned = warp(reflectance, native, target, Resampling.bilinear, nodata=np.nan)
        support = warp(valid.astype("float32"), native, target, Resampling.bilinear)
        aligned[support < 1 - 1e-6] = np.nan
        arrays[band] = aligned
    return target, scl, classi, arrays


def footprint_features(points, source, bands, product=PRODUCT + ".SAFE"):
    east, north = transform("EPSG:4326", "EPSG:32652", points.longitude.tolist(), points.latitude.tolist())
    bounds = (min(east) - 45, min(north) - 45, max(east) + 45, max(north) + 45)
    target, scl, classi, arrays = prepare_rasters(source, bounds, bands, product)
    grid = target["transform"]
    rows = []
    for point_id, x, y in zip(points.point_id, east, north):
        disk = Point(x, y).buffer(12.5, quad_segs=64)
        left, bottom, right, top = disk.bounds
        r0, c0 = rowcol(grid, left, top)
        r1, c1 = rowcol(grid, right, bottom)
        weights = np.zeros((r1 - r0 + 1, c1 - c0 + 1))
        for i, rr in enumerate(range(r0, r1 + 1)):
            for j, cc in enumerate(range(c0, c1 + 1)):
                west, north_edge = grid * (cc, rr)
                east_edge, south = grid * (cc + 1, rr + 1)
                weights[i, j] = disk.intersection(box(west, south, east_edge, north_edge)).area
        if not np.isclose(weights.sum(), disk.area, rtol=1e-10, atol=1e-7):
            raise ValueError(f"Footprint area mismatch for {point_id}")
        window = np.s_[r0:r1 + 1, c0:c1 + 1]
        values = {band: array[window] for band, array in arrays.items()}
        band_ok = np.isfinite(np.stack(list(values.values()))).all(axis=0)
        scl_ok = np.isin(scl[window], [4, 5, 6])
        classi_ok = (classi[:, window[0], window[1]] == 0).all(axis=0)
        valid = band_ok & scl_ok & classi_ok
        coverage = float(weights[valid].sum() / disk.area)
        passed = coverage >= 1 - 1e-6
        rr, cc = rowcol(grid, x, y)
        rows.append({
            "point_id": point_id, "easting_m": x, "northing_m": y,
            "original_scl_aligned_center": int(scl[rr, cc]),
            "footprint_pass": passed, "coverage": coverage,
            "band_coverage": float(weights[band_ok].sum() / disk.area),
            "scl_coverage": float(weights[scl_ok].sum() / disk.area),
            "classi_coverage": float(weights[classi_ok].sum() / disk.area),
            **{
                f"reflectance_{b}": float(np.average(a[valid], weights=weights[valid])) if passed else np.nan
                for b, a in values.items()
            },
        })
    return points.merge(pd.DataFrame(rows), on="point_id", validate="one_to_one")


def scores(actual, prediction):
    residual = actual - prediction
    return {
        "n": len(actual), "rmse": root_mean_squared_error(actual, prediction),
        "mae": mean_absolute_error(actual, prediction),
        "r2": r2_score(actual, prediction) if len(actual) > 1 else None,
        "mean_agbd_minus_prediction": float(residual.mean()),
        "underprediction_n": int((residual > 0).sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "source" / "january_selected")
    parser.add_argument("--scl-points", type=Path, default=HERE / "outputs" / "february_scl_points.csv")
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    args = parser.parse_args()
    train_path = APRIL / "rededge_forest184_input.csv"
    train = pd.read_csv(train_path)
    config_path = APRIL / "rededge_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    features = config["feature_sets"]["bands10"]
    bands = [feature.removeprefix("reflectance_") for feature in features]
    points = pd.read_csv(args.scl_points)
    points = points.loc[points.item_id.eq(ITEM) & points.scl_vegetation].copy()
    if len(points) != 60 or len(train) != 184 or points.point_id.isin(train.point_id).any():
        raise ValueError("Unexpected cohort or overlapping training and evaluation IDs")
    if not train.point_id.is_unique or not points.point_id.is_unique:
        raise ValueError("Point IDs must be unique within each cohort")
    data = footprint_features(points, args.source, bands)
    if not data.original_scl_aligned_center.eq(data.scl_aligned_10m).all():
        raise ValueError("Public COG and original SCL disagree at a selected point")
    usable = data.loc[data.footprint_pass].copy()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output_dir / "february_footprint_audit.csv", index=False, encoding="utf-8-sig")
    summary = {
        "scl_candidate_n": len(data), "n": len(usable),
        "high_agbd_n": int(usable.high_agbd_ge100.sum()),
        "train_n": len(train), "product": PRODUCT + ".SAFE",
        "source_sha256": {path.name: digest(path) for path in sorted(args.source.iterdir()) if path.is_file()},
        "training_sha256": digest(train_path), "config_sha256": digest(config_path),
        "metrics": {},
    }
    if len(usable) >= 2 and usable.high_agbd_ge100.any():
        model = RandomForestRegressor(**config["rf_params"])
        model.fit(train[features], train.agbd)
        usable["rf_bands10"] = model.predict(usable[features])
        usable["dummy_mean"] = train.agbd.mean()
        for scope, group in [("all", usable), ("agbd_ge100", usable.loc[usable.high_agbd_ge100])]:
            summary["metrics"][scope] = {
                name: scores(group.agbd_Mg_ha, group[name]) for name in ("rf_bands10", "dummy_mean")
            }
        distance = np.hypot(
            usable.easting_m.to_numpy()[:, None] - train.easting_m.to_numpy(),
            usable.northing_m.to_numpy()[:, None] - train.northing_m.to_numpy(),
        ).min(axis=1)
        usable["nearest_april_training_m"] = distance
        summary["less_than_500m_from_training_n"] = int((distance < 500).sum())
        summary["minimum_train_test_distance_m"] = float(distance.min())
        usable.to_csv(args.output_dir / "february_holdout_predictions.csv", index=False, encoding="utf-8-sig")
    payload = json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)
    (args.output_dir / "february_holdout_summary.json").write_text(payload + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if not k.endswith("sha256")}, indent=2))


if __name__ == "__main__":
    main()
