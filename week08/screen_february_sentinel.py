"""Screen the fixed February forest cohort using the catalogued L2A SCL rasters."""

import argparse
from pathlib import Path
import json

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform
from shapely.geometry import Point, shape

HERE = Path(__file__).resolve().parent


def screen(item, forest, source):
    path = source / f"{item['id']}.tif"
    footprint = shape(item["geometry"])
    covered = [
        footprint.covers(Point(lon, lat))
        for lon, lat in zip(forest.longitude, forest.latitude)
    ]
    with rasterio.open(path) as raster:
        if raster.res != (20, 20) or raster.count != 1:
            raise ValueError(f"Expected original 20 m SCL grid: {path}")
        x, y = transform("EPSG:4326", raster.crs, forest.longitude.tolist(), forest.latitude.tolist())
        vals = np.array([sample[0] for sample in raster.sample(zip(x, y))])
        pixel10 = raster.transform * rasterio.Affine.scale(0.5)
        col10, row10 = (~pixel10) * (np.array(x), np.array(y))
        center_x, center_y = pixel10 * (np.floor(col10) + 0.5, np.floor(row10) + 0.5)
        aligned = np.array([sample[0] for sample in raster.sample(zip(center_x, center_y))])
    if not np.isin(vals, np.arange(12)).all() or not np.array_equal(vals, aligned):
        raise ValueError(f"Unexpected SCL codes or nearest-neighbour grid mismatch: {path}")
    result = forest[["point_id", "latitude", "longitude", "agbd_Mg_ha", "high_agbd_ge100"]].copy()
    result["item_id"] = item["id"]
    result["product"] = item["properties"]["s2:product_uri"]
    result["sensing_utc"] = pd.to_datetime(result["product"].iloc[0].split("_")[2], utc=True)
    result["offset_hours"] = (
        result.sensing_utc - pd.to_datetime(forest.observed_at_utc, format="ISO8601", utc=True)
    ).dt.total_seconds() / 3600
    result["catalog_footprint"] = covered
    result["scl_native_20m"] = vals
    result["scl_aligned_10m"] = aligned
    result["scl_vegetation"] = result.catalog_footprint & (vals == 4) & (aligned == 4)
    result["scl_cloud_shadow"] = result.catalog_footprint & np.isin(vals, [2, 3, 8, 9, 10])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    args = parser.parse_args()
    source = HERE / "source" / "scl"
    forest = pd.read_csv(HERE / "outputs" / "february_landcover_points.csv")
    forest = forest.loc[forest.landcover_status.eq("forest")].copy()
    if len(forest) != 220 or not forest.point_id.is_unique:
        raise ValueError("Expected 220 distinct forest candidates")
    catalog = json.loads((HERE / "outputs" / "february_earth_search.json").read_text(encoding="utf-8-sig"))
    items = catalog["features"]
    if not items or len(items) != catalog["context"]["matched"]:
        raise ValueError("Missing or incomplete search result")
    points = pd.concat([screen(item, forest, source) for item in items], ignore_index=True)
    points["scl_vegetation_high_agbd"] = points.scl_vegetation & points.high_agbd_ge100
    summary = points.groupby(["item_id", "product", "sensing_utc"], as_index=False).agg(
        footprint_n=("catalog_footprint", "sum"),
        scl_vegetation_n=("scl_vegetation", "sum"),
        scl_vegetation_high_agbd_n=("scl_vegetation_high_agbd", "sum"),
        scl_cloud_shadow_n=("scl_cloud_shadow", "sum"),
        median_offset_hours=("offset_hours", "median"),
    )
    summary["abs_offset_hours"] = summary.median_offset_hours.abs()
    summary = summary.sort_values(
        ["scl_vegetation_n", "abs_offset_hours", "item_id"], ascending=[False, True, True]
    ).drop(columns="abs_offset_hours")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    points.to_csv(args.output_dir / "february_scl_points.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(args.output_dir / "february_scl_products.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
