"""Prioritize NASA V003 originals using the GFZ public TileDB records."""

import argparse
import json
import sys

import numpy as np
import pandas as pd

from inventory_spring import AOI, HERE, OUT, classify_points, load_polygons, load_search

CACHE = HERE / "source" / "gfz_multiyear_spring.csv"
ARRAY = "s3://dog-ext--gedidb--gedi-l2-l4-v003/array_uri"
COLUMNS = [
    "shot_number", "timestamp_ns", "agbd", "agbd_se",
    "l2a_quality_flag_rel3", "l4a_quality_flag_rel3", "l2_algrunflag", "degrade_flag",
]


def fetch_gfz():
    sys.path.append(str(HERE / "source" / "optional_python"))
    import tiledb

    cfg = tiledb.Config({
        "vfs.s3.endpoint_override": "s3.gfz.de",
        "vfs.s3.region": "eu-central-1",
        "vfs.s3.scheme": "https",
        "vfs.s3.use_virtual_addressing": "false",
        "vfs.s3.no_sign_request": "true",
        "vfs.s3.connect_timeout_ms": "10000",
        "vfs.s3.request_timeout_ms": "45000",
        "vfs.s3.connect_max_tries": "1",
        "sm.compute_concurrency_level": "4",
        "sm.io_concurrency_level": "4",
    })
    west, south, east, north = AOI
    with tiledb.open(ARRAY, ctx=tiledb.Ctx(cfg)) as array:
        raw = array.query(attrs=COLUMNS, coords=True).multi_index[south:north, west:east]
    frame = pd.DataFrame(raw)
    frame["observed_at_utc"] = pd.to_datetime(frame.timestamp_ns, utc=True, unit="ns")
    frame = frame.loc[
        frame.observed_at_utc.dt.month.isin([4, 5])
        & frame.latitude.lt(north) & frame.longitude.lt(east)
    ].copy()
    frame["orbit"] = "O" + frame.shot_number.astype(str).str.zfill(18).str[:5]
    frame.to_csv(CACHE, index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--sentinel", action="store_true", help="check SCL at promising GFZ forest points")
    parser.add_argument("--footprint-orbit", help="also check original imagery for one candidate")
    args = parser.parse_args()
    if args.fetch or not CACHE.is_file():
        fetch_gfz()
    frame = pd.read_csv(CACHE, dtype={"shot_number": "string"})
    frame["date_utc"] = pd.to_datetime(frame.observed_at_utc, format="ISO8601", utc=True).dt.strftime("%Y-%m-%d")
    inventory = pd.read_csv(OUT / "multiyear_spring_gedi_inventory.csv")
    stages = frame.loc[
        np.isfinite(frame.agbd) & frame.agbd.ge(0)
        & np.isfinite(frame.agbd_se) & frame.agbd_se.ge(0)
        & frame.l2a_quality_flag_rel3.eq(1) & frame.l4a_quality_flag_rel3.eq(1)
        & frame.l2_algrunflag.eq(1)
    ].copy()
    stages = stages.join(classify_points(stages, *load_polygons()))
    summary = inventory.merge(
        frame.groupby(["date_utc", "orbit"]).size().rename("gfz_reported_n"),
        left_on=["date_utc", "orbit"], right_index=True, how="left", validate="one_to_one",
    )
    stages["forest"] = stages.landcover_status.eq("forest")
    stages["forest_high_agbd"] = stages.forest & stages.agbd.ge(100)
    counts = stages.groupby(["date_utc", "orbit"]).agg(
        gfz_partial_quality_n=("agbd", "size"),
        gfz_forest_n=("forest", "sum"),
        gfz_forest_high_agbd_n=("forest_high_agbd", "sum"),
    )
    summary = summary.merge(counts, left_on=["date_utc", "orbit"], right_index=True,
                            how="left", validate="one_to_one")
    group_cols = ["gfz_partial_quality_n", "gfz_forest_n", "gfz_forest_high_agbd_n"]
    summary.loc[summary.gfz_reported_n.notna(), group_cols] = (
        summary.loc[summary.gfz_reported_n.notna(), group_cols].fillna(0)
    )
    if args.sentinel or args.footprint_orbit:
        from evaluate_spring import APRIL, choose_product, compatible_products, footprint_features, source_files

        candidates = pd.read_csv(OUT / "multiyear_spring_sentinel_candidates.csv")
        _, catalog = load_search(True)
        promising = summary.loc[
            summary.role.eq("holdout_candidate") & ~summary.source_present
            & summary.gfz_forest_n.ge(2)
        ]
        for index, candidate in promising.iterrows():
            orbit = candidate.orbit
            points = stages.loc[stages.orbit.eq(orbit) & stages.forest].copy()
            points["point_id"] = points.shot_number
            points["agbd_Mg_ha"] = points.agbd
            points["high_agbd_ge100"] = points.forest_high_agbd
            item_ids = set(candidates.loc[candidates.orbit.eq(orbit), "item"])
            products = compatible_products(catalog, item_ids)
            if not products:
                continue
            best, counts, selected = choose_product(points, products)
            counts.to_csv(OUT / f"multiyear_precheck_{orbit}_scl.csv", index=False)
            summary.loc[index, "gfz_best_scl_item"] = best
            summary.loc[index, "gfz_forest_scl_n"] = len(selected)
            summary.loc[index, "gfz_forest_high_agbd_scl_n"] = int(selected.high_agbd_ge100.sum())
            print(f"{orbit}: {best} forest SCL = {len(selected)}", flush=True)
            if args.footprint_orbit == orbit and not selected.empty:
                item = next(product for product in products if product["id"] == best)
                config = json.loads((APRIL / "rededge_config.json").read_text(encoding="utf-8"))
                bands = [feature.removeprefix("reflectance_") for feature in config["feature_sets"]["bands10"]]
                product = item["properties"]["s2:product_uri"]
                raw = source_files(item, bands)
                audited = footprint_features(selected, raw, bands, product)
                if not audited.original_scl_aligned_center.eq(audited.scl_aligned_10m).all():
                    raise ValueError("Public COG and original SCL disagree at a GFZ candidate point")
                usable = audited.loc[audited.footprint_pass]
                summary.loc[index, "gfz_full_sentinel_n"] = len(usable)
                summary.loc[index, "gfz_full_sentinel_high_agbd_n"] = int(usable.high_agbd_ge100.sum())
                print(f"{orbit}: original Sentinel-2 footprint quality = {len(usable)}", flush=True)
    summary.to_csv(OUT / "multiyear_spring_gfz_precheck.csv", index=False, encoding="utf-8-sig")
    print(summary[["date_utc", "orbit", "gfz_reported_n", *group_cols]].to_string(index=False))


if __name__ == "__main__":
    main()
