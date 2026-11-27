"""Screen available spring GEDI orbits and compare fixed models orbit by orbit."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
from urllib.parse import urlsplit

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

from compare_topography import APRIL, sample_dsm
from evaluate_february_holdout import footprint_features, scores
from inventory_spring import HERE, OUT, load_search
from screen_february_sentinel import screen

AWS = "https://sentinel-s2-l2a.s3.eu-central-1.amazonaws.com/"
SOURCE = HERE / "source"
TOPOGRAPHY = ["dsm_elevation_m", "dsm_slope_deg"]


def download(url, path):
    if path.is_file() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = path.with_suffix(path.suffix + ".part")
    subprocess.run(["curl.exe", "-fsSL", "--retry", "2", "--max-time", "200",
                    "-o", str(staging), url], check=True)
    staging.replace(path)


def source_files(item, bands):
    product = item["properties"]["s2:product_uri"]
    source = SOURCE / "spring_tiles" / item["id"]
    jp2 = urlsplit(item["assets"]["scl-jp2"]["href"])
    if jp2.scheme != "s3" or jp2.netloc != "sentinel-s2-l2a" or not jp2.path.endswith("/SCL.jp2"):
        raise ValueError(f"Unexpected original tile asset for {item['id']}")
    tile_prefix = jp2.path.lstrip("/").removesuffix("SCL.jp2").removesuffix("R20m/")
    base = AWS + tile_prefix
    download(base + "productInfo.json", source / "productInfo.json")
    info = json.loads((source / "productInfo.json").read_text(encoding="utf-8"))
    if info["name"] + ".SAFE" != product:
        raise ValueError(f"AWS tile does not match the selected product: {product}")
    tasks = [
        (base + "R20m/SCL.jp2", source / "SCL.jp2"),
        (base + "qi/CLASSI_B00.jp2", source / "CLASSI_B00.jp2"),
        (AWS + info["path"] + "/metadata.xml", source / "product_metadata.xml"),
    ]
    for band in bands:
        resolution = "R10m" if band in ("B02", "B03", "B04", "B08") else "R20m"
        for name in (f"{resolution}/{band}.jp2", f"qi/QUALIT_{band}.jp2", f"qi/DETFOO_{band}.jp2"):
            tasks.append((base + name, source / Path(name).name))
    with ThreadPoolExecutor(max_workers=5) as pool:
        list(pool.map(lambda task: download(*task), tasks))
    return source


def choose_product(points, items):
    scl_path = SOURCE / "spring_scl"
    with ThreadPoolExecutor(max_workers=5) as pool:
        list(pool.map(lambda item: download(item["assets"]["scl"]["href"], scl_path / f"{item['id']}.tif"),
                      items))
    screened = pd.concat([screen(item, points, scl_path) for item in items], ignore_index=True)
    screened["scl_vegetation"] &= screened.offset_hours.abs().le(168)
    screened["abs_offset_hours"] = screened.offset_hours.abs()
    screened["high_agbd_scl_pass"] = screened.scl_vegetation & screened.high_agbd_ge100
    counts = screened.groupby("item_id", as_index=False).agg(
        scl_n=("scl_vegetation", "sum"), high_agbd_scl_n=("high_agbd_scl_pass", "sum"),
        hours=("abs_offset_hours", "median"),
    ).sort_values(["scl_n", "hours", "item_id"], ascending=[False, True, True])
    best = counts.iloc[0].item_id
    return best, counts, screened.loc[screened.item_id.eq(best) & screened.scl_vegetation].copy()


def compatible_products(items, ids):
    """Select Collection-1 processing baselines with the same mask family as week08."""
    return [
        item for item in items
        if item["id"] in ids
        and float(item["properties"]["s2:processing_baseline"]) >= 5.0
    ]


def evaluate(usable, train, config):
    features = config["feature_sets"]["bands10"]
    coords = ["point_id", "easting_m", "northing_m"]
    if usable.point_id.isin(train.point_id).any():
        raise ValueError("Training and validation IDs overlap")
    dsm = sample_dsm(pd.concat([train[coords], usable[coords]], ignore_index=True))
    train = train.merge(dsm, on="point_id", validate="one_to_one")
    usable = usable.merge(dsm, on="point_id", validate="one_to_one")
    usable["dummy_mean"] = train.agbd.mean()
    for name, columns in (("bands10", features), ("bands10_with_dsm", features + TOPOGRAPHY)):
        model = RandomForestRegressor(**config["rf_params"])
        model.fit(train[columns], train.agbd)
        usable[name] = model.predict(usable[columns])
    distance = np.hypot(
        usable.easting_m.to_numpy()[:, None] - train.easting_m.to_numpy(),
        usable.northing_m.to_numpy()[:, None] - train.northing_m.to_numpy(),
    ).min(axis=1)
    usable["nearest_april_training_m"] = distance
    metrics = {}
    for scope, frame in (("all", usable), ("agbd_ge100", usable.loc[usable.high_agbd_ge100])):
        if not frame.empty:
            metrics[scope] = {
                name: scores(frame.agbd_Mg_ha, frame[name])
                for name in ("dummy_mean", "bands10", "bands10_with_dsm")
            }
    return usable, metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbits", nargs="+")
    parser.add_argument("--extended", action="store_true", help="use the multi-year inventory")
    args = parser.parse_args()
    prefix = "multiyear_spring" if args.extended else "spring"
    inventory = pd.read_csv(OUT / f"{prefix}_gedi_inventory.csv").set_index("orbit")
    requested = args.orbits if args.orbits else inventory.index.tolist()
    orbit_ids = [orbit for orbit in requested if orbit in inventory.index
                 and inventory.loc[orbit, "role"] == "holdout_candidate"
                 and inventory.loc[orbit, "source_present"]]
    if not orbit_ids:
        parser.exit(0, "No spring holdout files available. Run inventory_spring.py after download.\n")
    points = pd.read_csv(OUT / f"{prefix}_gedi_available_points.csv")
    candidates = pd.read_csv(OUT / f"{prefix}_sentinel_candidates.csv")
    _, catalog = load_search(args.extended)
    config = json.loads((APRIL / "rededge_config.json").read_text(encoding="utf-8"))
    train = pd.read_csv(APRIL / "rededge_forest184_input.csv")
    bands = [name.removeprefix("reflectance_") for name in config["feature_sets"]["bands10"]]
    for orbit in orbit_ids:
        points_for_orbit = points.loc[points.orbit.eq(orbit) & points.landcover_status.eq("forest")].copy()
        run = {
            "orbit": orbit, "forest_n": len(points_for_orbit), "usable_n": 0,
            "status": "no_gedi_quality_pass" if inventory.loc[orbit, "gedi_quality_n"] == 0 else "no_forest_points",
            "metrics": {},
        }
        candidate_items = set(candidates.loc[candidates.orbit.eq(orbit), "item"])
        items = compatible_products(catalog, candidate_items)
        if not points_for_orbit.empty and items:
            best, counts, selected = choose_product(points_for_orbit, items)
            counts.to_csv(OUT / f"{prefix}_{orbit}_scl_counts.csv", index=False)
            run.update({"item_id": best, "scl_n": len(selected),
                        "scl_high_agbd_n": int(selected.high_agbd_ge100.sum()),
                        "status": "no_scl_pass"})
            if not selected.empty:
                item = next(item for item in items if item["id"] == best)
                product = item["properties"]["s2:product_uri"]
                source = source_files(item, bands)
                data = footprint_features(selected, source, bands, product)
                if not data.original_scl_aligned_center.eq(data.scl_aligned_10m).all():
                    raise ValueError("Public COG and original SCL disagree at a selected point")
                data.to_csv(OUT / f"{prefix}_{orbit}_footprint_audit.csv", index=False, encoding="utf-8-sig")
                usable = data.loc[data.footprint_pass].copy()
                run.update({"product": product, "usable_n": len(usable),
                            "usable_high_agbd_n": int(usable.high_agbd_ge100.sum()),
                            "status": "insufficient_usable_points"})
                if len(usable) >= 2:
                    predictions, run["metrics"] = evaluate(usable, train, config)
                    predictions.to_csv(OUT / f"{prefix}_{orbit}_predictions.csv", index=False, encoding="utf-8-sig")
                    run["minimum_train_test_distance_m"] = float(predictions.nearest_april_training_m.min())
                    spatial_subset = predictions.loc[predictions.nearest_april_training_m.ge(500)]
                    run["n_within_500m_of_training"] = len(predictions) - len(spatial_subset)
                    if not spatial_subset.empty:
                        run["metrics_at_least_500m_from_training"] = {
                            name: scores(spatial_subset.agbd_Mg_ha, spatial_subset[name])
                            for name in ("dummy_mean", "bands10", "bands10_with_dsm")
                        }
                    run["status"] = "evaluated"
        elif not points_for_orbit.empty:
            run["status"] = "no_sentinel_products"
        payload = json.dumps(run, indent=2, allow_nan=False)
        (OUT / f"{prefix}_{orbit}_summary.json").write_text(payload + "\n", encoding="utf-8")
        inventory.loc[orbit, "sentinel_pixel_valid_n"] = run["usable_n"]
        inventory.sort_values("date_utc").reset_index().to_csv(
            OUT / f"{prefix}_gedi_inventory.csv", index=False, encoding="utf-8-sig"
        )
        print(payload, flush=True)


if __name__ == "__main__":
    main()
