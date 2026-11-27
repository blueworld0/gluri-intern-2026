"""Audit spring CMR granules and screen any available V003 HDF5 shots."""

import argparse
import hashlib
import json
from collections import Counter
from importlib.util import find_spec
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
# h5py lives in week03's Python 3.12 environment in this workspace.
if find_spec("h5py") is None:
    sys.path.append(str(ROOT / "week03" / ".venv" / "Lib" / "site-packages"))

import h5py
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, box, shape

OUT = HERE / "outputs"
sys.path[:0] = [str(ROOT / "week06" / "step2"), str(ROOT / "week08")]
from analyze_expansion import FLAGS, original_quality_functions
from screen_february_forest import classify_points, load_polygons

AOI = (128.75, 35, 129, 35.25)
REFERENCE_ORBIT = "O36071"


def shots_in_area(path, orbit, year):
    physical, quality_mask = original_quality_functions()
    columns = [
        "source_file", "orbit", "beam", "shot_number", "point_id",
        "observed_at_utc", "observation_date_utc", "latitude", "longitude", "agbd_Mg_ha",
    ]
    rows = []
    counts = Counter({"raw_in_aoi": 0})
    spring_start = pd.Timestamp(year=year, month=4, day=1, tz="UTC")
    spring_end = pd.Timestamp(year=year, month=6, day=1, tz="UTC")
    with h5py.File(path) as file:
        if file.attrs["short_name"] != "GEDI_L4A" or not path.name.endswith("_V003.h5"):
            raise ValueError(f"Not the expected GEDI L4A V003 product: {path}")
        for beam in sorted(key for key in file if key.startswith("BEAM")):
            group = file[beam]
            if group["delta_time"].attrs["units"] != "seconds since 2018-01-01":
                raise ValueError(f"Unexpected GEDI time units in {path}")
            x, y = group["lon_lowestmode"][:], group["lat_lowestmode"][:]
            west, south, east, north = AOI
            idx = np.flatnonzero((x >= west) & (x < east) & (y >= south) & (y < north))
            if not len(idx):
                continue
            counts["raw_in_aoi"] += len(idx)
            flags = {name: group[name][idx] for name in FLAGS}
            agbd = physical(group["agbd"][idx], nonnegative=True)
            se = physical(group["agbd_se"][idx], nonnegative=True)
            valid, checks = quality_mask(flags, agbd, se)
            cumulative = np.ones(len(idx), dtype=bool)
            for name, check in checks.items():
                cumulative &= check
                counts[name] += int(cumulative.sum())
            times = pd.to_datetime(group["delta_time"][idx], unit="s", origin="2018-01-01", utc=True)
            for shot, time, lat, lon, density in zip(
                group["shot_number"][idx][valid], times[valid], y[idx][valid], x[idx][valid], agbd[valid]
            ):
                if spring_start <= time < spring_end:
                    rows.append((path.name, orbit, beam, str(shot), f"{path.name}:{beam}:{shot}",
                                 time.isoformat(), time.date().isoformat(), lat, lon, density))
    return pd.DataFrame(rows, columns=columns), dict(counts)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def load_search(extended):
    cmr = read_json(OUT / "gedi_april_may_cmr.json")
    catalog = read_json(OUT / "spring_sentinel_catalog.json")
    if cmr["hits"] != len(cmr["items"]) or catalog["context"]["matched"] != len(catalog["features"]):
        raise ValueError("Incomplete 2025 satellite search")
    granules = cmr["items"].copy()
    images = catalog["features"].copy()
    if extended:
        for path in sorted(OUT.glob("spring_cmr_20??.json")):
            search = read_json(path)
            if search["hits"] != len(search["items"]):
                raise ValueError(f"Incomplete GEDI search: {path}")
            granules.extend(search["items"])
            if not search["items"]:
                continue
            year = path.stem[-4:]
            first = read_json(OUT / f"spring_sentinel_{year}.json")
            year_images = first["features"].copy()
            for page in sorted(OUT.glob(f"spring_sentinel_{year}_page*.json")):
                year_images.extend(read_json(page)["features"])
            if len({item["id"] for item in year_images}) != first["context"]["matched"]:
                raise ValueError(f"Incomplete or duplicate Sentinel-2 search: {year}")
            images.extend(year_images)
    return granules, images


def cmr_area_overlap(granule):
    geom = granule["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]
    area = box(*AOI)
    parts = []
    for item in geom.get("GPolygons", []):
        vertices = np.deg2rad([
            (point["Longitude"], point["Latitude"])
            for point in item["Boundary"]["Points"]
        ])
        lon, lat = vertices.T
        vectors = np.column_stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)))
        steps = []
        for first, last in zip(vectors, np.roll(vectors, -1, axis=0)):
            angle = np.arccos(np.clip(first @ last, -1, 1))
            if angle < 1e-7:
                continue
            fraction = np.linspace(0, 1, max(2, int(np.ceil(np.rad2deg(angle) / 0.05))), endpoint=False)
            arc = (np.sin((1 - fraction) * angle)[:, None] * first
                   + np.sin(fraction * angle)[:, None] * last) / np.sin(angle)
            steps.append(arc)
        x, y, z = np.concatenate(steps).T
        longitudes = np.rad2deg(np.unwrap(np.arctan2(y, x)))
        latitudes = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
        shift = 360 * round((area.centroid.x - np.median(longitudes)) / 360)
        parts.append(Polygon(np.column_stack((longitudes + shift, latitudes))).buffer(0))
    return sum(area.intersection(part).area for part in parts) / area.area * 100


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extended", action="store_true", help="include all saved April–May searches")
    args = parser.parse_args()
    prefix = "multiyear_spring" if args.extended else "spring"
    granules, catalog = load_search(args.extended)
    images = [
        (item, pd.to_datetime(item["properties"]["s2:product_uri"].split("_")[2], utc=True))
        for item in catalog
        if shape(item["geometry"]).intersects(box(*AOI))
    ]
    spatial_index = load_polygons()
    records, point_tables, image_rows, audit = [], [], [], []
    for item in granules:
        granule = item["umm"]
        name = granule["GranuleUR"]
        orbit = re.search(r"_(O\d+)_", name).group(1)
        start = pd.Timestamp(granule["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"])
        end = pd.Timestamp(granule["TemporalExtent"]["RangeDateTime"]["EndingDateTime"])
        ref = orbit == REFERENCE_ORBIT
        raw_dir = ROOT / "week03" / "GEDI04_A" if ref else HERE / "source"
        path = raw_dir / (name + ".h5")
        if not ref and not path.is_file():
            downloads = sorted(raw_dir.rglob(path.name))
            if len(downloads) > 1:
                raise ValueError(f"Multiple downloads found for {path.name}: {downloads}")
            prior_download_dir = ROOT / "week09" / "source"
            prior_path = prior_download_dir / path.name
            if downloads:
                path = downloads[0]
            elif prior_path.is_file():
                path = prior_path
        row = {
            "date_utc": start.date().isoformat(), "orbit": orbit, "granule": name,
            "role": "existing_training_orbit" if ref else "holdout_candidate",
            "source_present": path.is_file(),
            "https_h5": next(link["URL"] for link in granule["RelatedUrls"]
                             if link["URL"].startswith("https://") and link["URL"].endswith(".h5")),
            "gedi_in_aoi_n": None, "gedi_quality_n": None, "forest_n": None, "forest_high_agbd_n": None,
            "sentinel_pixel_valid_n": None,
        }
        if args.extended:
            row["cmr_aoi_overlap_percent_approx"] = round(cmr_area_overlap(granule), 2)
        nearby = [(image, time) for image, time in images
                  if start - pd.Timedelta(days=7) <= time <= end + pd.Timedelta(days=7)]
        row["sentinel_catalogue_products_7days"] = len(nearby)
        if args.extended:
            clouds = [image["properties"]["eo:cloud_cover"] for image, _ in nearby]
            row["lowest_tile_cloud_percent"] = min(clouds, default=None)
        for image, time in nearby:
            image_rows.append({
                "gedi_date_utc": row["date_utc"], "orbit": orbit, "item": image["id"],
                "product": image["properties"]["s2:product_uri"], "sensing_utc": time.isoformat(),
                "hours_from_granule_start": (time - start).total_seconds() / 3600,
                "tile_cloud_percent": image["properties"]["eo:cloud_cover"],
                "scl_url": image["assets"]["scl"]["href"],
            })
        if path.is_file():
            archived = next(entry for entry in granule["DataGranule"]["ArchiveAndDistributionInformation"]
                            if entry["Name"] == path.name)
            with path.open("rb") as handle:
                checksum = hashlib.file_digest(handle, "sha512").hexdigest()
            if path.stat().st_size != archived["SizeInBytes"] or checksum != archived["Checksum"]["Value"].lower():
                raise ValueError(f"File differs from NASA CMR checksum: {path}")
            points, checks = shots_in_area(path, orbit, start.year)
            if not points.point_id.is_unique:
                raise ValueError(f"Duplicate point IDs in {path}")
            if points.empty:
                points["landcover_status"] = pd.Series(dtype="str")
            else:
                points = points.join(classify_points(points, *spatial_index))
            points["high_agbd_ge100"] = points.agbd_Mg_ha.ge(100)
            forest = points.loc[points.landcover_status.eq("forest")]
            row.update(gedi_in_aoi_n=checks["raw_in_aoi"], gedi_quality_n=len(points), forest_n=len(forest),
                       forest_high_agbd_n=int(forest.high_agbd_ge100.sum()))
            previous = OUT / f"{prefix}_{orbit}_summary.json"
            if previous.is_file() and read_json(previous).get("forest_n") == len(forest):
                row["sentinel_pixel_valid_n"] = read_json(previous)["usable_n"]
            audit.append({"orbit": orbit, "sha512": checksum, "cumulative_aoi_quality_counts": checks})
            point_tables.append(points)
        records.append(row)
    inventory = pd.DataFrame(records).sort_values(["date_utc", "orbit"])
    ref = inventory.loc[inventory.orbit.eq(REFERENCE_ORBIT)]
    if ref.gedi_quality_n.tolist() != [1055]:
        raise ValueError("Spring inventory failed to reproduce week06's 1,055 April quality points")
    inventory.to_csv(OUT / f"{prefix}_gedi_inventory.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(image_rows).to_csv(OUT / f"{prefix}_sentinel_candidates.csv", index=False, encoding="utf-8-sig")
    pd.concat(point_tables, ignore_index=True).to_csv(OUT / f"{prefix}_gedi_available_points.csv",
                                                    index=False, encoding="utf-8-sig")
    (OUT / f"{prefix}_gedi_quality_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print(inventory.drop(columns="https_h5").to_string(index=False))


if __name__ == "__main__":
    main()
