"""Offline, chunked inventory of existing GEDI L4A V003 expansion candidates.

Run from any directory using the existing week03 Python environment.
No project modules are imported: exact pure quality functions are extracted via
AST so importing the old plotting module cannot write to Week03 caches.
"""
from pathlib import Path
import ast
import csv
import hashlib
import json
import os
import re
from collections import Counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["MPLCONFIGDIR"] = str(HERE / ".matplotlib")
import h5py
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform

CHUNK = 100_000
CELL = 0.25
# Geographic screening window only; NOT a national/admin/forest boundary.
WINDOW = {"west": 124.0, "south": 33.0, "east": 132.0, "north": 39.0}
FLAGS = ["l2_algrunflag", "l2a_quality_flag_rel3", "l4a_quality_flag_rel3",
         "degrade_include_flag", "elev_highestreturn_outlier_flag"]
SOURCE = ROOT / "week03/analyze_gedi.py"


def original_quality_functions():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef)
             and n.name in {"physical_values", "quality_mask"}]
    assert len(funcs) == 2
    namespace = {"np": np}
    exec(compile(ast.Module(body=funcs, type_ignores=[]), str(SOURCE), "exec"), namespace)
    return namespace["physical_values"], namespace["quality_mask"]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_counter(counter, keys, path):
    pd.DataFrame([dict(zip(keys, k), count=v) for k, v in sorted(counter.items())],
                 columns=keys + ["count"]).to_csv(HERE / path, index=False)


def run():
    physical, quality_mask = original_quality_functions()
    baseline_path = ROOT / "week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv"
    old97_path = ROOT / "week03/outputs/gedi/gedi_valid_points.csv"
    base = pd.read_csv(baseline_path, dtype={"shot_number": "string"})
    old97 = pd.read_csv(old97_path, dtype={"shot_number": "string"})
    baseline_keys = set(zip(base.gedi_file, base.beam, base.shot_number))
    old_keys = set(zip(old97.source_file, old97.beam, old97.shot_number))
    assert len(baseline_keys) == 90 and len(old_keys) == 97
    b04 = sorted((ROOT / "week03/Sentinel-2").glob("*B04*.tiff"))
    assert len(b04) == 1
    with rasterio.open(b04[0]) as ds:
        crs, bounds = ds.crs, ds.bounds
    files = sorted((ROOT / "week03/GEDI04_A").glob("*.h5"))
    assert len(files) == 5
    dimensions = Counter()
    stages = Counter()
    independent = Counter()
    total = 0
    qualified = 0
    quality_invalid_coords = 0
    quality_invalid_time = 0
    point_path = HERE / "quality_valid_points.csv"
    first = True
    manifests = []
    for path in files:
        orbit = re.search(r"_(O\d+)_", path.name).group(1)
        nfile = qfile = 0
        with h5py.File(path, "r") as f:
            assert f.attrs["short_name"] == "GEDI_L4A" and "_V003" in path.name
            for beam in sorted(k for k in f if k.startswith("BEAM")):
                group = f[beam]
                assert group["delta_time"].attrs["units"] == "seconds since 2018-01-01"
                size = len(group["shot_number"])
                for start in range(0, size, CHUNK):
                    sl = slice(start, min(start + CHUNK, size))
                    data = {k: group[k][sl] for k in FLAGS}
                    agbd = physical(group["agbd"][sl], nonnegative=True)
                    se = physical(group["agbd_se"][sl], nonnegative=True)
                    valid, checks = quality_mask(data, agbd, se)
                    times = pd.to_datetime(group["delta_time"][sl], unit="s",
                                           origin="2018-01-01", utc=True, errors="coerce")
                    dates = np.asarray(times.strftime("%Y-%m-%d").fillna("INVALID_TIME"))
                    cumulative = np.ones(len(valid), dtype=bool)
                    stages["all_observations"] += len(valid)
                    for name, check in checks.items():
                        cumulative &= check
                        stages[name] += int(cumulative.sum())
                        independent[name] += int(check.sum())
                    for day in np.unique(dates):
                        sel = dates == day
                        dimensions[(path.name, day, orbit, beam, "all")] += int(sel.sum())
                        dimensions[(path.name, day, orbit, beam, "quality_valid")] += int((sel & valid).sum())
                    total += len(valid)
                    nfile += len(valid)
                    qualified += int(valid.sum())
                    qfile += int(valid.sum())
                    if not valid.any():
                        continue
                    lat = group["lat_lowestmode"][sl][valid]
                    lon = group["lon_lowestmode"][sl][valid]
                    coord_ok = np.isfinite(lat) & np.isfinite(lon) & (abs(lat) <= 90) & (abs(lon) <= 180)
                    quality_invalid_coords += int((~coord_ok).sum())
                    quality_invalid_time += int(times[valid].isna().sum())
                    shots = group["shot_number"][sl][valid].astype(str)
                    ids = [(path.name, beam, s) for s in shots]
                    inside = np.zeros(len(lat), dtype=bool)
                    if coord_ok.any():
                        x, y = map(np.asarray, transform("EPSG:4326", crs, lon[coord_ok].tolist(), lat[coord_ok].tolist()))
                        inside[coord_ok] = ((x >= bounds.left) & (x <= bounds.right)
                                            & (y >= bounds.bottom) & (y <= bounds.top))
                    frame = pd.DataFrame({
                        "source_file": path.name, "orbit": orbit, "beam": beam,
                        "shot_number": shots,
                        "point_id": [":".join(k) for k in ids],
                        "observed_at_utc": times[valid].strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
                        "observation_date_utc": dates[valid], "latitude": lat, "longitude": lon,
                        "agbd_Mg_ha": agbd[valid], "agbd_se_Mg_ha": se[valid],
                        **{k: v[valid] for k, v in data.items()},
                        "degrade_flag": group["degrade_flag"][sl][valid],
                        "sensitivity": physical(group["sensitivity"][sl])[valid] if "sensitivity" in group else np.nan,
                        "coordinate_valid": coord_ok, "quality_valid": True,
                        "existing_week04_05_90": [k in baseline_keys for k in ids],
                        "existing_gedi_97": [k in old_keys for k in ids],
                        "within_original_sentinel_rectangle": inside,
                        "sentinel2_status": "not_evaluated_in_step2",
                        "forest_status": "not_evaluated_in_step2",
                    })
                    frame.to_csv(point_path, mode="w" if first else "a", header=first, index=False)
                    first = False
        manifests.append({"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                          "all_observations": nfile, "quality_valid": qfile})
        print(f"{path.name}: all={nfile:,}; quality={qfile:,}", flush=True)
    write_counter(dimensions, ["source_file", "observation_date_utc", "orbit", "beam", "stage"], "counts_granule_date_orbit_beam.csv")
    detail = pd.read_csv(HERE / "counts_granule_date_orbit_beam.csv")
    for fields, name in [(["source_file", "orbit"], "granule"), (["observation_date_utc"], "date"),
                         (["orbit"], "orbit"), (["beam"], "beam")]:
        table = detail.groupby(fields + ["stage"], dropna=False)["count"].sum().unstack("stage", fill_value=0).reset_index()
        table["quality_pass_percent"] = table.quality_valid / table["all"] * 100
        table.to_csv(HERE / f"counts_by_{name}.csv", index=False)
    pd.DataFrame([{"stage": k, "cumulative_count": v,
                   "independent_pass_count": independent.get(k, v)} for k, v in stages.items()]).to_csv(HERE / "quality_filter_counts.csv", index=False)
    # Only quality-valid rows are materialized, never the 6.35 million input rows.
    p = pd.read_csv(point_path, dtype={"shot_number": "string"})
    assert len(p) == qualified and not p.point_id.duplicated().any()
    assert int(p.existing_week04_05_90.sum()) == 90
    assert int(p.existing_gedi_97.sum()) == 97
    assert int(p.within_original_sentinel_rectangle.sum()) == 97
    assert total == 6_356_011
    geo = p.loc[p.coordinate_valid].copy()
    geo["cell_lon_index"] = np.floor(geo.longitude / CELL).astype(int)
    geo["cell_lat_index"] = np.floor(geo.latitude / CELL).astype(int)
    geo["region_id"] = [f"G025_{x}_{y}" for x, y in zip(geo.cell_lon_index, geo.cell_lat_index)]
    geo["in_korea_screening_window"] = (geo.longitude.ge(WINDOW["west"]) & geo.longitude.lt(WINDOW["east"])
        & geo.latitude.ge(WINDOW["south"]) & geo.latitude.lt(WINDOW["north"]))
    geo[["point_id", "region_id", "in_korea_screening_window"]].to_csv(HERE / "point_region_assignment.csv", index=False)
    # Connected occupied geographic cells, 8-neighbor topology; no distances.
    occupied = set(zip(geo.cell_lon_index, geo.cell_lat_index))
    components = {}
    component = 0
    for cell in sorted(occupied):
        if cell in components:
            continue
        component += 1
        stack = [cell]
        components[cell] = component
        while stack:
            x, y = stack.pop()
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    neighbor = (x + dx, y + dy)
                    if neighbor in occupied and neighbor not in components:
                        components[neighbor] = component
                        stack.append(neighbor)
    geo["occupied_grid_component"] = [components[c] for c in zip(geo.cell_lon_index, geo.cell_lat_index)]
    geo[["point_id", "occupied_grid_component"]].to_csv(HERE / "point_component_assignment.csv", index=False)
    def summarize(group):
        return pd.Series({"quality_valid": len(group),
            "existing_90": int(group.existing_week04_05_90.sum()),
            "additional_vs_90": int((~group.existing_week04_05_90).sum()),
            "existing_97": int(group.existing_gedi_97.sum()),
            "outside_original_rectangle": int((~group.within_original_sentinel_rectangle).sum()),
            "latitude_min": group.latitude.min(), "latitude_max": group.latitude.max(),
            "longitude_min": group.longitude.min(), "longitude_max": group.longitude.max(),
            "date_count": group.observation_date_utc.nunique(), "orbit_count": group.orbit.nunique(),
            "beam_count": group.beam.nunique(), "dates_utc": ";".join(sorted(group.observation_date_utc.unique())),
            "orbits": ";".join(sorted(group.orbit.unique())), "agbd_median_Mg_ha": group.agbd_Mg_ha.median()})
    def grouped_summary(data, key):
        return pd.DataFrame([dict({key: k}, **summarize(g).to_dict()) for k, g in data.groupby(key)])
    grids = grouped_summary(geo, "region_id")
    indexes = geo.groupby("region_id")[["cell_lon_index", "cell_lat_index"]].first().reset_index()
    grids = grids.merge(indexes, on="region_id", validate="one_to_one")
    grids["grid_west"] = grids.cell_lon_index * CELL
    grids["grid_south"] = grids.cell_lat_index * CELL
    grids["grid_east"] = grids.grid_west + CELL
    grids["grid_north"] = grids.grid_south + CELL
    grids = grids.sort_values("quality_valid", ascending=False)
    grids.to_csv(HERE / "regions_global_025deg.csv", index=False)
    regional = geo.loc[geo.in_korea_screening_window]
    region_ids = set(regional.region_id)
    kg = grids.loc[grids.region_id.isin(region_ids)].copy()
    kg.to_csv(HERE / "regions_korea_window_025deg.csv", index=False)
    grouped_summary(geo, "occupied_grid_component").sort_values("quality_valid", ascending=False).to_csv(HERE / "occupied_grid_components.csv", index=False)
    geo.groupby(["region_id", "observation_date_utc", "orbit", "beam"]).agg(
        quality_valid=("point_id", "size"), existing_90=("existing_week04_05_90", "sum"),
        existing_97=("existing_gedi_97", "sum")).reset_index().to_csv(HERE / "counts_region_date_orbit_beam.csv", index=False)
    summary = {
        "all_observations": total, "quality_valid": qualified,
        "quality_pass_percent": qualified / total * 100,
        "quality_valid_invalid_coordinates": quality_invalid_coords,
        "quality_valid_invalid_time": quality_invalid_time,
        "existing_week04_05_90": int(p.existing_week04_05_90.sum()),
        "existing_gedi_97": int(p.existing_gedi_97.sum()),
        "additional_quality_valid_vs_90": int((~p.existing_week04_05_90).sum()),
        "outside_original_rectangle": int((~p.within_original_sentinel_rectangle).sum()),
        "inside_original_rectangle_not_in_90": int((p.within_original_sentinel_rectangle & ~p.existing_week04_05_90).sum()),
        "geolocatable_quality_valid": len(geo), "bounds_wgs84": {
            "west": float(geo.longitude.min()), "south": float(geo.latitude.min()),
            "east": float(geo.longitude.max()), "north": float(geo.latitude.max())},
        "korea_screening_window": WINDOW, "korea_window_counts": summarize(regional).to_dict(),
        "occupied_025_degree_cells": len(grids), "occupied_grid_components": component,
        "final_training_ready_new_points": None,
        "final_training_ready_status": "Not evaluated: no new Sentinel-2 or forest mask checks",
        "chunk_size": CHUNK, "grid_cell_degrees": CELL,
        "source_quality_script_sha256": sha256(SOURCE),
        "baseline_90_sha256": sha256(baseline_path), "old_97_sha256": sha256(old97_path),
        "inputs": manifests,
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
    print("Korea screening cells:", flush=True)
    print(kg.head(25).to_string(index=False), flush=True)
    # Bounded deterministic plotting samples; all regional points shown if small.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    for ax, data, title, limit in [(axes[0], geo, "Global quality-valid GEDI shots", 160_000),
                                   (axes[1], regional, "Korea geographic screening window (not a boundary)", 100_000)]:
        display = data if len(data) <= limit else data.sample(limit, random_state=42)
        for orbit, group in display.groupby("orbit"):
            ax.scatter(group.longitude, group.latitude, s=1, alpha=.5, label=orbit, rasterized=True)
        if ax == axes[1]:
            ax.scatter(base.longitude, base.latitude, s=20, c="black", marker="x", label="Existing 90")
            ax.set_xlim(WINDOW["west"], WINDOW["east"])
            ax.set_ylim(WINDOW["south"], WINDOW["north"])
        ax.set(xlabel="Longitude (WGS84 degrees)", ylabel="Latitude (WGS84 degrees)",
               title=f"{title}\nn={len(data):,}; displayed={len(display):,}; no basemap")
        ax.grid(alpha=.25)
        ax.legend(markerscale=4, fontsize=8)
    fig.savefig(HERE / "spatial_distribution.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    run()
