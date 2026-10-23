"""Inspect GEDI04_A V003 shots and analyze those inside the Sentinel-2 extent."""

import argparse
from collections import Counter
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform, transform_bounds


ROOT = Path(__file__).resolve().parent
GEDI_DIR = ROOT / "GEDI04_A"
S2_DIR = ROOT / "Sentinel-2"
DEFAULT_OUTPUT = ROOT / "outputs" / "gedi"
TIME_ORIGIN = pd.Timestamp("2018-01-01")
REQUIRED = (
    "shot_number", "delta_time", "lat_lowestmode", "lon_lowestmode", "agbd", "agbd_se",
    "l2_algrunflag", "l2a_quality_flag_rel3", "l4a_quality_flag_rel3",
    "degrade_flag", "degrade_include_flag", "elev_highestreturn_outlier_flag",
)
OPTIONAL = ("agbd_pi_lower", "agbd_pi_upper", "predictor_limit_flag", "response_limit_flag", "sensitivity")
COLUMNS = [
    "source_file", "beam", "shot_number", "observed_at_utc", "latitude", "longitude",
    "agbd_Mg_ha", "agbd_se_Mg_ha", "agbd_pi_lower_Mg_ha", "agbd_pi_upper_Mg_ha",
    "l2_algrunflag", "l2a_quality_flag_rel3", "l4a_quality_flag_rel3", "degrade_flag",
    "degrade_include_flag", "elev_highestreturn_outlier_flag", "predictor_limit_flag",
    "response_limit_flag", "sensitivity", "within_sentinel_extent", "quality_valid",
]
REGION_COLUMNS = COLUMNS + ["easting_m", "northing_m"]
STATS_COLUMNS = ["metric", "value", "unit"]
COUNT_COLUMNS = ["source_file", "observation_date_utc", "all_points", "within_sentinel_extent", "quality_valid"]


def find_reference() -> tuple[Path, dict]:
    candidates = sorted(p for p in S2_DIR.rglob("*") if p.suffix.lower() in {".tif", ".tiff"} and "B04" in p.name.upper())
    if len(candidates) != 1:
        raise ValueError(f"Expected one identifiable Sentinel-2 B04 GeoTIFF, found {len(candidates)}")
    with rasterio.open(candidates[0]) as ds:
        if ds.crs is None:
            raise ValueError(f"Reference GeoTIFF has no CRS: {candidates[0]}")
        reference = {"crs": ds.crs, "bounds": ds.bounds, "transform": ds.transform,
                     "width": ds.width, "height": ds.height}
    return candidates[0], reference


def read_dataset(group: h5py.Group, name: str, count: int) -> np.ndarray | None:
    if name not in group:
        return None
    ds = group[name]
    if not isinstance(ds, h5py.Dataset) or ds.shape != (count,):
        raise ValueError(f"Unexpected shape at {ds.name}: {ds.shape}, expected {(count,)}")
    return ds[:]


def physical_values(values: np.ndarray | None, *, nonnegative: bool = False) -> np.ndarray | None:
    if values is None:
        return None
    result = values.astype(np.float64)
    result[~np.isfinite(result)] = np.nan
    # AGBD, standard error, and prediction interval bounds cannot use the observed -9999 sentinel.
    result[result == -9999] = np.nan
    if nonnegative:
        result[result < 0] = np.nan
    return result


def quality_mask(data: dict[str, np.ndarray], agbd: np.ndarray, se: np.ndarray) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    checks = {
        "finite_nonnegative_agbd": np.isfinite(agbd),
        "finite_nonnegative_agbd_se": np.isfinite(se),
        "algorithm_run": data["l2_algrunflag"] == 1,
        "l2a_quality_rel3": data["l2a_quality_flag_rel3"] == 1,
        "l4a_quality_rel3": data["l4a_quality_flag_rel3"] == 1,
        "degrade_include": data["degrade_include_flag"] == 1,
        "no_elevation_outlier": data["elev_highestreturn_outlier_flag"] == 0,
    }
    valid = np.ones(len(agbd), dtype=bool)
    for name, check in checks.items():
        valid &= check
    return valid, checks


def _background(ax, reference: dict) -> None:
    bands = {}
    for band in ("B04", "B03", "B02"):
        matches = [p for p in S2_DIR.iterdir() if p.suffix.lower() in {".tif", ".tiff"} and band in p.name.upper()]
        if len(matches) != 1:
            return
        with rasterio.open(matches[0]) as ds:
            if ds.crs != reference["crs"] or ds.transform != reference["transform"]:
                return
            array = ds.read(1).astype(np.float32)
            mask = ds.read_masks(1) > 0
            if ds.count == 2:
                mask &= ds.read(2) > 0
            array[~mask] = np.nan
            bands[band] = array
    raw = np.stack([bands[b] for b in ("B04", "B03", "B02")], axis=-1)
    valid = np.all(np.isfinite(raw), axis=-1)
    if not valid.any():
        return
    rgb = np.zeros_like(raw, dtype=np.float32)
    for channel in range(3):
        lo, hi = np.percentile(raw[:, :, channel][valid], [2, 98])
        if hi > lo:
            rgb[:, :, channel] = np.clip((raw[:, :, channel] - lo) / (hi - lo), 0, 1)
    rgb[~valid] = 0
    bounds = reference["bounds"]
    ax.imshow(rgb, extent=(bounds.left, bounds.right, bounds.bottom, bounds.top), origin="upper")


def plot_map(points: pd.DataFrame, reference: dict, path: Path, title: str,
             colored: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    _background(ax, reference)
    bounds = reference["bounds"]
    ax.set_xlim(bounds.left, bounds.right)
    ax.set_ylim(bounds.bottom, bounds.top)
    if len(points):
        if colored:
            colors = points["agbd_Mg_ha"].to_numpy()
            lo, hi = np.percentile(colors, [2, 98])
            if lo == hi:
                lo, hi = lo - 1, hi + 1
            artist = ax.scatter(points["easting_m"], points["northing_m"], c=colors,
                                s=20, cmap="YlOrRd", vmin=lo, vmax=hi,
                                edgecolors="black", linewidths=0.15, alpha=0.9)
            fig.colorbar(artist, ax=ax, label="GEDI predicted AGBD (Mg/ha)")
        else:
            ax.scatter(points["easting_m"], points["northing_m"], s=10,
                       color="cyan", edgecolors="black", linewidths=0.15, alpha=0.75)
    else:
        ax.text(0.5, 0.5, "No GEDI shots in this extent\nReview dates or analysis bounds",
                ha="center", va="center", transform=ax.transAxes,
                bbox={"facecolor": "white", "alpha": 0.9})
    ax.set_title(f"{title}\nSentinel-2 rectangle: part of Geumjeongsan | n={len(points):,}", fontsize=11)
    ax.set_xlabel(f"Easting (m, {reference['crs']})")
    ax.set_ylabel("Northing (m)")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_histogram(points: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    if len(points):
        ax.hist(points["agbd_Mg_ha"].dropna(), bins=30, color="forestgreen")
    else:
        ax.text(0.5, 0.5, "No quality-valid GEDI shots in this extent",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_title(f"GEDI04_A V003 predicted AGBD | valid shots n={len(points):,}")
    ax.set_xlabel("AGBD (Mg/ha)")
    ax.set_ylabel("Shot count")
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def statistics(all_count: int, region: pd.DataFrame, valid: pd.DataFrame,
               files: int) -> pd.DataFrame:
    rows = [("files", files, "files"), ("all_shots", all_count, "shots"),
            ("within_sentinel_extent", len(region), "shots"),
            ("quality_valid", len(valid), "shots")]
    for column, label in (("agbd_Mg_ha", "agbd"), ("agbd_se_Mg_ha", "agbd_se"),
                          ("agbd_pi_lower_Mg_ha", "agbd_pi_lower"),
                          ("agbd_pi_upper_Mg_ha", "agbd_pi_upper")):
        values = valid[column].dropna().to_numpy(dtype=float) if column in valid else np.array([])
        rows.append((f"{label}_count", len(values), "shots"))
        for key, func in (("min", np.min), ("max", np.max), ("mean", np.mean),
                          ("median", np.median), ("std", np.std)):
            rows.append((f"{label}_{key}", float(func(values)) if len(values) else np.nan, "Mg/ha"))
    return pd.DataFrame(rows, columns=STATS_COLUMNS)


def analyze(output: Path) -> None:
    files = sorted(GEDI_DIR.glob("*.h5")) + sorted(GEDI_DIR.glob("*.hdf5"))
    if not files:
        raise FileNotFoundError(f"No GEDI HDF5 files in {GEDI_DIR}")
    reference_path, reference = find_reference()
    geo_bounds = transform_bounds(reference["crs"], "EPSG:4326", *reference["bounds"], densify_pts=21)
    print(f"Sentinel reference: {reference_path.name}, {reference['crs']}, bounds={reference['bounds']}")
    print(f"Approximate WGS84 bounding box: {geo_bounds}")
    print(f"GEDI files: {len(files)}")
    targets = ["gedi_all_points.csv", "gedi_geumjeongsan_points.csv", "gedi_valid_points.csv",
               "gedi_statistics.csv", "gedi_counts_by_file_date.csv", "gedi_quality_counts.csv"]
    figures = ["gedi_footprints.png", "gedi_valid_footprints.png", "gedi_agbd_map.png", "gedi_agbd_histogram.png"]
    existing = [p for p in [*(output / name for name in targets), *(output / "figures" / name for name in figures)] if p.exists()]
    if existing:
        raise FileExistsError(f"GEDI output already exists; refusing to overwrite: {existing[0]}")
    (output / "figures").mkdir(parents=True, exist_ok=True)
    all_path = output / "gedi_all_points.csv"
    first_chunk = True
    all_count = 0
    region_parts: list[pd.DataFrame] = []
    file_date_counts: Counter = Counter()
    quality_counts: Counter = Counter()
    quality_order = ("finite_nonnegative_agbd", "finite_nonnegative_agbd_se", "algorithm_run",
                     "l2a_quality_rel3", "l4a_quality_rel3", "degrade_include", "no_elevation_outlier")
    for path in files:
        file_all = file_region = file_valid = 0
        with h5py.File(path, "r") as h5:
            if "_V003" not in path.name or h5.attrs.get("short_name") != "GEDI_L4A":
                raise ValueError(f"Not a recognized GEDI04_A V003 granule: {path.name}")
            beams = sorted(key for key in h5.keys() if key.startswith("BEAM") and isinstance(h5[key], h5py.Group))
            print(f"\n{path.name}: BEAM groups={beams}")
            if beams:
                example = h5[beams[0]]
                for name in ("lat_lowestmode", "lon_lowestmode", "agbd", "agbd_se",
                             "l4a_quality_flag_rel3", "degrade_flag", "degrade_include_flag"):
                    ds = example[name]
                    print(f"  {beams[0]}/{name}: shape={ds.shape}, dtype={ds.dtype}, "
                          f"HDF5 fill={ds.fillvalue}, attrs={dict(ds.attrs)}")
            for beam in beams:
                group = h5[beam]
                missing = [name for name in REQUIRED if name not in group]
                if missing:
                    raise ValueError(f"{path.name}/{beam}: missing V003 datasets: {missing}")
                count = len(group["shot_number"])
                data = {name: read_dataset(group, name, count) for name in REQUIRED + OPTIONAL}
                if group["delta_time"].attrs.get("units") != "seconds since 2018-01-01":
                    raise ValueError(f"Unexpected delta_time units in {path.name}/{beam}")
                lat, lon = data["lat_lowestmode"], data["lon_lowestmode"]
                coord_valid = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
                west, south, east, north = geo_bounds
                candidates = np.flatnonzero(coord_valid & (lon >= west) & (lon <= east) & (lat >= south) & (lat <= north))
                inside = np.zeros(count, dtype=bool)
                x = y = np.array([], dtype=float)
                if len(candidates):
                    x, y = (np.asarray(z) for z in transform("EPSG:4326", reference["crs"],
                                                            lon[candidates].tolist(), lat[candidates].tolist()))
                    bounds = reference["bounds"]
                    exact = (x >= bounds.left) & (x <= bounds.right) & (y >= bounds.bottom) & (y <= bounds.top)
                    inside[candidates[exact]] = True
                agbd = physical_values(data["agbd"], nonnegative=True)
                se = physical_values(data["agbd_se"], nonnegative=True)
                quality, checks = quality_mask(data, agbd, se)
                dates = pd.to_datetime(data["delta_time"], unit="s", origin=TIME_ORIGIN,
                                       utc=True, errors="coerce")
                date_str = dates.strftime("%Y-%m-%d")
                observed = dates.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                fields = {
                    "source_file": path.name, "beam": beam, "shot_number": data["shot_number"],
                    "observed_at_utc": observed, "latitude": lat, "longitude": lon,
                    "agbd_Mg_ha": agbd, "agbd_se_Mg_ha": se,
                    "agbd_pi_lower_Mg_ha": physical_values(data["agbd_pi_lower"]),
                    "agbd_pi_upper_Mg_ha": physical_values(data["agbd_pi_upper"]),
                    "l2_algrunflag": data["l2_algrunflag"],
                    "l2a_quality_flag_rel3": data["l2a_quality_flag_rel3"],
                    "l4a_quality_flag_rel3": data["l4a_quality_flag_rel3"],
                    "degrade_flag": data["degrade_flag"],
                    "degrade_include_flag": data["degrade_include_flag"],
                    "elev_highestreturn_outlier_flag": data["elev_highestreturn_outlier_flag"],
                    "predictor_limit_flag": data["predictor_limit_flag"],
                    "response_limit_flag": data["response_limit_flag"],
                    "sensitivity": physical_values(data["sensitivity"]),
                    "within_sentinel_extent": inside, "quality_valid": quality,
                }
                frame = pd.DataFrame(fields, columns=COLUMNS)
                frame.to_csv(all_path, mode="w" if first_chunk else "a", header=first_chunk,
                             index=False, na_rep="", float_format="%.8f")
                first_chunk = False
                selected = frame.loc[inside].copy()
                if len(selected):
                    selected["easting_m"] = np.nan
                    selected["northing_m"] = np.nan
                    lookup = {int(i): (float(px), float(py)) for i, px, py in zip(candidates, x, y)}
                    selected["easting_m"] = [lookup[int(i)][0] for i in selected.index]
                    selected["northing_m"] = [lookup[int(i)][1] for i in selected.index]
                    region_parts.append(selected[REGION_COLUMNS])
                for day in np.unique(date_str):
                    day_mask = date_str == day
                    file_date_counts[(path.name, day, "all")] += int(day_mask.sum())
                    file_date_counts[(path.name, day, "region")] += int((day_mask & inside).sum())
                    file_date_counts[(path.name, day, "valid")] += int((day_mask & inside & quality).sum())
                quality_counts["region_total"] += int(inside.sum())
                cumulative = inside.copy()
                for name in quality_order:
                    cumulative &= checks[name]
                    quality_counts[name] += int(cumulative.sum())
                all_count += count
                file_all += count
                file_region += int(inside.sum())
                file_valid += int((inside & quality).sum())
                print(f"  {beam}: {count:,} shots, region={int(inside.sum()):,}, valid={int((inside & quality).sum()):,}")
        print(f"  File totals: all={file_all:,}, region={file_region:,}, valid={file_valid:,}")
    region = pd.concat(region_parts, ignore_index=True) if region_parts else pd.DataFrame(columns=REGION_COLUMNS)
    valid = region.loc[region["quality_valid"] == True].copy() if len(region) else region.copy()
    region.to_csv(output / "gedi_geumjeongsan_points.csv", index=False)
    valid.to_csv(output / "gedi_valid_points.csv", index=False)
    statistics(all_count, region, valid, len(files)).to_csv(output / "gedi_statistics.csv", index=False)
    count_rows = []
    for file_name, day, kind in sorted(file_date_counts):
        if kind != "all":
            continue
        count_rows.append((file_name, day, file_date_counts[(file_name, day, "all")],
                           file_date_counts[(file_name, day, "region")],
                           file_date_counts[(file_name, day, "valid")]))
    pd.DataFrame(count_rows, columns=COUNT_COLUMNS).to_csv(output / "gedi_counts_by_file_date.csv", index=False)
    pd.DataFrame([{"stage": "within_sentinel_extent", "count": len(region)}] +
                 [{"stage": name, "count": quality_counts[name]} for name in quality_order] +
                 [{"stage": "quality_valid", "count": len(valid)}]).to_csv(output / "gedi_quality_counts.csv", index=False)
    fig_dir = output / "figures"
    plot_map(region, reference, fig_dir / "gedi_footprints.png", "All GEDI shots in analysis extent")
    plot_map(valid, reference, fig_dir / "gedi_valid_footprints.png", "Quality-valid GEDI shots")
    plot_map(valid, reference, fig_dir / "gedi_agbd_map.png", "Quality-valid GEDI AGBD", colored=True)
    plot_histogram(valid, fig_dir / "gedi_agbd_histogram.png")
    print(f"\nTOTAL: files={len(files)}, all={all_count:,}, in Sentinel rectangle={len(region):,}, quality-valid={len(valid):,}")
    print(statistics(all_count, region, valid, len(files)).to_string(index=False))
    if not len(region):
        print("No GEDI shots inside the Sentinel-2 extent. Review date coverage or analysis bounds.")
    elif not len(valid):
        print("No quality-valid GEDI shots. Review dates and quality flags before data fusion.")
    print(f"Output: {output}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT,
                        help="New GEDI output directory; existing results are never overwritten")
    args = parser.parse_args()
    try:
        analyze(args.output_dir.resolve())
        return 0
    except Exception as exc:
        print(f"GEDI analysis failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
