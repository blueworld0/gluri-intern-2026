"""Prepare the 2025 Sentinel-2 export and pair it with the same GEDI shots as 2026.

The TIFFs have no radiometric or acquisition-time tags. The 1/65535 scale and
2025-04-28 date are documented working assumptions, not verified metadata.
"""

import argparse
import json
import os
from pathlib import Path
import re
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol
from rasterio.warp import transform as project

from main import ROOT, save_raster, statistics
from merge_sentinel_gedi import project_and_sample, resolve_columns
from src.inspect_data import EXPECTED, inspect_directory
from src.preprocess import prepare_bands
from src.vegetation_indices import calculate
from src.visualize import save_index_map, save_ndvi_histogram, save_rgb


INPUT = ROOT / "Sentinel-2" / "2025-04-28"
GEDI_CSV = ROOT / "outputs" / "gedi" / "gedi_valid_points.csv"
OLD_CSV = ROOT / "outputs" / "merged" / "sentinel_gedi_merged.csv"
OUTPUT = ROOT / "outputs" / "sentinel_2025"
SCALE = 1 / 65535  # Provisional: same 16-bit Analytical Raw export as the 2026 inputs.
INDICES = ("ndvi", "evi", "ndmi")
DATE_PATTERN = re.compile(r"(?<!\d)(20\d{2}-\d{2}-\d{2})(?!\d)")


def checked_source_date(sources):
    labels = {source.date_label for source in sources.values()}
    if len(labels) != 1:
        raise ValueError(f"Source TIFF date intervals disagree: {labels}")
    label = next(iter(labels))
    if not label:
        raise ValueError("No date in TIFF metadata or filenames")
    dates = set(DATE_PATTERN.findall(label))
    if len(dates) != 1:
        raise ValueError(f"TIFF request spans multiple dates: {label}")
    return dates.pop(), all(source.date_is_exact for source in sources.values())


def source_diagnostics(sources, prepared):
    items = {}
    for name, source in sources.items():
        with rasterio.open(source.path) as ds:
            raw = ds.read(source.band_index)
            mask = ds.read(2) if ds.count == 2 else None
            if mask is None or not np.isin(mask, [0, 1]).all():
                raise ValueError(f"{source.path.name}: expected a binary dataMask in band 2")
            valid = (mask == 1) & (ds.read_masks(1) > 0)
            values = prepared.arrays[name]
            finite = values[np.isfinite(values)]
            items[name] = {
                "file": source.path.name, "size": [ds.width, ds.height],
                "crs": str(ds.crs), "resolution_m": list(ds.res),
                "bounds": list(ds.bounds), "transform": list(ds.transform)[:6],
                "dtype": ds.dtypes[0], "nodata": ds.nodatavals[0],
                "geotiff_scale": ds.scales[0], "geotiff_offset": ds.offsets[0],
                "data_mask_zero_pixels": int((mask == 0).sum()),
                "valid_pixels": int(valid.sum()),
                "saturated_65535_pixels": int(((raw == 65535) & valid).sum()),
                "reflectance_min": float(finite.min()),
                "reflectance_median": float(np.median(finite)),
                "reflectance_p99": float(np.percentile(finite, 99)),
                "reflectance_max": float(finite.max()),
            }
    return items


def sample_quality(frame, sources, reference):
    """Sample masks and raw saturation separately from index raster validity."""
    result = pd.DataFrame(index=frame.index)
    result["data_mask_valid"] = True
    result["raw_65535_any"] = False
    result["raw_zero_any"] = False
    inside = frame["inside_raster_extent"].to_numpy(dtype=bool)
    rows = frame.loc[inside, "raster_row"].to_numpy(dtype=int)
    cols = frame.loc[inside, "raster_col"].to_numpy(dtype=int)
    result.loc[~inside, "data_mask_valid"] = False
    for name, source in sources.items():
        with rasterio.open(source.path) as ds:
            if not (ds.crs == reference["crs"] and ds.transform.almost_equals(reference["transform"])
                    and ds.width == reference["width"] and ds.height == reference["height"]):
                raise ValueError(f"{name} grid changed; direct quality sampling is unsafe")
            data = ds.read(1)
            mask = ds.read(2)
            valid = np.zeros(len(frame), dtype=bool)
            saturated = np.zeros(len(frame), dtype=bool)
            zero = np.zeros(len(frame), dtype=bool)
            valid[inside] = (mask[rows, cols] == 1) & (ds.read_masks(1)[rows, cols] > 0)
            saturated[inside] = data[rows, cols] == 65535
            zero[inside] = data[rows, cols] == 0
            result[f"{name}_data_mask_valid"] = valid
            result[f"{name}_raw_65535"] = saturated
            result["data_mask_valid"] &= valid
            result["raw_65535_any"] |= saturated
            result["raw_zero_any"] |= zero
    return result


def optional_scl(frame, reference):
    candidates = sorted(p for p in INPUT.iterdir() if p.suffix.lower() in {".tif", ".tiff"}
                        and re.search(r"(?<![A-Za-z0-9])SCL(?![A-Za-z0-9])", p.stem, re.I))
    statuses = pd.Series("unassessed_no_scl", index=frame.index, dtype="string")
    classes = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    if not candidates:
        return statuses, classes, "No SCL raster supplied; dataMask cannot establish cloud or shadow absence"
    if len(candidates) != 1:
        raise ValueError(f"Expected at most one SCL raster, found {len(candidates)}")
    with rasterio.open(candidates[0]) as ds:
        if ds.crs is None or ds.count < 1:
            raise ValueError("SCL raster lacks CRS or data band")
        if ds.crs == reference["crs"]:
            x, y = frame.easting_m.to_numpy(), frame.northing_m.to_numpy()
        else:
            x, y = project(reference["crs"], ds.crs, frame.easting_m.tolist(), frame.northing_m.tolist())
        rr, cc = rowcol(ds.transform, x, y)
        rr, cc = np.asarray(rr), np.asarray(cc)
        inside = (rr >= 0) & (rr < ds.height) & (cc >= 0) & (cc < ds.width)
        values = ds.read(1)
        sampled = np.full(len(frame), -1, dtype=int)
        sampled[inside] = values[rr[inside], cc[inside]]
        classes.loc[inside] = sampled[inside]
        statuses.loc[inside] = "clear_or_other"
        statuses.loc[np.isin(sampled, [3, 8, 9, 10])] = "cloud_or_shadow"
        statuses.loc[np.isin(sampled, [0, 1])] = "scl_invalid"
        statuses.loc[sampled == 7] = "scl_unclassified"
        statuses.loc[sampled == 11] = "snow_ice"
        statuses.loc[~inside] = "outside_scl"
        return statuses, classes, f"SCL sampled from {candidates[0].name} with nearest pixel"


def paired_results(new, old, nominal_date, exact):
    keys = ["gedi_file", "beam", "shot_number"]
    old = old.copy()
    old["shot_number"] = old.shot_number.astype("string")
    new["shot_number"] = new.shot_number.astype("string")
    if old.duplicated(keys).any() or new.duplicated(keys).any():
        raise ValueError("Duplicate GEDI identity in a merge table")
    if set(map(tuple, old[keys].to_numpy())) != set(map(tuple, new[keys].to_numpy())):
        raise ValueError("The two merge tables do not contain the same GEDI points")
    old_columns = keys + ["merge_status", "date_gap_days"] + list(INDICES)
    old_subset = old[old_columns].rename(columns={
        "merge_status": "merge_status_2026", "date_gap_days": "date_gap_days_2026",
        **{name: f"{name}_2026" for name in INDICES},
    })
    paired = new.merge(old_subset, on=keys, validate="one_to_one", how="left")
    paired = paired.rename(columns={"merge_status": "merge_status_2025",
                                    **{name: f"{name}_2025" for name in INDICES}})
    paired["sentinel_2025_request_date"] = nominal_date
    paired["date_gap_days_2025_nominal"] = paired.pop("date_gap_days")
    paired = paired.drop(columns=["sentinel_acquisition_date"])
    paired["sentinel_2025_acquisition_verified"] = exact
    for name in INDICES:
        paired[f"{name}_2025_minus_2026"] = paired[f"{name}_2025"] - paired[f"{name}_2026"]
    return paired


def comparison_table(paired):
    rows = []
    for name in INDICES:
        common = paired[[f"{name}_2025", f"{name}_2026", "agbd"]].apply(pd.to_numeric, errors="coerce").dropna()
        for year in (2025, 2026):
            values = common[f"{name}_{year}"]
            rows.append({"index": name.upper(), "year": year, "common_n": len(common),
                         "min": values.min(), "max": values.max(), "mean": values.mean(),
                         "median": values.median(), "std": values.std(ddof=0),
                         "pearson_vs_agbd": values.corr(common.agbd, method="pearson"),
                         "spearman_vs_agbd": values.rank().corr(common.agbd.rank(), method="pearson")})
    return pd.DataFrame(rows)


def comparison_figure(paired, comparison, folder):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    for ax, name in zip(axes, INDICES):
        subset = paired[[f"{name}_2025", f"{name}_2026"]].apply(pd.to_numeric, errors="coerce").dropna()
        ax.scatter(subset[f"{name}_2026"], subset[f"{name}_2025"], s=22, alpha=.75)
        if len(subset):
            low = min(subset.min()); high = max(subset.max())
            ax.plot([low, high], [low, high], color="gray", linestyle="--", linewidth=1)
        ax.set_title(f"{name.upper()} (same GEDI points; n={len(subset)})")
        ax.set_xlabel("2026-09-12 index")
        ax.set_ylabel("2025-04-28 request index")
        ax.grid(alpha=.2)
    fig.savefig(folder / "paired_indices_2025_vs_2026.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    for ax, name in zip(axes, INDICES):
        for year, color in ((2025, "#228833"), (2026, "#4477aa")):
            pair = paired[[f"{name}_{year}", "agbd"]].apply(pd.to_numeric, errors="coerce").dropna()
            ax.scatter(pair[f"{name}_{year}"], pair.agbd, s=22, alpha=.65,
                       color=color, label=f"{year} (n={len(pair)})")
        ax.set_title(f"{name.upper()} vs GEDI predicted AGBD")
        ax.set_xlabel("Index")
        ax.set_ylabel("AGBD (Mg/ha)")
        ax.grid(alpha=.2)
        ax.legend()
    fig.savefig(folder / "indices_vs_agbd_by_year.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    x = np.arange(len(INDICES))
    for method, offset, color in (("pearson", -.19, "#4477aa"), ("spearman", .19, "#228833")):
        for year, dx, alpha in ((2025, -.07, .9), (2026, .07, .55)):
            values = [comparison.loc[(comparison["index"] == name.upper()) &
                                     (comparison.year == year), f"{method}_vs_agbd"].iloc[0]
                      for name in INDICES]
            ax.bar(x + offset + dx, values, width=.14, color=color, alpha=alpha,
                   label=f"{method.title()} {year}")
    ax.axhline(0, color="black", linewidth=.6)
    ax.set_xticks(x, [name.upper() for name in INDICES])
    ax.set_ylabel("Correlation with GEDI predicted AGBD")
    ax.set_title("Correlations on paired GEDI points (exploratory)")
    ax.legend(ncol=2)
    fig.savefig(folder / "correlations_by_year.png", dpi=160)
    plt.close(fig)


def footprint_figure(frame, raw_arrays, reference, folder):
    rgb = np.stack([raw_arrays[band] for band in ("B04", "B03", "B02")], axis=-1)
    valid = np.all(np.isfinite(rgb), axis=-1)
    display = np.zeros_like(rgb, dtype=np.float32)
    for channel in range(3):
        low, high = np.percentile(rgb[:, :, channel][valid], [2, 98])
        display[:, :, channel] = np.clip((rgb[:, :, channel] - low) / (high - low), 0, 1)
    display[~valid] = 0
    bounds = reference["bounds"]
    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    ax.imshow(display, extent=(bounds.left, bounds.right, bounds.bottom, bounds.top), origin="upper")
    for mask, color, label in (
        (~frame.quality_review_needed, "#00ff99", "No recorded quality issue"),
        (frame.quality_review_needed, "#ffcc00", "Quality review required"),
    ):
        if mask.any():
            ax.scatter(frame.loc[mask, "easting_m"], frame.loc[mask, "northing_m"],
                       s=22, facecolors="none", edgecolors=color, linewidths=.8,
                       label=f"{label} (n={int(mask.sum())})")
    ax.set_title("GEDI points on 2025 Sentinel-2 RGB\nSCL unavailable: cloud status unassessed")
    ax.set_xlabel(f"Easting (m, {reference['crs']})")
    ax.set_ylabel("Northing (m)")
    ax.legend(loc="lower left")
    fig.savefig(folder / "gedi_rgb_quality_overlay.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {out}")
    stage = "input inspection"
    try:
        if not GEDI_CSV.exists() or not OLD_CSV.exists():
            raise FileNotFoundError("Existing GEDI and 2026 merge CSVs are required")
        sources, warnings = inspect_directory(INPUT)
        if set(sources) != set(EXPECTED):
            raise ValueError(f"Expected six identifiable bands; found {sorted(sources)}")
        nominal_date, exact = checked_source_date(sources)
        if nominal_date != "2025-04-28":
            raise ValueError(f"Unexpected source date: {nominal_date}")
        print(f"Date: {nominal_date}; acquisition metadata verified={exact}")
        stage = "preprocessing and index calculation"
        prepared = prepare_bands(sources, SCALE, 0.0)
        diagnostics = source_diagnostics(sources, prepared)
        indices, skipped = calculate(prepared.arrays)
        if skipped:
            raise ValueError(f"Indices could not be calculated: {skipped}")
        stage = "GEDI merge"
        gedi = pd.read_csv(GEDI_CSV, dtype={"shot_number": "string"})
        old = pd.read_csv(OLD_CSV, dtype={"shot_number": "string"})
        if len(gedi) != 97 or len(old) != len(gedi) or not gedi.quality_valid.eq(True).all():
            raise ValueError("GEDI input no longer matches the 97 quality-valid points")
        ref = sources["B04"]
        with rasterio.open(ref.path) as ds:
            reference = {"crs": ds.crs, "transform": ds.transform,
                         "width": ds.width, "height": ds.height, "bounds": ds.bounds}
        arrays = {name: (data.astype(np.float64), np.isfinite(data)) for name, data in indices.items()}
        new = project_and_sample(gedi, resolve_columns(gedi), reference, arrays, nominal_date)
        quality = sample_quality(new, sources, reference)
        new = pd.concat([new, quality], axis=1)
        new["scl_status"], new["scl_class"], scl_note = optional_scl(new, reference)
        new["nodata_affected"] = ~new["data_mask_valid"] | new[list(INDICES)].isna().any(axis=1)
        new["quality_review_needed"] = (new["nodata_affected"] | new["raw_65535_any"]
                                        | new["scl_status"].ne("clear_or_other"))
        paired = paired_results(new, old, nominal_date, exact)
        comparison = comparison_table(paired)
        gap = (pd.Timestamp(nominal_date, tz="UTC") -
               pd.to_datetime(gedi.observed_at_utc, utc=True).dt.normalize()).dt.days
        if not gap.eq(3).all():
            raise ValueError(f"Unexpected GEDI to Sentinel request-date gaps: {gap.value_counts().to_dict()}")
        stage = "saving rasters, tables, and figures"
        raster_dir, figure_dir = out / "rasters", out / "figures"
        raster_dir.mkdir(parents=True)
        figure_dir.mkdir()
        for name, data in indices.items():
            path = raster_dir / f"{name}.tif"
            save_raster(path, data, prepared.profile, name, SCALE, 0.0,
                        nominal_date if exact else None)
            with rasterio.open(path, "r+") as ds:
                ds.update_tags(requested_date=nominal_date,
                               acquisition_date_verified=str(exact).lower(),
                               reflectance_scale_status="provisional_16bit_analytical_export")
            save_index_map(name, data, prepared.profile, figure_dir, nominal_date, exact)
        save_rgb(prepared.raw_arrays, prepared.profile, figure_dir, nominal_date, exact)
        save_ndvi_histogram(indices["ndvi"], figure_dir, nominal_date, exact)
        pd.DataFrame([statistics(name, data) for name, data in indices.items()]).to_csv(
            out / "raster_statistics.csv", index=False)
        paired.to_csv(out / "paired_2025_2026.csv", index=False)
        new.to_csv(out / "gedi_sentinel_2025_merged.csv", index=False)
        comparison.to_csv(out / "index_comparison.csv", index=False)
        new.loc[new.quality_review_needed].to_csv(out / "quality_review_points.csv", index=False)
        comparison_figure(paired, comparison, figure_dir)
        footprint_figure(new, prepared.raw_arrays, reference, figure_dir)
        with rasterio.open(ROOT / "outputs" / "rasters" / "ndvi.tif") as ds:
            old_grid = {"crs": str(ds.crs), "resolution_m": list(ds.res),
                        "bounds": list(ds.bounds), "shape": [ds.height, ds.width]}
        report = {
            "date": {"gedi_utc": "2025-04-25", "sentinel_request": nominal_date,
                     "acquisition_verified": exact, "nominal_gap_days": sorted(gap.unique().tolist())},
            "calibration": {"applied_scale": SCALE, "applied_offset": 0,
                            "status": "provisional; TIFF has no radiometric tags",
                            "basis": "same named 16-bit Analytical Raw export convention as 2026"},
            "grids": {"sentinel_2025": {"crs": str(reference["crs"]),
                                        "bounds": list(reference["bounds"]),
                                        "shape": [reference["height"], reference["width"]]},
                      "sentinel_2026": old_grid},
            "bands": diagnostics,
            "counts": {"gedi_input": len(new), "inside_2025_extent": int(new.inside_raster_extent.sum()),
                       "merge_success_2025": int(new.merge_status.eq("success").sum()),
                       "merge_success_2026": int(old.merge_status.eq("success").sum()),
                       "missing_any_index_2025": int(new[list(INDICES)].isna().any(axis=1).sum()),
                       "nodata_affected_shots": int(new.nodata_affected.sum()),
                       "saturated_65535_shots": int(new.raw_65535_any.sum()),
                       "cloud_shadow_confirmed_shots": int(new.scl_status.eq("cloud_or_shadow").sum()),
                       "cloud_shadow_unassessed_shots": int(new.scl_status.eq("unassessed_no_scl").sum())},
            "scl": scl_note,
            "index_extremes": {name: {"outside_minus1_plus1_pixels": int((np.abs(data[np.isfinite(data)]) > 1).sum()),
                                     "abs_gt2_pixels": int((np.abs(data[np.isfinite(data)]) > 2).sum())}
                               for name, data in indices.items()},
            "warnings": warnings + prepared.warnings,
        }
        (out / "quality_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"Saved {out}; merged={report['counts']['merge_success_2025']}/{len(new)}; "
              f"nominal date gap={report['date']['nominal_gap_days']}; SCL={scl_note}")
        print(comparison.to_string(index=False))
        return 0
    except Exception as exc:
        print(f"ERROR during {stage}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
