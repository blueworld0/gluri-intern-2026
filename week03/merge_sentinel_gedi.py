"""Sample Sentinel-2 indices at quality-valid GEDI04_A footprints."""

import argparse
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol
from rasterio.warp import transform


ROOT = Path(__file__).resolve().parent
INPUT_CSV = ROOT / "outputs" / "gedi" / "gedi_valid_points.csv"
RASTER_DIR = ROOT / "outputs" / "rasters"
DEFAULT_OUTPUT = ROOT / "outputs" / "merged"
INDICES = ("ndvi", "evi", "ndmi")
ALIASES = {
    "latitude": ("latitude", "lat", "lat_lowestmode"),
    "longitude": ("longitude", "lon", "lon_lowestmode"),
    "agbd": ("agbd", "agbd_Mg_ha"),
    "agbd_se": ("agbd_se", "agbd_se_Mg_ha"),
    "observed_at_utc": ("observed_at_utc", "observation_time_utc", "delta_time_utc"),
    "source_file": ("source_file", "gedi_file"),
    "beam": ("beam", "beam_id"),
    "shot_number": ("shot_number", "point_id"),
}


def resolve_columns(frame: pd.DataFrame) -> dict[str, str]:
    found = {}
    for logical, names in ALIASES.items():
        matches = [name for name in names if name in frame.columns]
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one column for {logical}; found {matches}")
        found[logical] = matches[0]
    for flag in ("l2a_quality_flag_rel3", "l4a_quality_flag_rel3", "degrade_include_flag",
                 "elev_highestreturn_outlier_flag", "quality_valid"):
        if flag not in frame.columns:
            raise ValueError(f"Missing GEDI quality column: {flag}")
    return found


def inspect_rasters() -> tuple[dict, dict, str]:
    arrays = {}
    metadata = {}
    reference = None
    dates = set()
    for name in INDICES:
        path = RASTER_DIR / f"{name}.tif"
        if not path.exists():
            raise FileNotFoundError(path)
        with rasterio.open(path) as ds:
            if ds.crs is None or ds.count != 1:
                raise ValueError(f"{path.name}: expected one band with a CRS")
            if reference is None:
                reference = {"crs": ds.crs, "transform": ds.transform,
                             "width": ds.width, "height": ds.height, "bounds": ds.bounds}
            else:
                same = (ds.crs == reference["crs"] and ds.transform.almost_equals(reference["transform"])
                        and ds.width == reference["width"] and ds.height == reference["height"]
                        and np.allclose(tuple(ds.bounds), tuple(reference["bounds"]), atol=1e-6))
                if not same:
                    raise ValueError(f"{path.name}: raster grid differs from NDVI")
            data = ds.read(1, masked=True)
            valid = (~np.ma.getmaskarray(data)) & np.isfinite(np.asarray(data))
            arrays[name] = (np.asarray(data, dtype=np.float64), valid)
            acquisition_date = ds.tags().get("acquisition_date")
            if not acquisition_date:
                raise ValueError(f"{path.name}: missing acquisition_date tag")
            dates.add(acquisition_date)
            metadata[name] = {
                "file": str(path.relative_to(ROOT)), "crs": str(ds.crs),
                "width": ds.width, "height": ds.height, "resolution_m": list(ds.res),
                "bounds": list(ds.bounds), "transform": list(ds.transform)[:6],
                "nodata": ds.nodata, "dtype": ds.dtypes[0],
                "valid_pixel_count": int(valid.sum()),
                "valid_pixel_fraction": float(valid.mean()),
                "acquisition_date": acquisition_date,
            }
            print(f"{name.upper()}: {ds.width}x{ds.height}, {ds.crs}, res={ds.res}, "
                  f"NoData={ds.nodata}, dtype={ds.dtypes[0]}, valid={valid.mean():.2%}")
    if len(dates) != 1:
        raise ValueError(f"Index raster acquisition dates differ: {dates}")
    return reference, arrays, metadata, next(iter(dates))


def project_and_sample(frame: pd.DataFrame, columns: dict, reference: dict,
                       arrays: dict, sentinel_date: str) -> pd.DataFrame:
    merged = frame.copy()
    merged = merged.rename(columns={columns["agbd"]: "agbd", columns["agbd_se"]: "agbd_se",
                                    columns["source_file"]: "gedi_file"})
    lat = pd.to_numeric(merged[columns["latitude"]], errors="coerce").to_numpy(dtype=float)
    lon = pd.to_numeric(merged[columns["longitude"]], errors="coerce").to_numpy(dtype=float)
    n = len(merged)
    coordinate_valid = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
    east = np.full(n, np.nan)
    north = np.full(n, np.nan)
    if coordinate_valid.any():
        x, y = transform("EPSG:4326", reference["crs"],
                         lon[coordinate_valid].tolist(), lat[coordinate_valid].tolist())
        east[coordinate_valid], north[coordinate_valid] = x, y
    merged["easting_m"], merged["northing_m"] = east, north
    bounds = reference["bounds"]
    inside = coordinate_valid & (east >= bounds.left) & (east < bounds.right)
    inside &= (north > bounds.bottom) & (north <= bounds.top)
    rows = np.full(n, -1, dtype=int)
    cols = np.full(n, -1, dtype=int)
    if inside.any():
        rr, cc = rowcol(reference["transform"], east[inside], north[inside])
        rows[inside], cols[inside] = rr, cc
    in_pixels = inside & (rows >= 0) & (rows < reference["height"])
    in_pixels &= (cols >= 0) & (cols < reference["width"])
    merged["raster_row"] = pd.Series(rows).where(in_pixels, pd.NA).astype("Int64")
    merged["raster_col"] = pd.Series(cols).where(in_pixels, pd.NA).astype("Int64")
    merged["inside_raster_extent"] = in_pixels
    if {"easting_m", "northing_m"}.issubset(frame.columns):
        previous_e = pd.to_numeric(frame["easting_m"], errors="coerce").to_numpy(dtype=float)
        previous_n = pd.to_numeric(frame["northing_m"], errors="coerce").to_numpy(dtype=float)
        merged["coordinate_difference_m"] = np.hypot(east - previous_e, north - previous_n)
    keys = pd.Series(list(zip(rows[in_pixels], cols[in_pixels])))
    counts = keys.value_counts()
    pixel_counts = np.zeros(n, dtype=int)
    pixel_counts[in_pixels] = [int(counts[(r, c)]) for r, c in zip(rows[in_pixels], cols[in_pixels])]
    merged["pixel_observation_count"] = pixel_counts
    merged["duplicate_pixel"] = pixel_counts > 1
    merged["point_id"] = (merged["gedi_file"].astype(str) + ":" + merged[columns["beam"]].astype(str)
                          + ":" + merged[columns["shot_number"]].astype(str))
    if merged["point_id"].duplicated().any():
        raise ValueError("GEDI point identifiers are not unique")
    observed = pd.to_datetime(merged[columns["observed_at_utc"]], utc=True, errors="coerce")
    sentinel = pd.Timestamp(sentinel_date, tz="UTC")
    merged["gedi_observation_date"] = observed.dt.strftime("%Y-%m-%d")
    merged["sentinel_acquisition_date"] = sentinel_date
    merged["date_gap_days"] = (sentinel - observed.dt.normalize()).dt.days.astype("Int64")
    reasons = [[] for _ in range(n)]
    for i in range(n):
        if not coordinate_valid[i]:
            reasons[i].append("invalid_coordinate")
        elif not in_pixels[i]:
            reasons[i].append("outside_raster_extent")
        if pd.isna(observed.iloc[i]):
            reasons[i].append("invalid_observation_time")
        agbd_value = pd.to_numeric(merged["agbd"], errors="coerce").iloc[i]
        se_value = pd.to_numeric(merged["agbd_se"], errors="coerce").iloc[i]
        if not np.isfinite(agbd_value) or agbd_value < 0:
            reasons[i].append("invalid_agbd")
        if not np.isfinite(se_value) or se_value < 0:
            reasons[i].append("invalid_agbd_se")
    for name, (data, valid) in arrays.items():
        sampled = np.full(n, np.nan)
        sampled[in_pixels] = data[rows[in_pixels], cols[in_pixels]]
        sample_valid = np.zeros(n, dtype=bool)
        sample_valid[in_pixels] = valid[rows[in_pixels], cols[in_pixels]]
        sample_valid &= np.isfinite(sampled)
        if name in ("ndvi", "ndmi"):
            sample_valid &= (sampled >= -1) & (sampled <= 1)
        for i in np.flatnonzero(in_pixels & ~sample_valid):
            reasons[i].append(f"{name}_nodata_or_invalid")
        sampled[~sample_valid] = np.nan
        merged[name] = sampled
        merged[f"{name}_sample_valid"] = sample_valid
    merged["failure_reason"] = [";".join(parts) if parts else "none" for parts in reasons]
    merged["merge_status"] = np.where(merged["failure_reason"] == "none", "success", "failed")
    return merged


def descriptive_and_correlations(clean: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows = []
    correlations = {}
    for name in ("ndvi", "evi", "ndmi", "agbd", "agbd_se"):
        values = pd.to_numeric(clean[name], errors="coerce").dropna()
        unit = "Mg/ha" if name in ("agbd", "agbd_se") else "unitless"
        for metric, value in (("count", len(values)), ("min", values.min()), ("max", values.max()),
                              ("mean", values.mean()), ("median", values.median()), ("std", values.std(ddof=0))):
            rows.append({"variable": name, "metric": metric, "value": value, "n": len(values), "unit": unit})
    for name in INDICES:
        pair = clean[[name, "agbd"]].apply(pd.to_numeric, errors="coerce").dropna()
        pearson = pair[name].corr(pair["agbd"], method="pearson") if len(pair) >= 2 else np.nan
        ranked = pair.rank()
        spearman = ranked[name].corr(ranked["agbd"], method="pearson") if len(pair) >= 2 else np.nan
        correlations[name] = {"n": len(pair), "pearson": float(pearson), "spearman": float(spearman)}
        rows.extend([
            {"variable": name, "metric": "pearson_vs_agbd", "value": pearson, "n": len(pair), "unit": "unitless"},
            {"variable": name, "metric": "spearman_vs_agbd", "value": spearman, "n": len(pair), "unit": "unitless"},
        ])
    return pd.DataFrame(rows), correlations


def figures(merged: pd.DataFrame, clean: pd.DataFrame, reference: dict,
            arrays: dict, output: Path) -> None:
    folder = output / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    bounds = reference["bounds"]
    extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
    ndvi, valid = arrays["ndvi"]
    cmap = plt.get_cmap("RdYlGn").copy()
    cmap.set_bad("0.75")
    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    image = ax.imshow(np.ma.masked_where(~valid, ndvi), extent=extent, origin="upper",
                      cmap=cmap, vmin=-1, vmax=1)
    fig.colorbar(image, ax=ax, label="Sentinel-2 NDVI")
    if len(merged):
        shown = merged.loc[merged["inside_raster_extent"] & pd.to_numeric(merged["agbd"], errors="coerce").notna()]
        if len(shown):
            colors = shown["agbd"].astype(float).to_numpy()
            points = ax.scatter(shown["easting_m"], shown["northing_m"], c=colors,
                                cmap="plasma", s=32, edgecolors="black", linewidths=0.35,
                                label=f"GEDI footprints (n={len(shown)})")
            fig.colorbar(points, ax=ax, label="GEDI predicted AGBD (Mg/ha)")
            ax.legend(loc="lower right")
    ax.set_xlim(bounds.left, bounds.right)
    ax.set_ylim(bounds.bottom, bounds.top)
    ax.set_title("GEDI footprints on Sentinel-2 NDVI\nGeumjeongsan analysis rectangle")
    ax.set_xlabel(f"Easting (m, {reference['crs']})")
    ax.set_ylabel("Northing (m)")
    fig.savefig(folder / "ndvi_gedi_overlay.png", dpi=160)
    plt.close(fig)
    for name in INDICES:
        fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        pair = clean[[name, "agbd"]].dropna()
        if len(pair):
            ax.scatter(pair[name], pair["agbd"], alpha=0.75, s=28, edgecolors="black", linewidths=0.2,
                       label=f"GEDI footprints (n={len(pair)})")
            ax.legend()
        else:
            ax.text(0.5, 0.5, "No valid paired samples", transform=ax.transAxes,
                    ha="center", va="center")
        ax.set_title(f"{name.upper()} vs GEDI predicted AGBD")
        ax.set_xlabel(f"Sentinel-2 {name.upper()} (unitless)")
        ax.set_ylabel("GEDI predicted AGBD (Mg/ha)")
        ax.grid(alpha=0.25)
        fig.savefig(folder / f"{name}_vs_agbd.png", dpi=160)
        plt.close(fig)
    variables = ["ndvi", "evi", "ndmi", "agbd"]
    correlation = clean[variables].corr(method="pearson") if len(clean) >= 2 else pd.DataFrame(np.nan, index=variables, columns=variables)
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    image = ax.imshow(correlation.to_numpy(dtype=float), cmap="coolwarm", vmin=-1, vmax=1)
    fig.colorbar(image, ax=ax, label="Pearson r")
    ax.set_xticks(range(len(variables)), [v.upper() for v in variables])
    ax.set_yticks(range(len(variables)), [v.upper() for v in variables])
    for i in range(len(variables)):
        for j in range(len(variables)):
            value = correlation.iloc[i, j]
            ax.text(j, i, f"{value:.2f}" if np.isfinite(value) else "NA", ha="center", va="center")
    ax.set_title(f"Pairwise Pearson correlations (n={len(clean)})")
    fig.savefig(folder / "correlation_heatmap.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    gaps = merged["date_gap_days"].dropna().astype(int)
    if len(gaps):
        counts = gaps.value_counts().sort_index()
        ax.bar(counts.index.astype(str), counts.values, color="steelblue", label=f"GEDI footprints (n={len(gaps)})")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "No valid observation dates", transform=ax.transAxes,
                ha="center", va="center")
    ax.set_title("Sentinel-2 minus GEDI acquisition date")
    ax.set_xlabel("Date gap (days)")
    ax.set_ylabel("Footprint count")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(folder / "observation_date_gap.png", dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT,
                        help="New output directory; existing merge products are never overwritten")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    stage = "input validation"
    try:
        if not INPUT_CSV.exists():
            raise FileNotFoundError(INPUT_CSV)
        frame = pd.read_csv(INPUT_CSV, dtype={"shot_number": "string"})
        columns = resolve_columns(frame)
        if not frame["quality_valid"].eq(True).all():
            raise ValueError("GEDI valid-points input contains shots without quality_valid=True")
        print(f"GEDI: {INPUT_CSV.relative_to(ROOT)}, rows={len(frame)}, columns={list(frame.columns)}")
        reference, arrays, raster_meta, sentinel_date = inspect_rasters()
        targets = [output / name for name in ("sentinel_gedi_merged.csv", "sentinel_gedi_clean.csv",
                                                "merge_statistics.csv", "merge_quality_report.json")]
        targets += [output / "figures" / name for name in
                    ("ndvi_gedi_overlay.png", "ndvi_vs_agbd.png", "evi_vs_agbd.png",
                     "ndmi_vs_agbd.png", "correlation_heatmap.png", "observation_date_gap.png")]
        existing = [path for path in targets if path.exists()]
        if existing:
            raise FileExistsError(f"Refusing to overwrite existing result: {existing[0]}")
        stage = "coordinate conversion and raster sampling"
        merged = project_and_sample(frame, columns, reference, arrays, sentinel_date)
        clean = merged.loc[merged["merge_status"] == "success"].copy()
        stats, correlations = descriptive_and_correlations(clean)
        output.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output / "sentinel_gedi_merged.csv", index=False)
        clean.to_csv(output / "sentinel_gedi_clean.csv", index=False)
        stats.to_csv(output / "merge_statistics.csv", index=False)
        gaps = merged["date_gap_days"].dropna().astype(float)
        valid_xy = merged.loc[merged["inside_raster_extent"], ["easting_m", "northing_m"]].to_numpy(dtype=float)
        near_count = 0
        min_distance = None
        if len(valid_xy) > 1:
            distances = np.linalg.norm(valid_xy[:, None, :] - valid_xy[None, :, :], axis=2)
            np.fill_diagonal(distances, np.inf)
            nearest = distances.min(axis=1)
            near_count = int((nearest < 60).sum())
            min_distance = float(nearest.min())
        dates = merged["gedi_observation_date"].value_counts(dropna=False).sort_index().to_dict()
        report = {
            "input": {"gedi_file": str(INPUT_CSV.relative_to(ROOT)), "gedi_rows": len(frame),
                      "column_mapping": columns, "rasters": raster_meta},
            "counts": {"gedi_valid_input": len(frame),
                       "inside_sentinel_extent": int(merged["inside_raster_extent"].sum()),
                       **{f"{name}_sample_success": int(merged[f"{name}_sample_valid"].sum()) for name in INDICES},
                       "all_three_indices_success": int(merged[[f"{name}_sample_valid" for name in INDICES]].all(axis=1).sum()),
                       "clean_rows": len(clean), "merge_failed": int((merged["merge_status"] == "failed").sum()),
                       "missing_index_rows": int(merged[list(INDICES)].isna().any(axis=1).sum()),
                       "duplicate_point_ids": int(merged["point_id"].duplicated().sum()),
                       "duplicate_pixel_rows": int(merged["duplicate_pixel"].sum()),
                       "duplicate_pixel_groups": int(merged.loc[merged["duplicate_pixel"], ["raster_row", "raster_col"]].drop_duplicates().shape[0]),
                       "nearby_point_rows_under_60m": near_count},
            "minimum_nearest_neighbor_distance_m": min_distance,
            "failure_reasons": merged.loc[merged["merge_status"] == "failed", "failure_reason"].value_counts().to_dict(),
            "missing_by_column": {name: int(value) for name, value in merged.isna().sum().items()},
            "clean_shape": {"rows": len(clean), "columns": len(clean.columns)},
            "date_gap_days": {"min": float(gaps.min()) if len(gaps) else None,
                              "max": float(gaps.max()) if len(gaps) else None,
                              "mean": float(gaps.mean()) if len(gaps) else None,
                              "median": float(gaps.median()) if len(gaps) else None},
            "gedi_observation_date_counts": dates,
            "correlations_with_agbd": correlations,
            "diagnostics": {"evi_abs_gt_2": int((merged["evi"].abs() > 2).sum()),
                            "ndvi_outside_minus1_plus1": int(((merged["ndvi"] < -1) | (merged["ndvi"] > 1)).sum()),
                            "ndmi_outside_minus1_plus1": int(((merged["ndmi"] < -1) | (merged["ndmi"] > 1)).sum()),
                            "max_coordinate_difference_m": float(merged["coordinate_difference_m"].max()) if "coordinate_difference_m" in merged else None},
            "interpretation": "GEDI AGBD is a model prediction, not field truth; the Sentinel rectangle is not a forest boundary. Correlation is exploratory and does not establish causation or model performance.",
        }
        (output / "merge_quality_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        stage = "figure generation"
        figures(merged, clean, reference, arrays, output)
        print(f"Merged={len(merged)}, clean={len(clean)}, failures={report['counts']['merge_failed']}")
        print(f"Date gap: {report['date_gap_days']}; date counts={dates}")
        print(f"Correlations: {correlations}")
        print(f"Outputs: {output}")
        return 0
    except Exception as exc:
        print(f"ERROR during {stage}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
