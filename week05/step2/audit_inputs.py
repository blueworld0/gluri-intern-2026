"""Read-only audit of the 2025 Sentinel-2 TIFFs used by Week04."""

import json
import os
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol


BANDS = ("B02", "B03", "B04", "B08", "B11", "B12")
SOURCE = ROOT / "week03/Sentinel-2/2025-04-28"
POINTS = ROOT / "week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv"
OLD_STATS = ROOT / "week04/outputs/reflectance_validation/band_statistics.csv"


def main():
    points = pd.read_csv(POINTS)
    if len(points) != 90 or points.point_id.duplicated().any():
        raise ValueError("Expected 90 unique Week04 input points")
    rows, cols = points.raster_row.to_numpy(int), points.raster_col.to_numpy(int)
    arrays, valid_masks, stats, metadata, point_raw = {}, {}, [], [], {}
    grid = None
    for band in BANDS:
        paths = list(SOURCE.glob(f"*_{band}_(Raw).tiff"))
        if len(paths) != 1:
            raise ValueError(f"Expected one 2025 {band} TIFF, found {len(paths)}")
        path = paths[0]
        with rasterio.open(path) as ds:
            current_grid = (ds.crs, ds.transform, ds.shape)
            if grid is None:
                grid = current_grid
            elif current_grid != grid:
                raise ValueError("2025 band TIFF grids differ")
            if ds.count != 2 or ds.dtypes != ("uint16", "uint16"):
                raise ValueError(f"Unexpected TIFF structure: {path}")
            raw, mask = ds.read(1), ds.read(2)
            if not np.isin(mask, (0, 1)).all():
                raise ValueError(f"Second band is not binary dataMask: {path}")
            valid = (mask == 1) & (ds.read_masks(1) > 0)
            if ds.nodatavals[0] is not None:
                valid &= raw != ds.nodatavals[0]
            vals = raw[valid]
            if not vals.size:
                raise ValueError(f"No valid pixels: {path}")
            if not valid[rows, cols].all():
                raise ValueError(f"GEDI point samples invalid in {band}")
            arrays[band] = raw
            valid_masks[band] = valid
            point_raw[band] = raw[rows, cols]
            stats.append({"band": band, "total_pixels": raw.size,
                          "valid_pixels": vals.size, "data_mask_zero_n": int((mask == 0).sum()),
                          "valid_zero_n": int((vals == 0).sum()),
                          "valid_65535_n": int((vals == 65535).sum()),
                          "valid_gt_32767_n": int((vals > 32767).sum()),
                          "raw_min": int(vals.min()), "raw_max": int(vals.max()),
                          "raw_mean": float(vals.mean()), "raw_median": float(np.median(vals)),
                          "if_div_65535_min": float(vals.min() / 65535),
                          "if_div_65535_max": float(vals.max() / 65535),
                          "if_div_65535_mean": float(vals.mean() / 65535),
                          "if_div_65535_median": float(np.median(vals) / 65535),
                          "geotiff_scale": float(ds.scales[0]),
                          "geotiff_offset": float(ds.offsets[0]),
                          "geotiff_nodata": ds.nodatavals[0],
                          "dtype": ds.dtypes[0]})
            metadata.append({"band": band, "path": str(path.relative_to(ROOT)),
                             "crs": str(ds.crs), "shape": list(ds.shape),
                             "resolution_m": list(ds.res), "bounds": list(ds.bounds),
                             "scales": list(ds.scales), "offsets": list(ds.offsets),
                             "nodatavals": list(ds.nodatavals), "units": list(ds.units),
                             "descriptions": list(ds.descriptions),
                             "dataset_tags": ds.tags(), "band_tags": [ds.tags(1), ds.tags(2)]})
    stat_df = pd.DataFrame(stats)
    old = pd.read_csv(OLD_STATS).query("year == 2025").set_index("band")
    for row in stat_df.itertuples():
        prior = old.loc[row.band]
        if not np.isclose([row.raw_min, row.raw_max, row.raw_mean, row.raw_median],
                          [prior["min"], prior["max"], prior["mean"], prior["median"]],
                          rtol=0, atol=1e-9).all():
            raise ValueError(f"Week04 statistics differ for {row.band}")
    stat_df.to_csv(OUT / "band_quality.csv", index=False)
    common_valid = np.logical_and.reduce(list(valid_masks.values()))
    full = np.column_stack([arrays[band][common_valid].astype(float) for band in BANDS])
    point = np.column_stack([point_raw[band].astype(float) for band in BANDS])
    corrs = []
    for name, matrix in (("all_common_valid_pixels", full), ("week04_90_points", point)):
        corr = np.corrcoef(matrix, rowvar=False)
        for i, a in enumerate(BANDS):
            for j, b in enumerate(BANDS):
                corrs.append({"sample": name, "band_a": a, "band_b": b,
                              "pearson_r": float(corr[i, j]), "n": len(matrix)})
    pd.DataFrame(corrs).to_csv(OUT / "band_correlations.csv", index=False)
    with rasterio.open(next(SOURCE.glob("*_B04_(Raw).tiff"))) as ds:
        xs, ys = ds.transform * (cols + .5, rows + .5)
        sample_rows, sample_cols = rowcol(ds.transform, points.easting_m.to_numpy(),
                                           points.northing_m.to_numpy())
    center_distance = np.hypot(points.easting_m.to_numpy() - xs, points.northing_m.to_numpy() - ys)
    if not np.array_equal(sample_rows, rows) or not np.array_equal(sample_cols, cols):
        raise ValueError("Saved raster row/column differs from current TIFF transform")
    index_stats = []
    for name in ("ndvi", "evi", "ndmi"):
        v = points[name].to_numpy(float)
        index_stats.append({"index": name.upper(), "n": len(v), "min": float(v.min()),
                            "max": float(v.max()), "mean": float(v.mean()),
                            "median": float(np.median(v)),
                            "outside_minus1_plus1_n": int(((v < -1) | (v > 1)).sum()),
                            "abs_gt_2_n": int((abs(v) > 2).sum())})
    pd.DataFrame(index_stats).to_csv(OUT / "index_quality.csv", index=False)
    metadata_summary = {"files": metadata,
                        "common_valid_pixel_n": int(common_valid.sum()),
                        "point_count": len(points),
                        "point_data_mask_invalid_n": 0,
                        "point_raw_65535_any_n": int(np.any(point == 65535, axis=1).sum()),
                        "point_raw_zero_any_n": int(np.any(point == 0, axis=1).sum()),
                        "point_raster_pixel_duplicates_n": int(points.duplicated(["raster_row", "raster_col"]).sum()),
                        "coordinate_difference_m_max": float(points.coordinate_difference_m.max()),
                        "point_to_pixel_center_distance_m": {
                            "min": float(center_distance.min()),
                            "median": float(np.median(center_distance)),
                            "max": float(center_distance.max())},
                        "nominal_date_gap_days": sorted(points.date_gap_days.unique().astype(int).tolist()),
                        "scl_native_values": points.scl_native_20m.value_counts().to_dict(),
                        "scl_aligned_values": points.scl_aligned_10m.value_counts().to_dict(),
                        "scl_cloud_or_shadow_n": int(points.scl_cloud_or_shadow.sum()),
                        "scl_boundary_disagreement_n": int(points.scl_boundary_disagreement.sum())}
    (OUT / "audit_metadata.json").write_text(json.dumps(metadata_summary, indent=2), encoding="utf-8")
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for ax, band in zip(axes.flat, BANDS):
        vals = arrays[band][common_valid].astype(float) / 65535
        ax.hist(vals, bins=np.linspace(0, 1, 51), color="#3a719e", edgecolor="white")
        ax.set(title=band, xlabel="Raw TIFF value / 65535 (assumption)", ylabel="Pixels")
        ax.grid(axis="y", alpha=.2)
    fig.savefig(OUT / "band_distributions.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, (name, matrix) in zip(axes, (("All common valid pixels", full), ("Week04 points", point))):
        corr = np.corrcoef(matrix, rowvar=False)
        im = ax.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
        ax.set(title=f"{name} (n={len(matrix):,})", xticks=range(6), yticks=range(6),
               xticklabels=BANDS, yticklabels=BANDS)
        for i in range(6):
            for j in range(6):
                ax.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=axes, shrink=.8, label="Pearson r")
    fig.savefig(OUT / "band_correlations.png", dpi=160)
    plt.close(fig)
    print(stat_df[["band", "valid_pixels", "raw_min", "raw_max", "raw_mean", "raw_median", "valid_zero_n", "valid_65535_n"]].to_string(index=False))
    print(json.dumps({k: v for k, v in metadata_summary.items() if k != "files"}, indent=2))


if __name__ == "__main__":
    main()
