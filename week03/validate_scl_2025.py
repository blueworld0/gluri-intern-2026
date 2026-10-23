"""Validate GEDI footprints using the native 20 m SCL and the 2025 index grid."""

import json
import os
from pathlib import Path
import re
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import rowcol
from rasterio.warp import reproject


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "Sentinel-2" / "2025-04-28"
SOURCE = ROOT / "outputs" / "sentinel_2025"
OUTPUT = SOURCE / "scl_validation"
EXPECTED_TILE = "T52SDE"
EXPECTED_SENSING = "20250428T020711"
PRODUCT_TIME = "2025-04-28T02:07:11.025Z"  # User-confirmed product metadata.
CLASSES = {
    0: "No data", 1: "Saturated/defective", 2: "Dark area", 3: "Cloud shadow",
    4: "Vegetation", 5: "Bare soil", 6: "Water", 7: "Unclassified/low cloud probability",
    8: "Medium cloud probability", 9: "High cloud probability", 10: "Cirrus",
    11: "Snow/ice",
}
CLOUD_SHADOW = {3, 8, 9, 10}
INDICES = ("ndvi", "evi", "ndmi")


def check_sources():
    scl_files = list(INPUT.rglob("*.jp2"))
    if len(scl_files) != 1:
        raise ValueError(f"Expected exactly one SCL JP2, found {len(scl_files)}")
    scl_path = scl_files[0]
    match = re.fullmatch(r"(T\d{2}[A-Z]{3})_(\d{8}T\d{6})_SCL_20m\.jp2", scl_path.name, re.I)
    if not match:
        raise ValueError(f"Unrecognized SCL product filename: {scl_path.name}")
    tile, sensing = match.group(1).upper(), match.group(2)
    if tile != EXPECTED_TILE or sensing != EXPECTED_SENSING:
        raise ValueError(f"SCL product mismatch: tile={tile}, sensing={sensing}")
    tiffs = sorted(p for p in INPUT.iterdir() if p.suffix.lower() in {".tif", ".tiff"})
    if len(tiffs) != 6:
        raise ValueError(f"Expected six Analytical TIFFs, found {len(tiffs)}")
    for path in tiffs:
        if "2025-04-28" not in path.name or "Sentinel-2_L2A" not in path.name:
            raise ValueError(f"TIFF date/collection contradicts SCL: {path.name}")
    with rasterio.open(scl_path) as scl:
        if scl.count != 1 or scl.crs is None or not np.allclose(scl.res, (20, 20)):
            raise ValueError("SCL is not a one-band georeferenced 20 m raster")
        scl_meta = {"file": scl_path.name, "driver": scl.driver, "dtype": scl.dtypes[0],
                    "crs": str(scl.crs), "resolution_m": list(scl.res),
                    "shape": [scl.height, scl.width], "bounds": list(scl.bounds),
                    "nodata": scl.nodata, "tile": tile, "sensing_from_filename": sensing}
        for path in tiffs:
            with rasterio.open(path) as ds:
                if ds.crs != scl.crs or ds.count != 2:
                    raise ValueError(f"TIFF CRS or band structure contradicts SCL: {path.name}")
                if not (scl.bounds.left <= ds.bounds.left < ds.bounds.right <= scl.bounds.right
                        and scl.bounds.bottom <= ds.bounds.bottom < ds.bounds.top <= scl.bounds.top):
                    raise ValueError(f"TIFF footprint outside SCL tile: {path.name}")
    # Analytical TIFFs lack a tile, exact sensing time, and product ID. User
    # provenance plus matching date/coverage supports use, but not proof.
    return scl_path, scl_meta, tiffs


def reference_and_indices():
    reference = None
    arrays = {}
    for name in INDICES:
        path = SOURCE / "rasters" / f"{name}.tif"
        with rasterio.open(path) as ds:
            if reference is None:
                reference = {"crs": ds.crs, "transform": ds.transform,
                             "height": ds.height, "width": ds.width, "bounds": ds.bounds}
            elif (ds.crs != reference["crs"] or ds.shape != (reference["height"], reference["width"])
                  or not ds.transform.almost_equals(reference["transform"])):
                raise ValueError(f"Index raster does not share 2025 grid: {path}")
            arrays[name] = ds.read(1, masked=True)
    return reference, arrays


def align_and_sample(scl_path, reference, points):
    with rasterio.open(scl_path) as src:
        if src.crs != reference["crs"]:
            raise ValueError("SCL CRS differs from target; product association needs review")
        aligned = np.zeros((reference["height"], reference["width"]), dtype=np.uint8)
        reproject(source=rasterio.band(src, 1), destination=aligned,
                  src_transform=src.transform, src_crs=src.crs,
                  dst_transform=reference["transform"], dst_crs=reference["crs"],
                  src_nodata=0, dst_nodata=0, resampling=Resampling.nearest)
        east = pd.to_numeric(points.easting_m, errors="coerce").to_numpy(float)
        north = pd.to_numeric(points.northing_m, errors="coerce").to_numpy(float)
        if not np.isfinite(east).all() or not np.isfinite(north).all():
            raise ValueError("GEDI projected coordinates contain nonfinite values")
        sr, sc = rowcol(src.transform, east, north)
        rr, cc = rowcol(reference["transform"], east, north)
        sr, sc, rr, cc = (np.asarray(a, dtype=int) for a in (sr, sc, rr, cc))
        if (np.any(sr < 0) or np.any(sr >= src.height) or np.any(sc < 0) or np.any(sc >= src.width)
                or np.any(rr < 0) or np.any(rr >= reference["height"])
                or np.any(cc < 0) or np.any(cc >= reference["width"])):
            raise ValueError("Some GEDI points lie outside the SCL or index raster")
        native = src.read(1)[sr, sc]
        at_index_pixel = aligned[rr, cc]
        if not set(np.unique(native)).issubset(CLASSES) or not set(np.unique(at_index_pixel)).issubset(CLASSES):
            raise ValueError("SCL contains unknown class codes")
        return aligned, native, at_index_pixel


def correlations(points):
    rows = []
    for name in INDICES:
        for stage, subset in (("before_scl", points),
                              ("after_scl_vegetation_both_grids", points.loc[points.scl_valid])):
            pair = subset[[name, "agbd"]].apply(pd.to_numeric, errors="coerce").dropna()
            x, y = pair[name], pair.agbd
            rows.append({"index": name.upper(), "stage": stage, "n": len(pair),
                         "index_mean": x.mean(), "index_median": x.median(),
                         "agbd_mean_Mg_ha": y.mean(),
                         "pearson_vs_agbd": x.corr(y) if len(pair) >= 2 else np.nan,
                         "spearman_vs_agbd": x.rank().corr(y.rank()) if len(pair) >= 2 else np.nan})
    return pd.DataFrame(rows)


def draw_maps(aligned, reference, points, folder):
    colors = ["#444444", "#d7191c", "#706e7b", "#5741a9", "#188037", "#dbaf5b",
              "#3264ba", "#aa99aa", "#ec8b34", "#d64226", "#ec56b9", "#d9f0fb"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-.5, 12.5, 1), len(colors))
    b = reference["bounds"]
    extent = (b.left, b.right, b.bottom, b.top)
    fig, ax = plt.subplots(figsize=(9, 10), constrained_layout=True)
    im = ax.imshow(aligned, cmap=cmap, norm=norm, extent=extent, origin="upper")
    cb = fig.colorbar(im, ax=ax, ticks=range(12), shrink=.75)
    cb.ax.set_yticklabels([f"{i}: {CLASSES[i]}" for i in range(12)])
    ax.set_title("Sentinel-2C SCL, 2025-04-28 02:07 UTC\n20 m to 2025 index grid by nearest neighbour")
    ax.set_xlabel("Easting (m, EPSG:32652)")
    ax.set_ylabel("Northing (m)")
    fig.savefig(folder / "scl_map.png", dpi=160)
    plt.close(fig)

    with rasterio.open(SOURCE / "rasters" / "ndvi.tif") as ds:
        ndvi = ds.read(1, masked=True)
    fig, ax = plt.subplots(figsize=(9, 10), constrained_layout=True)
    im = ax.imshow(ndvi, cmap="Greens", vmin=0, vmax=1, extent=extent, origin="upper")
    fig.colorbar(im, ax=ax, label="2025 NDVI (provisional reflectance scale)")
    groups = [
        (points.scl_valid, "#00bcd4", "SCL=4 on both grids"),
        (points.scl_boundary_disagreement, "#ff00aa", "20 m / 10 m SCL disagreement"),
        (~points.scl_valid & ~points.scl_boundary_disagreement, "#ff9d00", "Other SCL class"),
    ]
    for select, color, label in groups:
        if select.any():
            ax.scatter(points.loc[select, "easting_m"], points.loc[select, "northing_m"],
                       s=30, facecolors="none", edgecolors=color, linewidths=1.0,
                       label=f"{label} (n={int(select.sum())})")
    ax.set_title("GEDI footprints classified by SCL\nVegetation candidates require agreement at point and 10 m pixel")
    ax.set_xlabel("Easting (m, EPSG:32652)")
    ax.set_ylabel("Northing (m)")
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(folder / "gedi_scl_overlay.png", dpi=160)
    plt.close(fig)


def main():
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing result: {OUTPUT}")
    stage = "source and product checks"
    try:
        scl_path, scl_meta, tiffs = check_sources()
        reference, _ = reference_and_indices()
        points = pd.read_csv(SOURCE / "gedi_sentinel_2025_merged.csv",
                             dtype={"shot_number": "string"})
        if len(points) != 97 or not points.quality_valid.eq(True).all():
            raise ValueError("Expected exactly 97 unchanged GEDI quality-valid points")
        stage = "20 m SCL alignment and point extraction"
        aligned, native, at_index_pixel = align_and_sample(scl_path, reference, points)
        points["scl_native_20m"] = native.astype(int)
        points["scl_aligned_10m"] = at_index_pixel.astype(int)
        points["scl_native_label"] = points.scl_native_20m.map(CLASSES)
        points["scl_aligned_label"] = points.scl_aligned_10m.map(CLASSES)
        points["scl_boundary_disagreement"] = native != at_index_pixel
        points["scl_cloud_or_shadow"] = np.isin(native, list(CLOUD_SHADOW)) | np.isin(at_index_pixel, list(CLOUD_SHADOW))
        points["scl_vegetation_native"] = native == 4
        points["scl_vegetation_aligned"] = at_index_pixel == 4
        points["scl_valid"] = (points.scl_vegetation_native & points.scl_vegetation_aligned
                               & points.merge_status.eq("success") & ~points.scl_cloud_or_shadow)
        points["scl_exclusion_reason"] = np.select(
            [points.scl_valid, points.scl_cloud_or_shadow, points.scl_boundary_disagreement],
            ["accepted_vegetation", "cloud_or_shadow", "native_aligned_boundary_disagreement"],
            default="other_scl_class")
        count_rows = []
        for code, label in CLASSES.items():
            count_rows.append({"scl_class": code, "label": label,
                               "native_20m_point_count": int((native == code).sum()),
                               "aligned_10m_index_pixel_count": int((at_index_pixel == code).sum())})
        counts = pd.DataFrame(count_rows)
        correlation = correlations(points)
        with rasterio.open(scl_path) as ds:
            full = ds.read(1)
            product_cloud_percent_from_scl = 100 * np.isin(full, [8, 9, 10]).mean()
        observed = pd.to_datetime(points.observed_at_utc, utc=True)
        sensing = pd.Timestamp(PRODUCT_TIME)
        timedeltas = sensing - observed
        if observed.dt.date.nunique() != 1 or str(observed.dt.date.iloc[0]) != "2025-04-25":
            raise ValueError("GEDI UTC date differs from expected 2025-04-25")
        stage = "saving tables, figures, and report"
        OUTPUT.mkdir(parents=True)
        counts.to_csv(OUTPUT / "scl_class_counts.csv", index=False)
        points.to_csv(OUTPUT / "gedi_scl_points.csv", index=False)
        points.loc[points.scl_valid].to_csv(OUTPUT / "gedi_scl_valid_points.csv", index=False)
        correlation.to_csv(OUTPUT / "correlation_before_after.csv", index=False)
        draw_maps(aligned, reference, points, OUTPUT)
        report = {
            "source": scl_meta,
            "product": {"user_confirmed_satellite": "Sentinel-2C", "user_confirmed_product": "S2C_MSIL2A",
                        "user_confirmed_tile": EXPECTED_TILE, "user_confirmed_relative_orbit": "R103",
                        "user_confirmed_processing_baseline": "N05.11",
                        "user_confirmed_sensing_time_utc": PRODUCT_TIME,
                        "scl_cloud_classes_8_9_10_percent_full_tile": product_cloud_percent_from_scl},
            "tiff_association": {"status": "consistent_with_user_provenance_not_independently_proven",
                                 "reason": "TIFFs match L2A request date, CRS and SCL coverage, but carry no tile, sensing-time or product-ID tags",
                                 "analytical_tiffs": [p.name for p in tiffs]},
            "alignment": {"method": "nearest", "source_resolution_m": [20, 20],
                          "target_shape": [reference["height"], reference["width"]],
                          "target_crs": str(reference["crs"]),
                          "target_bounds": list(reference["bounds"]),
                          "native_vs_aligned_disagreement_points": int(points.scl_boundary_disagreement.sum())},
            "counts": {"gedi_quality_valid_before_scl": len(points),
                       "vegetation_native_20m": int(points.scl_vegetation_native.sum()),
                       "vegetation_aligned_10m": int(points.scl_vegetation_aligned.sum()),
                       "cloud_or_shadow_excluded": int(points.scl_cloud_or_shadow.sum()),
                       "other_or_boundary_excluded": int((~points.scl_valid & ~points.scl_cloud_or_shadow).sum()),
                       "final_vegetation_candidates": int(points.scl_valid.sum()),
                       "exclusion_reasons": points.scl_exclusion_reason.value_counts().to_dict()},
            "observation_gap": {"gedi_date_utc": "2025-04-25", "sentinel_date_utc": "2025-04-28",
                                "calendar_days": 3, "elapsed_hours_min": timedeltas.min().total_seconds() / 3600,
                                "elapsed_hours_max": timedeltas.max().total_seconds() / 3600},
            "reflectance": {"index_source": "existing 2025 index GeoTIFFs",
                            "scale_status": "provisional; Analytical TIFFs have no verified scale/offset",
                            "analytical_tiff_pixel_divisor_used": 65535,
                            "scl_scale_applied": False,
                            "note": "Source-product DN and Browser Analytical UINT16 values are not interchangeable"},
            "limits": ["SCL=4 means vegetation, not a confirmed forest boundary or tree species",
                       "SCL is an automated classification and may miss clouds or shadows",
                       "User product provenance cannot be independently proven from Analytical TIFF metadata"],
        }
        (OUTPUT / "scl_quality_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"SCL: {scl_path.name}; {scl_meta['crs']}, {scl_meta['resolution_m']} m, shape={scl_meta['shape']}")
        print(counts.to_string(index=False))
        print(f"Valid: {len(points)} -> {int(points.scl_valid.sum())}; cloud/shadow={int(points.scl_cloud_or_shadow.sum())}; boundary disagreements={int(points.scl_boundary_disagreement.sum())}")
        print(correlation.to_string(index=False))
        print(f"Saved: {OUTPUT}")
        return 0
    except Exception as exc:
        print(f"ERROR during {stage}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
