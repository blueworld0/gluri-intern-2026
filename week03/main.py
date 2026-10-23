"""Inspect and prepare Sentinel-2 inputs for later AGB modelling."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import rasterio

from src.inspect_data import inspect_directory
from src.preprocess import prepare_bands
from src.vegetation_indices import calculate
from src.visualize import save_index_map, save_ndvi_histogram, save_rgb


ROOT = Path(__file__).resolve().parent
FIELDS = ["index", "valid_pixels", "min", "max", "mean", "std", "median"]


def save_raster(path: Path, data: np.ndarray, profile: dict, name: str,
                scale: float | None, offset: float, acquisition_date: str | None) -> None:
    output_profile = profile.copy()
    output_profile.update(driver="GTiff", count=1, dtype="float32", nodata=-9999.0,
                          compress="deflate", predictor=3)
    output_profile.pop("photometric", None)
    with rasterio.open(path, "w", **output_profile) as dst:
        dst.write(np.where(np.isfinite(data), data, -9999.0).astype(np.float32), 1)
        dst.set_band_description(1, name.upper())
        dst.update_tags(1, units="unitless", formula=name.upper())
        dst.update_tags(reflectance_scale=str(scale) if scale is not None else "GeoTIFF metadata",
                        reflectance_offset=str(offset), acquisition_date=acquisition_date or "unknown")


def configured_calibration(sources: dict) -> tuple[float | None, float, str | None, list[str]]:
    config_path = ROOT / "reflectance_config.json"
    if not config_path.exists():
        return None, 0.0, None, ["No reflectance_config.json; using GeoTIFF metadata only"]
    config = json.loads(config_path.read_text(encoding="utf-8"))
    expected = config["sha256"]
    mismatches = []
    for band, source in sources.items():
        digest = hashlib.sha256(source.path.read_bytes()).hexdigest()
        if expected.get(band) != digest:
            mismatches.append(band)
    if mismatches:
        return None, 0.0, None, [f"Calibration config does not match file hash for {', '.join(mismatches)}; no configured scale applied"]
    print(f"Calibration: {config['source']}; SHA-256 matched for {', '.join(sorted(sources))}")
    return (float(config["reflectance_scale"]), float(config["reflectance_offset"]),
            config.get("acquisition_date"), [])


def diagnose_bands(prepared) -> None:
    print("\nReflectance checks (valid pixels only)")
    for band, values in prepared.arrays.items():
        valid = np.isfinite(values)
        sample = values[valid]
        raw = prepared.raw_arrays[band][valid]
        print(f"{band}: n={sample.size:,}, min={sample.min():.5f}, p01={np.percentile(sample, 1):.5f}, "
              f"median={np.median(sample):.5f}, p99={np.percentile(sample, 99):.5f}, max={sample.max():.5f}, "
              f"raw=65535: {int(np.count_nonzero(raw == 65535))}, reflectance outside [0,1]: "
              f"{int(np.count_nonzero((sample < 0) | (sample > 1)))}")


def diagnose_indices(indices: dict[str, np.ndarray], prepared) -> None:
    print("\nIndex plausibility checks")
    for name, data in indices.items():
        finite = data[np.isfinite(data)]
        outside = int(np.count_nonzero((finite < -1) | (finite > 1)))
        print(f"{name.upper()}: outside [-1,1]: {outside:,}; abs(value)>2: {int(np.count_nonzero(np.abs(finite)>2)):,}")
    if "evi" in indices:
        evi = indices["evi"]
        extreme = np.isfinite(evi) & (np.abs(evi) > 2)
        if extreme.any():
            for label, position in (("minimum", np.nanargmin(evi)), ("maximum", np.nanargmax(evi))):
                row, col = np.unravel_index(position, evi.shape)
                nir = float(prepared.arrays["B08"][row, col])
                red = float(prepared.arrays["B04"][row, col])
                blue = float(prepared.arrays["B02"][row, col])
                denominator = nir + 6 * red - 7.5 * blue + 1
                x, y = prepared.profile["transform"] * (col + 0.5, row + 0.5)
                print(f"EVI {label}: row={row}, col={col}, EPSG:32652=({x:.2f}, {y:.2f}), "
                      f"EVI={float(evi[row,col]):.4f}, B08={nir:.6f}, B04={red:.6f}, "
                      f"B02={blue:.6f}, denominator={denominator:.8f}")


def statistics(name: str, data: np.ndarray) -> dict:
    values = data[np.isfinite(data)]
    return {"index": name.upper(), "valid_pixels": int(values.size),
            "min": float(values.min()) if values.size else np.nan,
            "max": float(values.max()) if values.size else np.nan,
            "mean": float(values.mean()) if values.size else np.nan,
            "std": float(values.std()) if values.size else np.nan,
            "median": float(np.median(values)) if values.size else np.nan}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reflectance-scale", type=float,
                        help="Verified multiplier in reflectance = pixel * scale + offset; required when GeoTIFF metadata lacks calibration")
    parser.add_argument("--reflectance-offset", type=float, default=0.0,
                        help="Verified additive offset (default 0); applies with --reflectance-scale")
    parser.add_argument("--acquisition-date", help="Verified acquisition date YYYY-MM-DD for a manual scale override")
    args = parser.parse_args()
    if args.reflectance_scale is None and args.reflectance_offset != 0:
        parser.error("--reflectance-offset requires --reflectance-scale")
    output = ROOT / "outputs"
    raster_dir, figure_dir = output / "rasters", output / "figures"
    stage = "file inspection"
    try:
        sources, warnings = inspect_directory(ROOT / "Sentinel-2")
        config_scale, config_offset, config_date, config_warnings = configured_calibration(sources)
        warnings.extend(config_warnings)
        scale = args.reflectance_scale if args.reflectance_scale is not None else config_scale
        offset = args.reflectance_offset if args.reflectance_scale is not None else config_offset
        stage = "preprocessing"
        prepared = prepare_bands(sources, scale, offset)
        if args.acquisition_date:
            if any(src.date_label and args.acquisition_date not in src.date_label for src in sources.values()):
                raise ValueError("--acquisition-date conflicts with source date labels")
            prepared.date_label, prepared.date_is_exact = args.acquisition_date, True
        elif config_date and scale == config_scale:
            prepared.date_label, prepared.date_is_exact = config_date, True
        warnings.extend(w for w in prepared.warnings
                        if not (prepared.date_is_exact and w.startswith("Filename date is")))
        if prepared.arrays:
            diagnose_bands(prepared)
        stage = "vegetation indices"
        indices, skipped = calculate(prepared.arrays)
        diagnose_indices(indices, prepared)
        raster_dir.mkdir(parents=True, exist_ok=True)
        figure_dir.mkdir(parents=True, exist_ok=True)
        # Clear prior products so a later incomplete run cannot leave misleading index files.
        for name in ("ndvi", "evi", "ndmi"):
            for path in (raster_dir / f"{name}.tif", figure_dir / f"{name}.png"):
                path.unlink(missing_ok=True)
        for path in (figure_dir / "rgb.png", figure_dir / "ndvi_histogram.png"):
            path.unlink(missing_ok=True)
        stage = "raster and figure output"
        created: list[Path] = []
        rgb = save_rgb(prepared.raw_arrays, prepared.profile, figure_dir,
                       prepared.date_label, prepared.date_is_exact)
        if rgb:
            created.append(rgb)
        for name, data in indices.items():
            raster_path = raster_dir / f"{name}.tif"
            save_raster(raster_path, data, prepared.profile, name, scale, offset,
                        prepared.date_label if prepared.date_is_exact else None)
            created.append(raster_path)
            created.append(save_index_map(name, data, prepared.profile, figure_dir,
                                          prepared.date_label, prepared.date_is_exact))
        if "ndvi" in indices:
            created.append(save_ndvi_histogram(indices["ndvi"], figure_dir,
                                               prepared.date_label, prepared.date_is_exact))
        stage = "statistics"
        rows = [statistics(name, data) for name, data in indices.items()]
        csv_path = output / "statistics.csv"
        pd.DataFrame(rows, columns=FIELDS).to_csv(csv_path, index=False)
        created.append(csv_path)
        print("\nSummary")
        for warning in warnings:
            print(f"WARNING: {warning}")
        for name, missing in skipped.items():
            print(f"SKIPPED {name.upper()}: missing calibrated reflectance for {', '.join(missing)}")
        for row in rows:
            print(f"{row['index']}: valid={row['valid_pixels']:,}, min={row['min']:.4f}, max={row['max']:.4f}, mean={row['mean']:.4f}, std={row['std']:.4f}, median={row['median']:.4f}")
        print("Created: " + (", ".join(str(p.relative_to(ROOT)) for p in created) or "none"))
        return 0
    except Exception as exc:
        print(f"ERROR during {stage}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
