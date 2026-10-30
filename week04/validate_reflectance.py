"""Audit Week03 TIFF radiometry and compare *hypothetical* scaling scenarios.

Neither scale is asserted to be the true conversion for these exports.
"""

import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "outputs/reflectance_validation"
os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio


WEEK03 = ROOT.parent / "week03"
SOURCE = WEEK03 / "Sentinel-2"
POINTS = WEEK03 / "outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv"
BANDS = ("B02", "B03", "B04", "B08", "B11", "B12")
SCALES = {"browser_65535_assumption": 65535.0, "dn_10000_sensitivity": 10000.0}
REFERENCES = {
    "browser_staff": "https://forum.dataspace.copernicus.eu/t/issues-with-download-image-analytical/3566",
    "evalscript": "https://docs.sentinel-hub.com/api/latest/evalscript/v3/",
    "l2a_examples": "https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Process/Examples/S2L2A.html",
}


def source_files() -> dict[tuple[int, str], Path]:
    result = {}
    for year, folder in ((2025, SOURCE / "2025-04-28"), (2026, SOURCE)):
        for band in BANDS:
            matches = [p for p in folder.glob("*.tif*")
                       if re.search(rf"_{band}_\(Raw\)\.(tif|tiff)$", p.name, re.I)]
            if len(matches) != 1:
                raise ValueError(f"Expected one {year} {band} TIFF; found {len(matches)}")
            result[year, band] = matches[0]
    return result


def audit_file(year: int, band: str, path: Path) -> tuple[dict, dict]:
    with rasterio.open(path) as ds:
        if ds.count != 2 or ds.dtypes != ("uint16", "uint16"):
            raise ValueError(f"Unexpected TIFF band structure: {path}")
        raw = ds.read(1)
        mask = ds.read(2)
        if not np.isin(mask, [0, 1]).all():
            raise ValueError(f"Nonbinary second band in {path}")
        valid = (mask == 1) & (ds.read_masks(1) > 0)
        if ds.nodatavals[0] is not None:
            valid &= raw != ds.nodatavals[0]
        values = raw[valid]
        if not len(values):
            raise ValueError(f"No valid pixels in {path}")
        metadata = {
            "year": year, "band": band, "file": str(path.relative_to(ROOT.parent)),
            "driver": ds.driver, "shape": list(ds.shape), "crs": str(ds.crs),
            "resolution": list(ds.res), "transform": list(ds.transform)[:6],
            "dtypes": list(ds.dtypes), "nodata": list(ds.nodatavals),
            "scales": list(ds.scales), "offsets": list(ds.offsets),
            "units": list(ds.units), "descriptions": list(ds.descriptions),
            "dataset_tags": ds.tags(), "band_tags": [ds.tags(1), ds.tags(2)],
            "mask_flags": [[str(flag) for flag in flags] for flags in ds.mask_flag_enums],
        }
        stats = {
            "year": year, "band": band, "total_pixels": raw.size,
            "valid_pixels": int(valid.sum()), "data_mask_zero": int((mask == 0).sum()),
            "valid_zero": int((values == 0).sum()),
            "valid_65535": int((values == 65535).sum()),
            "min": int(values.min()), "p25": float(np.percentile(values, 25)),
            "median": float(np.median(values)), "mean": float(values.mean()),
            "p75": float(np.percentile(values, 75)),
            "p99": float(np.percentile(values, 99)), "max": int(values.max()),
        }
        return metadata, stats


def index_values(bands: dict[str, np.ndarray], divisor: float) -> dict[str, np.ndarray]:
    blue, red, nir, swir = (bands[name] / divisor for name in ("B02", "B04", "B08", "B11"))
    with np.errstate(divide="ignore", invalid="ignore"):
        return {
            "ndvi": (nir - red) / (nir + red),
            "evi": 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1),
            "ndmi": (nir - swir) / (nir + swir),
        }


def point_comparison(files: dict[tuple[int, str], Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    points = pd.read_csv(POINTS)
    required = {"point_id", "raster_row", "raster_col", "ndvi", "evi", "ndmi"}
    if not required <= set(points):
        raise ValueError(f"Missing point columns: {sorted(required - set(points))}")
    if len(points) != 90 or points.point_id.duplicated().any():
        raise ValueError("Expected 90 unique SCL-valid candidates")
    rows = points.raster_row.to_numpy(dtype=int)
    cols = points.raster_col.to_numpy(dtype=int)
    raw = {}
    grid = None
    for band in BANDS:
        with rasterio.open(files[2025, band]) as ds:
            current = (ds.crs, ds.transform, ds.shape)
            if grid is None:
                grid = current
            elif current != grid:
                raise ValueError("2025 TIFF grids differ; direct sampling is invalid")
            if ((rows < 0) | (rows >= ds.height) | (cols < 0) | (cols >= ds.width)).any():
                raise ValueError("A point falls outside the TIFF grid")
            values = ds.read(1)[rows, cols].astype(float)
            valid = (ds.read(2)[rows, cols] == 1) & (ds.read_masks(1)[rows, cols] > 0)
            if not valid.all():
                raise ValueError(f"A candidate is masked in {band}")
            raw[band] = values

    result = points[["point_id", "raster_row", "raster_col", "agbd", "ndvi", "evi", "ndmi"]].copy()
    summary = []
    for scenario, divisor in SCALES.items():
        indices = index_values(raw, divisor)
        for name, values in indices.items():
            if not np.isfinite(values).all():
                raise ValueError(f"Nonfinite {name} for {scenario}")
            result[f"{name}_{scenario}"] = values
            summary.append({"scenario": scenario, "index": name.upper(), "n": len(values),
                            "min": values.min(), "median": np.median(values),
                            "mean": values.mean(), "max": values.max()})
    for name in ("ndvi", "evi", "ndmi"):
        delta = result[f"{name}_browser_65535_assumption"] - result[name]
        result[f"{name}_reproduction_error"] = delta
        result[f"{name}_scenario_difference"] = (
            result[f"{name}_dn_10000_sensitivity"] - result[f"{name}_browser_65535_assumption"])
        summary.append({"scenario": "week03_saved", "index": name.upper(), "n": len(result),
                        "min": result[name].min(), "median": result[name].median(),
                        "mean": result[name].mean(), "max": result[name].max()})
    return result, pd.DataFrame(summary)


def make_plot(points: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for ax, name in zip(axes, ("ndvi", "evi", "ndmi")):
        x = points[f"{name}_browser_65535_assumption"]
        y = points[f"{name}_dn_10000_sensitivity"]
        ax.scatter(x, y, s=17, alpha=.75)
        low, high = min(x.min(), y.min()), max(x.max(), y.max())
        ax.plot([low, high], [low, high], "k--", lw=1)
        ax.set(title=name.upper(), xlabel="TIFF / 65535 (Week03 assumption)",
               ylabel="TIFF / 10000 (sensitivity only)")
        ax.grid(alpha=.2)
    fig.savefig(OUTPUT / "index_scale_sensitivity.png", dpi=160)
    plt.close(fig)


def main() -> None:
    files = source_files()
    metadata, stats = [], []
    for (year, band), path in sorted(files.items()):
        item, row = audit_file(year, band, path)
        metadata.append(item)
        stats.append(row)
    points, summary = point_comparison(files)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "tiff_metadata.json").write_text(json.dumps({
        "sources": REFERENCES, "interpretation": "No per-file radiometric provenance or scale tag; both conversions are assumptions.",
        "files": metadata}, indent=2), encoding="utf-8")
    pd.DataFrame(stats).to_csv(OUTPUT / "band_statistics.csv", index=False)
    points.to_csv(OUTPUT / "candidate_index_comparison.csv", index=False)
    summary.to_csv(OUTPUT / "index_comparison_statistics.csv", index=False)
    make_plot(points)
    print(f"Audited {len(metadata)} TIFFs; compared {len(points)} SCL-valid 2025 candidates")
    print(f"Output: {OUTPUT}")
    print("Maximum absolute difference from Week03 saved indices:")
    for name in ("ndvi", "evi", "ndmi"):
        print(f"  {name.upper()}: {points[f'{name}_reproduction_error'].abs().max():.8g}")
    print("Mean scenario differences (TIFF / 10000 minus TIFF / 65535):")
    for name in ("ndvi", "evi", "ndmi"):
        print(f"  {name.upper()}: {points[f'{name}_scenario_difference'].mean():.8g}")
    print("Radiometric conversion remains unverified for these specific TIFFs.")


if __name__ == "__main__":
    main()
