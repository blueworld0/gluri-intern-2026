"""Compute vegetation indices only from calibrated reflectance."""

import numpy as np


NEEDED = {"ndvi": ("B08", "B04"), "evi": ("B08", "B04", "B02"),
          "ndmi": ("B08", "B11")}


def calculate(arrays: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, list[str]]]:
    results: dict[str, np.ndarray] = {}
    skipped: dict[str, list[str]] = {}
    for name, needed in NEEDED.items():
        missing = [band for band in needed if band not in arrays]
        if missing:
            skipped[name] = missing
            continue
        nir = arrays["B08"].astype(np.float64)
        red = arrays["B04"].astype(np.float64)
        if name == "ndvi":
            numerator, denominator = nir - red, nir + red
        elif name == "evi":
            blue = arrays["B02"].astype(np.float64)
            numerator, denominator = 2.5 * (nir - red), nir + 6 * red - 7.5 * blue + 1
        else:
            swir = arrays["B11"].astype(np.float64)
            numerator, denominator = nir - swir, nir + swir
        valid = np.isfinite(numerator) & np.isfinite(denominator) & (np.abs(denominator) > 1e-12)
        out = np.full(nir.shape, np.nan, dtype=np.float32)
        np.divide(numerator, denominator, out=out, where=valid)
        out[~np.isfinite(out)] = np.nan
        results[name] = out
    return results, skipped
