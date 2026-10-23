"""Read, mask, align, and calibrate Sentinel-2 bands."""

from dataclasses import dataclass
import math

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

from .inspect_data import BandSource


@dataclass
class Prepared:
    arrays: dict[str, np.ndarray]
    raw_arrays: dict[str, np.ndarray]
    profile: dict
    date_label: str | None
    date_is_exact: bool
    warnings: list[str]


def _conversion(ds, index: int, user_scale: float | None, user_offset: float) -> tuple[float, float] | None:
    if user_scale is not None:
        return user_scale, user_offset
    scale, offset = ds.scales[index - 1], ds.offsets[index - 1]
    unit = (ds.units[index - 1] or "").lower()
    tags = {k.lower(): v for k, v in (ds.tags() | ds.tags(index)).items()}
    if "normalizationfactor" in tags:
        return float(tags["normalizationfactor"]), 0.0
    if "reflectance" in unit or tags.get("units", "").lower() == "reflectance":
        return scale, offset
    if not math.isclose(scale, 1.0) or not math.isclose(offset, 0.0):
        return scale, offset
    return None


def _mask_index(ds, source: BandSource) -> int | None:
    if ds.count != 2 or source.band_index != 1:
        return None
    description = (ds.descriptions[1] or "").lower()
    tags = " ".join(f"{k}={v}" for k, v in ds.tags(2).items()).lower()
    if "datamask" in description + tags:
        return 2
    return 2 if np.isin(ds.read(2), [0, 1]).all() else None


def prepare_bands(sources: dict[str, BandSource], user_scale: float | None = None,
                  user_offset: float = 0.0) -> Prepared:
    if user_scale is not None and (not math.isfinite(user_scale) or user_scale <= 0):
        raise ValueError("--reflectance-scale must be finite and positive")
    if not math.isfinite(user_offset):
        raise ValueError("--reflectance-offset must be finite")
    if not sources:
        raise ValueError("No identifiable spectral bands")
    warnings: list[str] = []
    dates = {src.date_label for src in sources.values()}
    if len(dates) > 1:
        raise ValueError(f"Different or unknown date labels across bands: {dates}; refusing to combine")
    date_label = next(iter(dates))
    date_is_exact = all(src.date_is_exact for src in sources.values())
    if not date_is_exact:
        warnings.append("Filename date is a request interval, not verified acquisition time")
    reference = sources.get("B04") or next(iter(sources.values()))
    with rasterio.open(reference.path) as ref:
        if ref.crs is None:
            raise ValueError(f"Missing CRS: {reference.path.name}")
        profile = ref.profile.copy()
        ref_crs, ref_transform = ref.crs, ref.transform
        shape = (ref.height, ref.width)
        if ref.crs.to_epsg() != 32652:
            warnings.append(f"Reference CRS is {ref.crs}, not expected EPSG:32652")
        if any(abs(x - 10) > 0.1 for x in ref.res):
            warnings.append(f"Reference resolution {ref.res} differs from target 10 m")
    arrays: dict[str, np.ndarray] = {}
    raw_arrays: dict[str, np.ndarray] = {}
    for band, source in sources.items():
        with rasterio.open(source.path) as ds:
            if ds.crs is None:
                warnings.append(f"{band}: missing CRS; skipped")
                continue
            original = ds.read(source.band_index).astype(np.float32)
            valid = (ds.read_masks(source.band_index) > 0) & np.isfinite(original)
            nodata = ds.nodatavals[source.band_index - 1]
            if nodata is not None:
                valid &= original != nodata
            mask_index = _mask_index(ds, source)
            if mask_index:
                valid &= ds.read(mask_index) > 0
            original[~valid] = np.nan
            same_grid = ds.crs == ref_crs and ds.transform.almost_equals(ref_transform) and original.shape == shape
            if not same_grid:
                warnings.append(f"{band}: CRS/size/resolution/extent differs; reprojected to {reference.band} grid")
                aligned = np.full(shape, np.nan, dtype=np.float32)
                aligned_mask = np.zeros(shape, dtype=np.uint8)
                reproject(original, aligned, src_transform=ds.transform, src_crs=ds.crs,
                          dst_transform=ref_transform, dst_crs=ref_crs, src_nodata=np.nan,
                          dst_nodata=np.nan, resampling=Resampling.bilinear)
                reproject(valid.astype(np.uint8), aligned_mask, src_transform=ds.transform,
                          src_crs=ds.crs, dst_transform=ref_transform, dst_crs=ref_crs,
                          src_nodata=0, dst_nodata=0, resampling=Resampling.nearest)
                aligned[aligned_mask == 0] = np.nan
                original = aligned
            raw_arrays[band] = original
            conversion = _conversion(ds, source.band_index, user_scale, user_offset)
            if conversion is None:
                warnings.append(f"{band}: no reliable reflectance scale/offset in GeoTIFF; index calculation disabled for this band")
                continue
            scale, offset = conversion
            reflectance = original * scale + offset
            reflectance[~np.isfinite(reflectance)] = np.nan
            arrays[band] = reflectance
            print(f"{band}: reflectance = pixel * {scale:g} + {offset:g}; valid={int(np.isfinite(reflectance).sum())}")
    return Prepared(arrays, raw_arrays, profile, date_label, date_is_exact, warnings)
