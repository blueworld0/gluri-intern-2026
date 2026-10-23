"""Save georeferenced map figures and a histogram."""

from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rasterio.transform import array_bounds


def _extent(profile: dict) -> tuple[float, float, float, float]:
    west, south, east, north = array_bounds(profile["height"], profile["width"], profile["transform"])
    return west, east, south, north


def _title(name: str, date_label: str | None, exact: bool) -> str:
    return f"{name}\nAcquisition: {date_label}" if exact and date_label else name


def save_rgb(raw_arrays: dict[str, np.ndarray], profile: dict, folder: Path,
             date_label: str | None, exact: bool) -> Path | None:
    if not all(band in raw_arrays for band in ("B04", "B03", "B02")):
        return None
    rgb = np.stack([raw_arrays[b] for b in ("B04", "B03", "B02")], axis=-1)
    valid = np.all(np.isfinite(rgb), axis=-1)
    if not valid.any():
        return None
    display = np.zeros_like(rgb, dtype=np.float32)
    for channel in range(3):
        low, high = np.percentile(rgb[:, :, channel][valid], [2, 98])
        if high > low:
            display[:, :, channel] = np.clip((rgb[:, :, channel] - low) / (high - low), 0, 1)
    display[~valid] = 0
    rgba = np.dstack([display, valid.astype(np.float32)])
    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    ax.set_facecolor("0.8")
    ax.imshow(rgba, extent=_extent(profile), origin="upper")
    ax.set_title(_title("Sentinel-2 RGB (2–98% raw-value stretch)", date_label, exact))
    ax.set_xlabel(f"Easting ({profile['crs']})")
    ax.set_ylabel("Northing (m)")
    # RGB has three independently stretched channels, so a single numeric colorbar is not meaningful.
    ax.text(0.01, 0.01, "R=B04, G=B03, B=B02 | NoData=gray | per-channel range: 0–1 display",
            transform=ax.transAxes, fontsize=8, color="white", bbox={"facecolor": "black", "alpha": 0.6})
    path = folder / "rgb.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def save_index_map(name: str, data: np.ndarray, profile: dict, folder: Path,
                   date_label: str | None, exact: bool) -> Path:
    valid = np.isfinite(data)
    masked = np.ma.masked_invalid(data)
    cmap = plt.get_cmap("RdYlGn" if name == "ndvi" else "BrBG" if name == "ndmi" else "viridis").copy()
    cmap.set_bad("0.8")
    actual_min = float(data[valid].min()) if valid.any() else -1.0
    actual_max = float(data[valid].max()) if valid.any() else 1.0
    vmin, vmax = (float(v) for v in np.percentile(data[valid], [2, 98])) if valid.any() else (-1.0, 1.0)
    if vmin == vmax:
        vmin, vmax = vmin - 0.01, vmax + 0.01
    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    im = ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=_extent(profile), origin="upper")
    fig.colorbar(im, ax=ax, label=name.upper())
    ax.set_title(_title(f"{name.upper()} | data range {actual_min:.3f} to {actual_max:.3f}\nNoData=gray", date_label, exact), fontsize=11)
    ax.set_xlabel(f"Easting ({profile['crs']})")
    ax.set_ylabel("Northing (m)")
    path = folder / f"{name}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def save_ndvi_histogram(data: np.ndarray, folder: Path, date_label: str | None,
                        exact: bool) -> Path:
    values = data[np.isfinite(data)]
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    if values.size:
        ax.hist(values, bins=60, color="forestgreen")
    ax.set_title(_title("NDVI distribution", date_label, exact))
    ax.set_xlabel("NDVI")
    ax.set_ylabel("Valid pixel count")
    ax.grid(alpha=0.25)
    path = folder / "ndvi_histogram.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path
