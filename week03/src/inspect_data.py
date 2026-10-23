"""Find Sentinel-2 GeoTIFFs and report their actual contents."""

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import rasterio


EXPECTED = ("B02", "B03", "B04", "B08", "B11", "B12")
BAND_PATTERN = re.compile(r"(?<![A-Za-z0-9])B(02|03|04|08|11|12)(?![A-Za-z0-9])", re.I)
DATE_PATTERN = re.compile(r"(?<!\d)(20\d{2}-\d{2}-\d{2})(?!\d)")


@dataclass(frozen=True)
class BandSource:
    path: Path
    band_index: int
    band: str
    date_label: str | None
    date_is_exact: bool


def _band_from_text(text: str | None) -> str | None:
    matches = {"B" + match for match in BAND_PATTERN.findall(text or "")}
    return next(iter(matches)) if len(matches) == 1 else None


def _date_label(dataset, path: Path) -> tuple[str | None, bool]:
    for tags in (dataset.tags(), dataset.tags(1)):
        for key, value in tags.items():
            if key.lower() in {"acquisition_date", "acquisition_time", "sensing_time", "datetime"}:
                dates = DATE_PATTERN.findall(value)
                if len(set(dates)) == 1:
                    return dates[0], True
    dates = list(dict.fromkeys(DATE_PATTERN.findall(path.name)))
    if len(dates) == 1:
        return dates[0], False
    if dates:
        return " to ".join(dates), False
    return None, False


def inspect_directory(folder: Path) -> tuple[dict[str, BandSource], list[str]]:
    files = sorted(p for p in folder.rglob("*") if p.is_file() and p.suffix.lower() in {".tif", ".tiff"})
    if not files:
        raise FileNotFoundError(f"No TIFF files in {folder}")
    print(f"Found {len(files)} TIFF file(s) in {folder}")
    sources: dict[str, BandSource] = {}
    duplicates: set[str] = set()
    warnings: list[str] = []
    for path in files:
        with rasterio.open(path) as ds:
            print(f"\nFile: {path.name}")
            print(f"  size={ds.width}x{ds.height}, count={ds.count}, dtypes={ds.dtypes}")
            print(f"  CRS={ds.crs}, resolution={ds.res}, bounds={tuple(ds.bounds)}")
            print(f"  NoData={ds.nodatavals}, scales={ds.scales}, offsets={ds.offsets}, units={ds.units}")
            print(f"  descriptions={ds.descriptions}, mask flags={ds.mask_flag_enums}")
            print(f"  dataset tags={ds.tags()}")
            filename_band = _band_from_text(path.stem)
            date_label, date_is_exact = _date_label(ds, path)
            data_mask_index = None
            if ds.count == 2:
                second_name = (ds.descriptions[1] or "").lower()
                second_tags = " ".join(f"{k}={v}" for k, v in ds.tags(2).items()).lower()
                second = ds.read(2)
                if "datamask" in second_name + second_tags or (filename_band and np.isin(second, [0, 1]).all()):
                    data_mask_index = 2
                    print("  band 2: dataMask (0=invalid, 1=valid); inferred from binary values" if "datamask" not in second_name + second_tags else "  band 2: dataMask (metadata)")
            for i in ds.indexes:
                tags = ds.tags(i)
                values = ds.read(i)
                valid = ds.read_masks(i) > 0
                nodata = ds.nodatavals[i - 1]
                if nodata is not None:
                    valid &= values != nodata
                if data_mask_index and i != data_mask_index:
                    valid &= ds.read(data_mask_index) != 0
                valid &= np.isfinite(values)
                band_name = _band_from_text(ds.descriptions[i - 1]) or _band_from_text(" ".join(tags.values()))
                if i == 1 and not band_name:
                    band_name = filename_band
                low = float(values[valid].min()) if valid.any() else None
                high = float(values[valid].max()) if valid.any() else None
                print(f"  band {i}: identified={band_name or ('dataMask' if i == data_mask_index else 'unknown')}, description={ds.descriptions[i - 1]!r}, tags={tags}, valid={int(valid.sum())}, range=[{low}, {high}]")
                if i == data_mask_index or band_name is None:
                    if i != data_mask_index:
                        warnings.append(f"Could not identify {path.name} band {i}; skipped")
                    continue
                if filename_band and i == 1 and band_name != filename_band:
                    warnings.append(f"Conflicting band identification in {path.name}; skipped")
                    continue
                if band_name in duplicates:
                    warnings.append(f"Additional duplicate {band_name}: {path.name}; skipped")
                    continue
                if band_name in sources:
                    warnings.append(f"Duplicate {band_name}: {sources[band_name].path.name} and {path.name}; both excluded")
                    sources.pop(band_name)
                    duplicates.add(band_name)
                    continue
                sources[band_name] = BandSource(path, i, band_name, date_label, date_is_exact)
    missing = sorted(set(EXPECTED) - set(sources))
    print(f"\nIdentified bands: {', '.join(sorted(sources)) or 'none'}")
    print(f"Missing bands: {', '.join(missing) or 'none'}")
    return sources, warnings
