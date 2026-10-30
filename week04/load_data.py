"""Load and summarize the Week03 SCL-validated vegetation candidates."""

from pathlib import Path

import pandas as pd


DATA_PATH = (
    Path(__file__).resolve().parent
    / ".."
    / "week03"
    / "outputs"
    / "sentinel_2025"
    / "scl_validation"
    / "gedi_scl_valid_points.csv"
)
STAT_COLUMNS = ["agbd", "ndvi", "evi", "ndmi"]
COORD_COLUMNS = ["latitude", "longitude", "easting_m", "northing_m"]


def load_data() -> pd.DataFrame:
    """Read the original Week03 CSV without copying or changing it."""
    if not DATA_PATH.is_file():
        raise FileNotFoundError(f"Week03 final CSV not found: {DATA_PATH.resolve()}")

    data = pd.read_csv(DATA_PATH)
    missing = set(STAT_COLUMNS + COORD_COLUMNS) - set(data.columns)
    if missing:
        raise ValueError(f"Required columns missing from Week03 CSV: {sorted(missing)}")
    return data


def main() -> None:
    data = load_data()
    print(f"Data file: {DATA_PATH.resolve()}")
    print(f"Total rows: {len(data)}")
    print("Columns:")
    print(", ".join(data.columns))
    print("\nAGBD and vegetation-index summary:")
    print(data[STAT_COLUMNS].describe().to_string())
    print("\nMissing values by column:")
    print(data.isna().sum().to_string())


if __name__ == "__main__":
    main()
