"""Reproducible, exploratory Random Forest AGBD baseline on 90 Week03 points."""

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "outputs" / "baseline_model"
os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split

from load_data import DATA_PATH, load_data


FEATURES = ("ndvi", "evi", "ndmi")
TARGET = "agbd"
SEED = 42
BLOCK_M = 1000.0
BUFFER_M = 500.0


def validate_data(data: pd.DataFrame) -> tuple[np.ndarray, dict]:
    needed = [*FEATURES, TARGET, "easting_m", "northing_m", "point_id",
              "gedi_file", "beam", "shot_number", "scl_valid"]
    missing_columns = sorted(set(needed) - set(data.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    if len(data) != 90:
        raise ValueError(f"Expected 90 SCL-valid observations; found {len(data)}")
    if not data.scl_valid.astype(str).str.lower().eq("true").all():
        raise ValueError("CSV includes a point that did not pass SCL validation")
    numeric = data[[*FEATURES, TARGET, "easting_m", "northing_m"]].apply(
        pd.to_numeric, errors="coerce")
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("Required feature, target, or coordinate contains missing/nonfinite values")
    if (numeric[TARGET] < 0).any():
        raise ValueError("Negative AGBD found")
    duplicate_ids = int(data.point_id.duplicated().sum())
    duplicate_shots = int(data.duplicated(["gedi_file", "beam", "shot_number"]).sum())
    duplicate_coords = int(numeric.duplicated(["easting_m", "northing_m"]).sum())
    if duplicate_ids or duplicate_shots or duplicate_coords:
        raise ValueError("Duplicate observation IDs, GEDI shots, or coordinates detected")

    coords = numeric[["easting_m", "northing_m"]].to_numpy()
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    np.fill_diagonal(distances, np.inf)
    nearest = distances.min(axis=1)
    outliers = {}
    for column in [TARGET, *FEATURES]:
        q1, q3 = numeric[column].quantile([.25, .75])
        iqr = q3 - q1
        low, high = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        outliers[column] = {"lower_fence": float(low), "upper_fence": float(high),
                            "flagged_n": int(((numeric[column] < low) | (numeric[column] > high)).sum())}
    diagnostics = {
        "source": str(DATA_PATH.resolve()), "observations": len(data),
        "missing_required": numeric.isna().sum().to_dict(),
        "duplicate_point_id": duplicate_ids, "duplicate_gedi_shot": duplicate_shots,
        "duplicate_coordinates": duplicate_coords,
        "gedi_files": data.gedi_file.value_counts().to_dict(),
        "coordinate_extent_utm52n_m": {
            "easting_min": float(coords[:, 0].min()), "easting_max": float(coords[:, 0].max()),
            "northing_min": float(coords[:, 1].min()), "northing_max": float(coords[:, 1].max())},
        "nearest_neighbor_m": {"min": float(nearest.min()),
                               "median": float(np.median(nearest)), "max": float(nearest.max())},
        "iqr_outlier_flags_not_removed": outliers,
        "numeric_summary": numeric[[TARGET, *FEATURES]].describe().to_dict(),
    }
    return coords, diagnostics


def spatial_split(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    # Easting-only blocks follow the long axis of the observed GEDI track.
    # The holdout is chosen using coordinates alone, before inspecting AGBD.
    blocks = np.floor((coords[:, 0] - coords[:, 0].min()) / BLOCK_M).astype(int)
    test = np.flatnonzero(blocks == blocks.max())
    candidates = np.flatnonzero(blocks != blocks.max())
    if len(test) < 10 or len(candidates) < 30:
        raise ValueError("Spatial block split has too few test or train candidates")
    distances = np.linalg.norm(coords[candidates, None, :] - coords[None, test, :], axis=2)
    keep = distances.min(axis=1) >= BUFFER_M
    train, excluded = candidates[keep], candidates[~keep]
    if len(train) < 30:
        raise ValueError("500 m buffer leaves too few training observations")
    minimum = float(np.linalg.norm(coords[train, None, :] - coords[None, test, :], axis=2).min())
    if minimum < BUFFER_M:
        raise AssertionError("Spatial train/test distance is below buffer")
    details = {"method": "easternmost 1 km easting block, 500 m test buffer",
               "block_origin_easting_m": float(coords[:, 0].min()),
               "block_width_m": BLOCK_M, "test_block_id": int(blocks.max()),
               "block_counts": {str(k): int((blocks == k).sum()) for k in np.unique(blocks)},
               "buffer_m": BUFFER_M, "buffer_excluded_n": len(excluded),
               "min_train_test_distance_m": minimum}
    return train, test, excluded, details


def evaluate(name: str, train: np.ndarray, test: np.ndarray, data: pd.DataFrame,
             coords: np.ndarray, excluded_n: int = 0) -> tuple[dict, pd.DataFrame, np.ndarray]:
    if np.intersect1d(train, test).size:
        raise AssertionError("Train/test observations overlap")
    if set(data.iloc[train].point_id) & set(data.iloc[test].point_id):
        raise AssertionError("Train/test point IDs overlap")
    model = RandomForestRegressor(n_estimators=100, random_state=SEED, n_jobs=1)
    model.fit(data.iloc[train][list(FEATURES)], data.iloc[train][TARGET])
    actual = data.iloc[test][TARGET].to_numpy(float)
    predicted = model.predict(data.iloc[test][list(FEATURES)])
    min_distance = float(np.linalg.norm(coords[train, None, :] - coords[None, test, :], axis=2).min())
    metrics = {"validation": name, "train_n": len(train), "test_n": len(test),
               "buffer_excluded_n": excluded_n,
               "min_train_test_distance_m": min_distance,
               "MAE": mean_absolute_error(actual, predicted),
               "RMSE": np.sqrt(mean_squared_error(actual, predicted)),
               "R2": r2_score(actual, predicted)}
    predictions = pd.DataFrame({"validation": name, "point_id": data.iloc[test].point_id.to_numpy(),
                                "actual_agbd": actual, "predicted_agbd": predicted,
                                "easting_m": coords[test, 0], "northing_m": coords[test, 1]})
    return metrics, predictions, model.feature_importances_


def figures(metrics: pd.DataFrame, predictions: pd.DataFrame,
            importances: dict[str, np.ndarray], data: pd.DataFrame,
            train_spatial: np.ndarray, test_spatial: np.ndarray,
            excluded: np.ndarray) -> None:
    colors = {"random_80_20": "#2563a6", "spatial_block_500m_buffer": "#bd5a35"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for ax, (name, group) in zip(axes, predictions.groupby("validation", sort=False)):
        ax.scatter(group.actual_agbd, group.predicted_agbd, alpha=.8, s=32, c=colors[name])
        bound = max(group.actual_agbd.max(), group.predicted_agbd.max()) * 1.05
        ax.plot([0, bound], [0, bound], "k--", lw=1)
        ax.set(xlim=(0, bound), ylim=(0, bound), title=f"{name} (n={len(group)})",
               xlabel="Observed GEDI AGBD (Mg/ha)", ylabel="Predicted AGBD (Mg/ha)")
        ax.grid(alpha=.2)
    fig.savefig(OUTPUT / "actual_vs_predicted.png", dpi=160)
    plt.close(fig)

    x = np.arange(len(FEATURES))
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for offset, (name, values) in zip((-.2, .2), importances.items()):
        ax.bar(x + offset, values, width=.4, label=name, color=colors[name])
    ax.set(xticks=x, xticklabels=[s.upper() for s in FEATURES], ylabel="Impurity-based importance",
           title="Training-set feature importance")
    ax.legend(fontsize=8)
    fig.savefig(OUTPUT / "feature_importance.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, metric in zip(axes, ("MAE", "RMSE", "R2")):
        ax.bar(range(len(metrics)), metrics[metric], color=[colors[n] for n in metrics.validation])
        ax.set(xticks=range(len(metrics)), xticklabels=["Random", "Spatial"],
               title=metric, ylabel="Mg/ha" if metric != "R2" else "R²")
        ax.axhline(0, color="black", lw=.7)
        for i, value in enumerate(metrics[metric]):
            ax.annotate(f"{value:.2f}", (i, value), xytext=(0, 4 if value >= 0 else -12),
                        textcoords="offset points", ha="center", fontsize=9)
    fig.savefig(OUTPUT / "validation_comparison.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for label, ids, color in (("Train", train_spatial, "#2563a6"),
                              ("Buffer excluded", excluded, "#888888"),
                              ("Validation", test_spatial, "#bd5a35")):
        ax.scatter(data.iloc[ids].easting_m, data.iloc[ids].northing_m,
                   label=f"{label} (n={len(ids)})", s=28, c=color)
    ax.set(xlabel="UTM 52N easting (m)", ylabel="UTM 52N northing (m)",
           title="Spatial holdout and buffer")
    ax.legend()
    ax.grid(alpha=.2)
    fig.savefig(OUTPUT / "spatial_split.png", dpi=160)
    plt.close(fig)


def main() -> None:
    data = load_data()
    coords, diagnostics = validate_data(data)
    all_indices = np.arange(len(data))
    random_train, random_test = train_test_split(all_indices, test_size=.2, random_state=SEED)
    spatial_train, spatial_test, excluded, spatial_details = spatial_split(coords)
    results = [evaluate("random_80_20", random_train, random_test, data, coords),
               evaluate("spatial_block_500m_buffer", spatial_train, spatial_test,
                        data, coords, len(excluded))]
    metrics = pd.DataFrame([item[0] for item in results])
    predictions = pd.concat([item[1] for item in results], ignore_index=True)
    importances = {name: item[2] for name, item in zip(metrics.validation, results)}
    diagnostics["splits"] = {"random_state": SEED, "random_test_fraction": .2,
                             "random_min_train_test_distance_m": float(metrics.iloc[0].min_train_test_distance_m),
                             "spatial": spatial_details,
                             "spatial_train_point_ids": data.iloc[spatial_train].point_id.tolist(),
                             "spatial_test_point_ids": data.iloc[spatial_test].point_id.tolist(),
                             "spatial_buffer_excluded_point_ids": data.iloc[excluded].point_id.tolist()}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUTPUT / "model_metrics.csv", index=False)
    predictions.to_csv(OUTPUT / "validation_predictions.csv", index=False)
    (OUTPUT / "data_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
    figures(metrics, predictions, importances, data, spatial_train, spatial_test, excluded)
    print(metrics.to_string(index=False))
    print(f"Saved outputs to {OUTPUT}")


if __name__ == "__main__":
    main()
