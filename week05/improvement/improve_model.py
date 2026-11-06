"""Exploratory Week05 A/B/C experiments on unchanged Week04 holdouts."""

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WEEK04 = ROOT / "week04"
SOURCE = ROOT / "week03/Sentinel-2/2025-04-28"
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
sys.path.insert(0, str(WEEK04))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold

from compare_baselines import splits
from train_baseline import FEATURES, SEED, TARGET, load_data, validate_data


VALIDATIONS = ("random_80_20", "spatial_block_500m_buffer")
BASE = list(FEATURES)
FEATURE_SETS = {
    "base_plus_nbr": BASE + ["nbr_raw_ratio"],
    "base_plus_unused_raw_bands": BASE + ["raw_B03", "raw_B12"],
    "base_plus_nbr_and_raw": BASE + ["nbr_raw_ratio", "raw_B03", "raw_B12"],
}
RF_BASE = {"n_estimators": 100, "random_state": SEED, "n_jobs": 1}
COMPLEXITY = {
    "depth4_leaf2_all": {"max_depth": 4, "min_samples_leaf": 2, "max_features": 1.0},
    "depth6_leaf3_two_thirds": {"max_depth": 6, "min_samples_leaf": 3, "max_features": 0.67},
    "depth3_leaf4_all": {"max_depth": 3, "min_samples_leaf": 4, "max_features": 1.0},
}
QUALITY_THRESHOLD = 0.95
MIN_QUALITY_TRAIN_N = 30


def add_existing_tiff_features(data):
    """Sample available 2025 raw TIFFs at already verified Week03 row/column."""
    result = data.copy()
    rr, cc = data.raster_row.to_numpy(int), data.raster_col.to_numpy(int)
    raw = {}
    grid = None
    for band in ("B03", "B08", "B12"):
        files = list(SOURCE.glob(f"*_{band}_(Raw).tiff"))
        if len(files) != 1:
            raise FileNotFoundError(f"Required existing TIFF absent: {band}")
        with rasterio.open(files[0]) as ds:
            current = (ds.crs, ds.transform, ds.shape)
            if grid is None:
                grid = current
            elif current != grid:
                raise ValueError("Candidate band grids differ")
            if ds.count != 2 or ds.dtypes[0] != "uint16":
                raise ValueError(f"Unexpected TIFF bands: {band}")
            values = ds.read(1)[rr, cc].astype(float)
            valid = (ds.read(2)[rr, cc] == 1) & (ds.read_masks(1)[rr, cc] > 0)
            if ds.nodatavals[0] is not None:
                valid &= values != ds.nodatavals[0]
            if not valid.all():
                raise ValueError(f"Invalid existing TIFF point sample: {band}")
            raw[band] = values
    denominator = raw["B08"] + raw["B12"]
    if (denominator == 0).any():
        raise ValueError("NBR denominator zero at a Week04 point")
    result["nbr_raw_ratio"] = (raw["B08"] - raw["B12"]) / denominator
    result["raw_B03"] = raw["B03"]
    result["raw_B12"] = raw["B12"]
    if not np.isfinite(result[["nbr_raw_ratio", "raw_B03", "raw_B12"]].to_numpy()).all():
        raise ValueError("Nonfinite candidate feature")
    return result


def metrics(actual, predicted):
    return {"R2": float(r2_score(actual, predicted)),
            "RMSE": float(np.sqrt(mean_squared_error(actual, predicted))),
            "MAE": float(mean_absolute_error(actual, predicted))}


def inner_cv(data, train_indices, coords, features, rf_params):
    """Pooled out-of-fold RMSE using only outer training points and 500 m groups."""
    x = data.iloc[train_indices][features]
    y = data.iloc[train_indices][TARGET].to_numpy(float)
    groups = np.floor((coords[train_indices, 0] - coords[:, 0].min()) / 500).astype(int)
    if len(np.unique(groups)) < 3:
        raise ValueError("Fewer than three train-only spatial groups")
    oof = np.full(len(train_indices), np.nan)
    folds = []
    for fold, (fit, val) in enumerate(GroupKFold(n_splits=3).split(x, y, groups)):
        if len(fit) < 15:
            raise ValueError("Inner spatial fold leaves fewer than 15 training points")
        model = RandomForestRegressor(**RF_BASE, **rf_params)
        model.fit(x.iloc[fit], y[fit])
        oof[val] = model.predict(x.iloc[val])
        folds.append({"fold": fold + 1, "train_n": len(fit), "validation_n": len(val),
                      "validation_group_ids": sorted(np.unique(groups[val]).tolist())})
    if not np.isfinite(oof).all():
        raise ValueError("Incomplete inner CV prediction")
    return metrics(y, oof), folds


def score_holdout(data, train, test, features, rf_params, label, validation,
                  experiment, excluded_n, cv_rmse=np.nan, quality_condition="all"):
    model = (DummyRegressor(strategy="mean") if label == "DummyMean" else
             RandomForestRegressor(**RF_BASE, **rf_params))
    model.fit(data.iloc[train][features], data.iloc[train][TARGET])
    actual = data.iloc[test][TARGET].to_numpy(float)
    predicted = model.predict(data.iloc[test][features])
    row = {"validation": validation, "experiment": experiment, "model": label,
           "features": ",".join(features), "train_n": len(train), "test_n": len(test),
           "buffer_excluded_n": excluded_n, "quality_condition": quality_condition,
           "train_agbd_mean": float(data.iloc[train][TARGET].mean()),
           "test_agbd_mean": float(actual.mean()),
           "inner_cv_rmse": cv_rmse,
           "max_depth": rf_params.get("max_depth"),
           "min_samples_leaf": rf_params.get("min_samples_leaf", 1),
           "max_features": rf_params.get("max_features", 1.0),
           **metrics(actual, predicted)}
    predictions = pd.DataFrame({"validation": validation, "experiment": experiment,
                                "model": label, "point_id": data.iloc[test].point_id.to_numpy(),
                                "actual_agbd": actual, "predicted_agbd": predicted})
    return row, predictions


def validate_week04(data, split_map, results, predictions):
    existing_metrics = pd.read_csv(WEEK04 / "outputs/baseline_comparison/comparison_metrics.csv")
    existing_predictions = pd.read_csv(WEEK04 / "outputs/baseline_model/validation_predictions.csv")
    results = pd.DataFrame(results)
    predictions = pd.concat(predictions, ignore_index=True)
    for validation in VALIDATIONS:
        train, test, excluded = split_map[validation]
        for model in ("RandomForest", "DummyMean"):
            new = results.loc[(results.validation == validation) &
                              (results.experiment == "week04_baseline") &
                              (results.model == model)].iloc[0]
            old = existing_metrics.loc[(existing_metrics.validation == validation) &
                                       (existing_metrics.model == model)].iloc[0]
            if (new.train_n, new.test_n, new.buffer_excluded_n) != (old.train_n, old.test_n, old.buffer_excluded_n):
                raise ValueError("Week04 split counts differ")
            if not np.allclose([new.R2, new.RMSE, new.MAE], [old.R2, old.RMSE, old.MAE], atol=1e-8, rtol=0):
                raise ValueError(f"Week04 metrics not reproduced: {validation}, {model}")
        rf = predictions.loc[(predictions.validation == validation) &
                             (predictions.experiment == "week04_baseline") &
                             (predictions.model == "RandomForest")]
        old = existing_predictions.loc[existing_predictions.validation == validation]
        if not np.array_equal(rf.point_id.to_numpy(), old.point_id.to_numpy()):
            raise ValueError("Week04 RF test IDs/order differ")
        if not np.allclose(rf.predicted_agbd, old.predicted_agbd, atol=1e-8, rtol=0):
            raise ValueError("Week04 RF predictions differ")
    return results, predictions


def quality_summary(data, split_map):
    rows = []
    for condition, mask in (("all", np.ones(len(data), dtype=bool)),
                            ("sensitivity_ge_0.95", data.sensitivity.to_numpy(float) >= QUALITY_THRESHOLD)):
        subset = data.iloc[np.flatnonzero(mask)]
        rows.append({"validation": "overall", "partition": "all", "condition": condition,
                     "n": len(subset), "agbd_min": subset.agbd.min(),
                     "agbd_mean": subset.agbd.mean(), "agbd_median": subset.agbd.median(),
                     "agbd_max": subset.agbd.max(),
                     "sensitivity_min": subset.sensitivity.min(),
                     "sensitivity_median": subset.sensitivity.median(),
                     "sensitivity_max": subset.sensitivity.max()})
    for validation, (train, test, excluded) in split_map.items():
        for partition, indices in (("train", train), ("test", test), ("buffer_excluded", excluded)):
            for condition, mask in (("all", np.ones(len(indices), dtype=bool)),
                                    ("sensitivity_ge_0.95", data.iloc[indices].sensitivity.to_numpy(float) >= QUALITY_THRESHOLD)):
                subset = data.iloc[indices[np.asarray(mask)]]
                rows.append({"validation": validation, "partition": partition,
                             "condition": condition, "n": len(subset),
                             "agbd_min": subset.agbd.min(), "agbd_mean": subset.agbd.mean(),
                             "agbd_median": subset.agbd.median(), "agbd_max": subset.agbd.max(),
                             "sensitivity_min": subset.sensitivity.min(),
                             "sensitivity_median": subset.sensitivity.median(),
                             "sensitivity_max": subset.sensitivity.max()})
    return pd.DataFrame(rows)


def quality_flags(data):
    columns = ("l2_algrunflag", "l2a_quality_flag_rel3", "l4a_quality_flag_rel3",
               "degrade_flag", "degrade_include_flag", "elev_highestreturn_outlier_flag",
               "predictor_limit_flag", "response_limit_flag", "quality_valid", "scl_valid")
    rows = []
    for column in columns:
        for value, section in data.groupby(column, dropna=False):
            rows.append({"column": column, "value": str(value), "n": len(section),
                         "agbd_mean": section.agbd.mean(), "agbd_median": section.agbd.median()})
    return pd.DataFrame(rows)


def plot_results(results):
    shown = results.loc[(results.experiment.isin(["week04_baseline", "feature_selected_inner_cv",
                                                  "complexity_selected_inner_cv", "quality_train_sens_ge_095"])) &
                        (results.model == "RandomForest")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    labels = {"week04_baseline": "Week04 RF", "feature_selected_inner_cv": "A: features",
              "complexity_selected_inner_cv": "B: simpler RF",
              "quality_train_sens_ge_095": "C: quality train"}
    for ax, validation, title in zip(axes, VALIDATIONS, ("Random holdout", "Spatial holdout")):
        section = shown.loc[shown.validation.eq(validation)]
        dummy_rmse = results.loc[(results.validation == validation) &
                                 (results.experiment == "week04_baseline") &
                                 (results.model == "DummyMean"), "RMSE"].iloc[0]
        names = [labels[name] for name in section.experiment]
        ax.bar(np.arange(len(section)), section.RMSE, color="#3b7198")
        ax.axhline(dummy_rmse, color="#ba5935", linestyle="--", linewidth=1.5,
                   label=f"Week04 Dummy {dummy_rmse:.1f}")
        ax.set(xticks=np.arange(len(section)), xticklabels=names, ylabel="RMSE (Mg/ha)", title=title)
        ax.tick_params(axis="x", labelrotation=25)
        for i, value in enumerate(section.RMSE):
            ax.text(i, value + 1, f"{value:.1f}", ha="center", fontsize=8)
        ax.grid(axis="y", alpha=.2)
        ax.legend(fontsize=8)
    fig.savefig(HERE / "performance_comparison.png", dpi=160)
    plt.close(fig)


def main():
    data = load_data()
    coords, _ = validate_data(data)
    data = add_existing_tiff_features(data)
    split_map = splits(data, coords)
    results, prediction_tables, cv_rows, cv_details = [], [], [], {}
    for validation, (train, test, excluded) in split_map.items():
        if np.intersect1d(train, test).size or np.intersect1d(excluded, train).size or np.intersect1d(excluded, test).size:
            raise ValueError("Outer split overlaps")
        # Baselines: exact Week04 settings, assessed before any new selection.
        for model in ("RandomForest", "DummyMean"):
            row, pred = score_holdout(data, train, test, BASE, {}, model, validation,
                                      "week04_baseline", len(excluded))
            results.append(row)
            prediction_tables.append(pred)
        base_cv, base_folds = inner_cv(data, train, coords, BASE, {})
        cv_rows.append({"validation": validation, "experiment": "baseline_reference",
                        "candidate": "week04_rf", "features": ",".join(BASE), **base_cv})
        cv_details[f"{validation}/baseline_reference"] = base_folds
        # A: choose one additional-input set solely from outer training CV.
        feature_cv = []
        for name, features in FEATURE_SETS.items():
            cv, folds = inner_cv(data, train, coords, features, {})
            feature_cv.append((name, cv["RMSE"]))
            cv_rows.append({"validation": validation, "experiment": "A_features", "candidate": name,
                            "features": ",".join(features), **cv})
            cv_details[f"{validation}/A/{name}"] = folds
        feature_name = min(feature_cv, key=lambda item: item[1])[0]
        feature_rmse = dict(feature_cv)[feature_name]
        row, pred = score_holdout(data, train, test, FEATURE_SETS[feature_name], {},
                                  "RandomForest", validation, "feature_selected_inner_cv",
                                  len(excluded), feature_rmse)
        row["selected_candidate"] = feature_name
        results.append(row)
        prediction_tables.append(pred)
        # B: a deliberately small, lower-complexity grid on original features.
        complexity_cv = []
        for name, params in COMPLEXITY.items():
            cv, folds = inner_cv(data, train, coords, BASE, params)
            complexity_cv.append((name, cv["RMSE"]))
            cv_rows.append({"validation": validation, "experiment": "B_complexity", "candidate": name,
                            "features": ",".join(BASE), **cv})
            cv_details[f"{validation}/B/{name}"] = folds
        complexity_name = min(complexity_cv, key=lambda item: item[1])[0]
        complexity_rmse = dict(complexity_cv)[complexity_name]
        row, pred = score_holdout(data, train, test, BASE, COMPLEXITY[complexity_name],
                                  "RandomForest", validation, "complexity_selected_inner_cv",
                                  len(excluded), complexity_rmse)
        row["selected_candidate"] = complexity_name
        results.append(row)
        prediction_tables.append(pred)
        # C: predefined V3-guided sensitivity screen on TRAIN only. Same test IDs.
        quality_train = train[data.iloc[train].sensitivity.to_numpy(float) >= QUALITY_THRESHOLD]
        if len(quality_train) >= MIN_QUALITY_TRAIN_N:
            for model in ("RandomForest", "DummyMean"):
                row, pred = score_holdout(data, quality_train, test, BASE, {}, model, validation,
                                          "quality_train_sens_ge_095", len(excluded),
                                          quality_condition="train_sensitivity_ge_0.95_test_all")
                results.append(row)
                prediction_tables.append(pred)
        else:
            cv_details[f"{validation}/C/skipped"] = {
                "reason": f"quality train n={len(quality_train)} below preset minimum {MIN_QUALITY_TRAIN_N}"}
    results, predictions = validate_week04(data, split_map, results, prediction_tables)
    q = quality_summary(data, split_map)
    results.to_csv(HERE / "experiment_results.csv", index=False)
    predictions.to_csv(HERE / "holdout_predictions.csv", index=False)
    pd.DataFrame(cv_rows).to_csv(HERE / "inner_cv_results.csv", index=False)
    q.to_csv(HERE / "quality_distribution.csv", index=False)
    quality_flags(data).to_csv(HERE / "quality_flag_summary.csv", index=False)
    (HERE / "experiment_design.json").write_text(json.dumps({
        "input_source": str(SOURCE.relative_to(ROOT)), "n_total": len(data),
        "base_features": BASE, "feature_sets": FEATURE_SETS,
        "rf_base": RF_BASE, "complexity_candidates": COMPLEXITY,
        "quality_threshold": QUALITY_THRESHOLD, "minimum_quality_train_n": MIN_QUALITY_TRAIN_N,
        "inner_cv": "3-fold GroupKFold on 500 m easting groups within each outer train; no inner buffer",
        "inner_folds": cv_details}, indent=2), encoding="utf-8")
    plot_results(results)
    print(results[["validation", "experiment", "model", "selected_candidate", "train_n", "test_n", "R2", "RMSE", "MAE"]].to_string(index=False))
    print(f"Output: {HERE}")


if __name__ == "__main__":
    main()
