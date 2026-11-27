"""Compare high-AGBD augmentation with an equally sized representative addition."""

import argparse
import json

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

from compare_topography import APRIL, HOLDOUT, OUT, sample_dsm
from evaluate_february_holdout import scores

ADDITIONS = OUT / "multiyear_spring_O07991_predictions.csv"
EVALUATION_INPUTS = {
    "february": HOLDOUT,
    "may6": OUT / "multiyear_spring_O07930_predictions.csv",
    "may28": OUT / "multiyear_spring_O08261_predictions.csv",
}
COORDS = ["point_id", "easting_m", "northing_m"]
TOPOGRAPHY = ["dsm_elevation_m", "dsm_slope_deg"]
TARGET = "agbd_Mg_ha"
THRESHOLD = 150


def cohorts(bands):
    columns = COORDS + [TARGET, *bands]
    april = pd.read_csv(APRIL / "rededge_forest184_input.csv").rename(columns={"agbd": TARGET})
    may = pd.read_csv(ADDITIONS)
    high = may.loc[may[TARGET].ge(THRESHOLD)]
    if (len(april), len(may), len(high)) != (184, 339, 58):
        raise ValueError("The fixed April and May 10 cohorts have changed")
    # Match the May 10 pool's proportion of observations above the cutoff.
    ordinary = may.loc[may[TARGET].lt(THRESHOLD)]
    n_high = round(len(high) * len(high) / len(may))
    representative = pd.concat([
        high.sample(n=n_high, random_state=42),
        ordinary.sample(n=len(high) - n_high, random_state=42),
    ]).sort_index()
    return {
        "original_rf": april[columns],
        "representative_addition_rf": pd.concat([april[columns], representative[columns]], ignore_index=True),
        "high_agbd_addition_rf": pd.concat([april[columns], high[columns]], ignore_index=True),
    }


def with_topography(df, topo):
    return df.drop(columns=TOPOGRAPHY, errors="ignore").merge(topo, on="point_id", validate="one_to_one")


def evaluate(test_path, key, training, config):
    used = pd.concat([frame[COORDS] for frame in training.values()], ignore_index=True).drop_duplicates()
    test = pd.read_csv(test_path)
    train_images = set(pd.read_csv(ADDITIONS, usecols=["product"])["product"])
    test_images = set(test["product"])
    shared_images = sorted(train_images & test_images)
    if key == "may28" and shared_images:
        raise ValueError("The May 28 holdout must use imagery distinct from the May 10 additions")
    if not used.point_id.is_unique or not test.point_id.is_unique or test.point_id.isin(used.point_id).any():
        raise ValueError("Training and evaluation IDs must be distinct")
    topo = sample_dsm(pd.concat([used, test[COORDS]], ignore_index=True))
    features = config["feature_sets"]["bands10"] + TOPOGRAPHY
    test = with_topography(test, topo)
    test["april_mean"] = training["original_rf"][TARGET].mean()
    names = ["april_mean", *training]
    for name, frame in training.items():
        inputs = with_topography(frame, topo)
        model = RandomForestRegressor(**config["rf_params"])
        model.fit(inputs[features], inputs[TARGET])
        test[name] = model.predict(test[features])
    if key == "february":
        reference = pd.read_csv(OUT / "topography_predictions.csv").set_index("point_id")
        np.testing.assert_allclose(
            test.original_rf, reference.loc[test.point_id, "bands10_with_dsm"], atol=1e-10, rtol=0,
        )
    else:
        np.testing.assert_allclose(test.original_rf, test.bands10_with_dsm, atol=1e-10, rtol=0)
    test["nearest_any_training_m"] = np.hypot(
        test.easting_m.to_numpy()[:, None] - used.easting_m.to_numpy(),
        test.northing_m.to_numpy()[:, None] - used.northing_m.to_numpy(),
    ).min(axis=1)
    separated = test.loc[test.nearest_any_training_m.ge(500)]
    groups = {
        "all": test,
        "ge100": test.loc[test[TARGET].ge(100)],
        "ge150": test.loc[test[TARGET].ge(150)],
        "lt100": test.loc[test[TARGET].lt(100)],
        "separated_500m": separated,
        "separated_500m_ge100": separated.loc[separated[TARGET].ge(100)],
    }
    result = {
        "evaluation_cohort": key, "test_n": len(test),
        "sentinel_products": sorted(test_images),
        "sentinel_products_shared_with_additions": shared_images,
        "n_at_least_500m_from_all_training": len(separated),
        "train": {
            name: {"n": len(frame), "ge150_n": int(frame[TARGET].ge(150).sum())}
            for name, frame in training.items()
        },
        "metrics": {
            scope: {name: scores(frame[TARGET], frame[name]) for name in names}
            for scope, frame in groups.items() if not frame.empty
        },
    }
    test[["point_id", TARGET, "nearest_any_training_m", *names]].to_csv(
        OUT / f"high_agbd_training_{key}_predictions.csv", index=False, encoding="utf-8-sig",
    )
    payload = json.dumps(result, indent=2, allow_nan=False)
    (OUT / f"high_agbd_training_{key}_summary.json").write_text(payload + "\n", encoding="utf-8")
    for scope, metrics in result["metrics"].items():
        print(key, scope, {name: round(val["rmse"], 2) for name, val in metrics.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=EVALUATION_INPUTS, default="february")
    args = parser.parse_args()
    config = json.loads((APRIL / "rededge_config.json").read_text(encoding="utf-8"))
    training = cohorts(config["feature_sets"]["bands10"])
    if args.cohort == "may28":
        del training["high_agbd_addition_rf"]
    path = EVALUATION_INPUTS[args.cohort]
    if not path.is_file():
        parser.error(f"Run spring quality screening first. Missing {path.name}")
    evaluate(path, args.cohort, training, config)


if __name__ == "__main__":
    main()
