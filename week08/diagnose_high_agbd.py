"""Compare observed AGBD and saved predictions for the spring holdout."""

import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplconfig"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = HERE / "outputs"
TRAIN = HERE.parent / "week07/landcover/forest_evaluation/rededge_forest184_input.csv"
HOLDOUT = OUT / "multiyear_spring_O07991_predictions.csv"
MODELS = {"bands10": "RF, 10 bands", "bands10_with_dsm": "RF, 10 bands + elevation/slope"}
EDGES = [0, 50, 100, 150, 200, np.inf]
LABELS = ["0–<50", "50–<100", "100–<150", "150–<200", "≥200"]


def bin_diagnostics(training, test):
    training_bins = pd.cut(training.agbd, EDGES, labels=LABELS, right=False)
    test_bins = pd.cut(test.agbd_Mg_ha, EDGES, labels=LABELS, right=False)
    rows = []
    for label in LABELS:
        group = test.loc[test_bins.eq(label)]
        row = {
            "agbd_bin": label,
            "train_n": int(training_bins.eq(label).sum()),
            "holdout_n": len(group),
            "observed_median": group.agbd_Mg_ha.median(),
        }
        for model in MODELS:
            residual = group.agbd_Mg_ha - group[model]
            row[f"{model}_median"] = group[model].median()
            row[f"{model}_max"] = group[model].max()
            row[f"{model}_mean_underprediction"] = residual.mean()
            row[f"{model}_underprediction_n"] = int(residual.gt(0).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def threshold_diagnostics(training, test):
    rows = []
    for threshold in (100, 150, 200):
        row = {"agbd_threshold": threshold}
        for cohort, actual in (("train", training.agbd), ("holdout", test.agbd_Mg_ha)):
            row[f"{cohort}_n"] = int(actual.ge(threshold).sum())
            row[f"{cohort}_percent"] = 100 * actual.ge(threshold).mean()
        for model in MODELS:
            row[f"{model}_prediction_n"] = int(test[model].ge(threshold).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def plot_diagnostics(training, test, bins):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), layout="constrained")
    fig.suptitle("Spring holdout  |  GEDI aboveground biomass density", fontsize=15)
    x = np.arange(len(LABELS))
    axes[0].bar(x - 0.19, bins.train_n / len(training) * 100, 0.37, color="#365a7c", label="Train, April 2025 (184)")
    axes[0].bar(x + 0.19, bins.holdout_n / len(test) * 100, 0.37, color="#c67537", label="Holdout, May 2020 (339)")
    axes[0].set(title="Observed distributions", ylabel="Share of cohort (%)",
                xlabel="GEDI AGBD (Mg/ha)", xticks=x,
                xticklabels=["<50", "50–<100", "100–<150", "150–<200", "≥200"])
    axes[0].tick_params(axis="x", labelsize=8)
    axes[0].legend(frameon=False, fontsize=8)
    high = test.agbd_Mg_ha.ge(100)
    upper = float(np.ceil(test.agbd_Mg_ha.max() / 25) * 25)
    for ax, (model, title) in zip(axes[1:], MODELS.items()):
        ax.plot([0, upper], [0, upper], ls="--", color="#6b7280", lw=1, label="Perfect agreement")
        for mask, color, label in ((~high, "#365a7c", "<100 Mg/ha"), (high, "#c67537", "≥100 Mg/ha")):
            ax.scatter(test.loc[mask, "agbd_Mg_ha"], test.loc[mask, model],
                       s=13, alpha=0.6, color=color, linewidths=0, label=label)
        ax.plot(bins.observed_median, bins[f"{model}_median"], "o-", color="#18222e",
                lw=1.5, ms=4, label="Bin medians")
        ax.set(title=title, xlabel="Observed GEDI AGBD (Mg/ha)", ylabel="Predicted AGBD (Mg/ha)",
               xlim=(0, upper), ylim=(0, upper), aspect="equal")
        ax.legend(frameon=False, fontsize=7, loc="upper left")
    for ax in axes:
        ax.grid(axis="y", alpha=0.18)
        ax.set_axisbelow(True)
    fig.savefig(OUT / "high_agbd_diagnostic.png", dpi=180, facecolor="white")
    plt.close(fig)


def main():
    train, test = pd.read_csv(TRAIN), pd.read_csv(HOLDOUT)
    if (len(train) != 184 or len(test) != 339 or not train.point_id.is_unique
            or not test.point_id.is_unique or test.point_id.isin(train.point_id).any()):
        raise ValueError("Expected distinct, fixed April training and spring holdout cohorts")
    if not np.isfinite(train.agbd).all() or not np.isfinite(test[["agbd_Mg_ha", *MODELS]]).all().all():
        raise ValueError("Unexpected missing observations or predictions")
    bins = bin_diagnostics(train, test)
    if bins.train_n.sum() != len(train) or bins.holdout_n.sum() != len(test):
        raise ValueError("AGBD bins do not cover the fixed cohorts")
    thresholds = threshold_diagnostics(train, test)
    bins.to_csv(OUT / "high_agbd_prediction_bins.csv", index=False, encoding="utf-8-sig")
    thresholds.to_csv(OUT / "high_agbd_thresholds.csv", index=False, encoding="utf-8-sig")
    plot_diagnostics(train, test, bins)
    print(thresholds.to_string(index=False))
    print("Maximums:", train.agbd.max(), test[["agbd_Mg_ha", *MODELS]].max().to_dict())


if __name__ == "__main__":
    main()
