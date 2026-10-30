"""Compare RF and train-mean Dummy on the exact Week04 baseline splits."""

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "outputs" / "baseline_comparison"
BASELINE = ROOT / "outputs" / "baseline_model"
os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split

from train_baseline import FEATURES, SEED, TARGET, load_data, spatial_split, validate_data


VALIDATIONS = ("random_80_20", "spatial_block_500m_buffer")
MODELS = ("RandomForest", "DummyMean")


def splits(data: pd.DataFrame, coords: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    all_indices = np.arange(len(data))
    random_train, random_test = train_test_split(all_indices, test_size=.2, random_state=SEED)
    spatial_train, spatial_test, spatial_excluded, _ = spatial_split(coords)
    return {
        VALIDATIONS[0]: (random_train, random_test, np.array([], dtype=int)),
        VALIDATIONS[1]: (spatial_train, spatial_test, spatial_excluded),
    }


def score_model(model, data: pd.DataFrame, train: np.ndarray, test: np.ndarray) -> tuple[dict, np.ndarray]:
    x_train = data.iloc[train][list(FEATURES)]
    y_train = data.iloc[train][TARGET].to_numpy(float)
    x_test = data.iloc[test][list(FEATURES)]
    actual = data.iloc[test][TARGET].to_numpy(float)
    model.fit(x_train, y_train)
    predicted = model.predict(x_test)
    return {
        "MAE": float(mean_absolute_error(actual, predicted)),
        "RMSE": float(np.sqrt(mean_squared_error(actual, predicted))),
        "R2": float(r2_score(actual, predicted)),
    }, predicted


def verify_existing(data: pd.DataFrame, split_map: dict,
                    metrics: pd.DataFrame, predictions: pd.DataFrame) -> tuple[bool, list[str]]:
    """Check RF scores, test IDs and predictions against immutable prior outputs."""
    issues = []
    prior_metrics = pd.read_csv(BASELINE / "model_metrics.csv").set_index("validation")
    prior_predictions = pd.read_csv(BASELINE / "validation_predictions.csv")
    prior_diagnostics = json.loads((BASELINE / "data_diagnostics.json").read_text(encoding="utf-8"))
    for validation, (train, test, excluded) in split_map.items():
        old = prior_metrics.loc[validation]
        new = metrics.loc[(metrics.validation == validation) & (metrics.model == "RandomForest")].iloc[0]
        for column, expected in (("train_n", len(train)), ("test_n", len(test)),
                                 ("buffer_excluded_n", len(excluded))):
            if int(old[column]) != expected:
                issues.append(f"{validation}: old {column}={old[column]}, current={expected}")
        for column in ("MAE", "RMSE", "R2"):
            if not np.isclose(float(old[column]), float(new[column]), atol=1e-8, rtol=0):
                issues.append(f"{validation}: old {column}={old[column]}, current={new[column]}")
        old_points = prior_predictions.loc[prior_predictions.validation == validation]
        new_points = predictions.loc[(predictions.validation == validation)
                                     & (predictions.model == "RandomForest")]
        expected_ids = data.iloc[test].point_id.astype(str).to_numpy()
        if not np.array_equal(old_points.point_id.astype(str).to_numpy(), expected_ids):
            issues.append(f"{validation}: prior test point IDs/order differ from reproduced split")
        if not np.array_equal(new_points.point_id.astype(str).to_numpy(), expected_ids):
            issues.append(f"{validation}: comparison test IDs/order differ from reproduced split")
        if not np.allclose(old_points.actual_agbd, new_points.actual_agbd, atol=1e-8, rtol=0):
            issues.append(f"{validation}: prior test AGBD values differ")
        if not np.allclose(old_points.predicted_agbd, new_points.predicted_agbd,
                           atol=1e-8, rtol=0):
            issues.append(f"{validation}: prior RF predictions differ")
        if validation == VALIDATIONS[1]:
            for label, indices in (("train", train), ("test", test),
                                   ("buffer_excluded", excluded)):
                key = f"spatial_{label}_point_ids"
                if prior_diagnostics["splits"][key] != data.iloc[indices].point_id.tolist():
                    issues.append(f"{validation}: prior {label} point IDs/order differ")
    return not issues, issues


def plot_comparison(metrics: pd.DataFrame, predictions: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(15, 9), constrained_layout=True)
    grid = fig.add_gridspec(2, 6, height_ratios=(1.3, 1))
    colors = {"RandomForest": "#2563a6", "DummyMean": "#bd5a35"}
    for validation, columns, title in ((VALIDATIONS[0], slice(0, 3), "Random 80/20"),
                                       (VALIDATIONS[1], slice(3, 6), "Spatial block + 500 m buffer")):
        ax = fig.add_subplot(grid[0, columns])
        section = predictions.loc[predictions.validation == validation]
        extent = max(section.actual_agbd.max(), section.predicted_agbd.max()) * 1.07
        for model in MODELS:
            rows = section.loc[section.model == model]
            ax.scatter(rows.actual_agbd, rows.predicted_agbd, s=31, alpha=.8,
                       label=model, color=colors[model])
        ax.plot([0, extent], [0, extent], "k--", lw=1, label="1:1")
        ax.set(xlim=(0, extent), ylim=(0, extent), title=title,
               xlabel="GEDI estimated AGBD (Mg/ha)", ylabel="Predicted AGBD (Mg/ha)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)

    for i, metric in enumerate(("MAE", "RMSE", "R2")):
        ax = fig.add_subplot(grid[1, slice(i * 2, (i + 1) * 2)])
        x = np.arange(2)
        for offset, model in ((-.18, MODELS[0]), (.18, MODELS[1])):
            values = [metrics.loc[(metrics.validation == validation)
                                  & (metrics.model == model), metric].iloc[0]
                      for validation in VALIDATIONS]
            ax.bar(x + offset, values, width=.36, color=colors[model], label=model)
        ax.set(xticks=x, xticklabels=["Random", "Spatial"], title=metric,
               ylabel="Mg/ha" if metric != "R2" else "R²")
        ax.axhline(0, color="black", lw=.7)
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("Held-out AGBD: Random Forest vs training-mean Dummy", fontsize=14)
    fig.savefig(OUTPUT / "comparison_chart.png", dpi=160)
    plt.close(fig)


def summary_text(metrics: pd.DataFrame, data: pd.DataFrame, split_map: dict,
                 reproduced: bool, issues: list[str]) -> str:
    rows = [
        "# Random Forest와 학습 평균 Dummy 비교",
        "",
        "## 방법과 재현 확인",
        "",
        "- 입력: Week03 SCL 통과 90개 점의 `ndvi`, `evi`, `ndmi`; 목표: GEDI L4A 추정 `agbd`(Mg/ha).",
        "- 기존 `train_baseline.py`의 `SEED=42`, 무작위 80/20 분할과 `spatial_split` 함수를 그대로 사용했다. 각 검증에서 두 모델에 동일한 학습·검증 **인덱스와 순서**를 전달했다.",
        "- 공간 분할의 500m 버퍼 제외점은 두 모델 모두 학습·검증에 사용하지 않았다. Dummy 평균은 전체 90개가 아닌 해당 **학습 집합의 AGBD 평균**으로 학습했다.",
        f"- 실행 환경: scikit-learn {sklearn.__version__}; Random Forest는 `n_estimators=100, random_state=42, n_jobs=1`, Dummy는 `strategy='mean'`.",
        "- 기존 `baseline_model/model_metrics.csv` 및 `validation_predictions.csv`와 RF의 검증점 ID·순서·실제값·예측값·지표를 대조했다.",
        f"- **기존 RF 결과 재현: {'성공' if reproduced else '불일치'}** (지표·예측값 허용오차 1e-8). 기존 결과 파일은 변경하지 않았다.",
    ]
    if issues:
        rows.extend(["", "불일치 항목:", *[f"- {item}" for item in issues]])
    rows.extend(["", "## 검증 데이터 성능", "",
                 "| 분할 | 모델 | 학습 | 검증 | MAE (Mg/ha) | RMSE (Mg/ha) | R² |",
                 "|---|---|---:|---:|---:|---:|---:|"])
    for validation in VALIDATIONS:
        for model in MODELS:
            row = metrics.loc[(metrics.validation == validation) & (metrics.model == model)].iloc[0]
            rows.append(f"| {validation} | {model} | {int(row.train_n)} | {int(row.test_n)} | "
                        f"{row.MAE:.2f} | {row.RMSE:.2f} | {row.R2:.3f} |")
    rows.extend(["", "## 해석", ""])
    rows.append(
        "**결론:** 무작위 분할에서는 RF의 MAE만 근소하게 낮고 RMSE·R²는 Dummy보다 나쁘다. "
        "공간 분리에서는 RF가 세 지표 모두 Dummy보다 나쁘다. "
        "현재 표본과 단일 분할에서는 RF가 단순 학습 평균보다 일관되게 낫다는 근거가 없다."
    )
    rows.append("")
    for validation in VALIDATIONS:
        rf = metrics.loc[(metrics.validation == validation) & (metrics.model == MODELS[0])].iloc[0]
        dummy = metrics.loc[(metrics.validation == validation) & (metrics.model == MODELS[1])].iloc[0]
        train, test, excluded = split_map[validation]
        delta_mae, delta_rmse = rf.MAE - dummy.MAE, rf.RMSE - dummy.RMSE
        rows.append(
            f"- **{validation}:** RF − Dummy MAE {delta_mae:+.2f} Mg/ha, "
            f"RMSE {delta_rmse:+.2f} Mg/ha, R² {rf.R2 - dummy.R2:+.3f}. "
            f"학습 AGBD 평균 {data.iloc[train][TARGET].mean():.2f} Mg/ha; "
            f"검증 AGBD 평균 {data.iloc[test][TARGET].mean():.2f} Mg/ha. "
            f"버퍼 제외 {len(excluded)}개."
        )
    rows.extend([
        "- Dummy의 R²가 반드시 0인 것은 아니다. `strategy='mean'`은 **학습** 평균을 예측하지만 R²의 분모 기준은 **검증** 집합 평균이기 때문이다.",
        "- 이 한 번의 분할 결과는 통계적 유의성 검정이 아니다. RF와 Dummy의 우열을 다른 위치·기간에 일반화할 수 없다.",
        "- 가능한 원인: 90개 중 한 GEDI 궤도의 점만 사용해 공간 대표성이 낮고, NDVI·EVI·NDMI 세 변수만으로 수고·수종·밀도·지형 차이를 담기 어렵다. 공간 홀드아웃에서는 학습 55개와 동쪽 검증 20개의 분포가 다를 수 있다. 반사도 `/65535` 변환은 미확정이며 특히 EVI에 영향을 줄 수 있다. 이 설명들은 진단 가능한 가설이지 인과적 결론은 아니다.",
        "- GEDI L4A AGBD는 기존 모델의 추정값이며 현장 실측값이 아니다. 이 점수는 독립 현장 검증 성능이나 금정산 전체 정확도가 아니다.",
        "- CSV에는 네 실험의 검증 점수와 분할 크기, 최소 학습·검증 거리, 버퍼 제외 수가 있다. 그래프 위쪽은 **동일 축의 검증점 실제값·두 모델 예측값**, 아래쪽은 지표 비교다.",
        "",
    ])
    return "\n".join(rows)


def main() -> None:
    data = load_data()
    coords, _ = validate_data(data)
    split_map = splits(data, coords)
    rows, prediction_rows = [], []
    for validation, (train, test, excluded) in split_map.items():
        if len(np.intersect1d(train, test)) or len(np.intersect1d(excluded, train)) or len(np.intersect1d(excluded, test)):
            raise AssertionError("Split overlap detected")
        distance = float(np.linalg.norm(coords[train, None, :] - coords[None, test, :], axis=2).min())
        for model_name, model in ((MODELS[0], RandomForestRegressor(n_estimators=100, random_state=SEED, n_jobs=1)),
                                  (MODELS[1], DummyRegressor(strategy="mean"))):
            scores, predicted = score_model(model, data, train, test)
            rows.append({"validation": validation, "model": model_name, "train_n": len(train),
                         "test_n": len(test), "buffer_excluded_n": len(excluded),
                         "min_train_test_distance_m": distance, **scores})
            prediction_rows.append(pd.DataFrame({
                "validation": validation, "model": model_name,
                "point_id": data.iloc[test].point_id.to_numpy(),
                "actual_agbd": data.iloc[test][TARGET].to_numpy(float),
                "predicted_agbd": predicted,
            }))
    metrics = pd.DataFrame(rows)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    reproduced, issues = verify_existing(data, split_map, metrics, predictions)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUTPUT / "comparison_metrics.csv", index=False)
    plot_comparison(metrics, predictions)
    (OUTPUT / "comparison_summary.md").write_text(
        summary_text(metrics, data, split_map, reproduced, issues), encoding="utf-8")
    print(metrics.to_string(index=False))
    print(f"Existing RF scores and predictions reproduced: {reproduced}")
    if issues:
        print("Differences: " + "; ".join(issues))
    print(f"Saved outputs to {OUTPUT}")


if __name__ == "__main__":
    main()
