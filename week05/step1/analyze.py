"""Week05 Step1 audit of the unchanged Week04 validation results."""

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WEEK04 = ROOT / "week04"
OUT = Path(__file__).resolve().parent / "outputs"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
sys.path.insert(0, str(WEEK04))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from compare_baselines import splits
from train_baseline import FEATURES, TARGET, load_data, validate_data


VALIDATIONS = ("random_80_20", "spatial_block_500m_buffer")
LABELS = ("Random 80/20", "Spatial block + 500 m buffer")


def check_metrics(actual, predicted, saved):
    calculated = (r2_score(actual, predicted),
                  np.sqrt(mean_squared_error(actual, predicted)),
                  mean_absolute_error(actual, predicted))
    expected = (saved.R2, saved.RMSE, saved.MAE)
    if not np.allclose(calculated, expected, rtol=0, atol=1e-8):
        raise ValueError(f"Week04 metric mismatch: {calculated} versus {expected}")


def make_predictions(data, coords, metrics):
    prior = pd.read_csv(WEEK04 / "outputs/baseline_model/validation_predictions.csv")
    split_map = splits(data, coords)
    result = []
    checks = {}
    for validation in VALIDATIONS:
        train, test, excluded = split_map[validation]
        if any(np.intersect1d(a, b).size for a, b in ((train, test), (train, excluded), (test, excluded))):
            raise ValueError("Overlapping spatial partitions")
        actual = data.iloc[test][TARGET].to_numpy(float)
        ids = data.iloc[test].point_id.astype(str).to_numpy()
        rf = prior.loc[prior.validation.eq(validation)].copy()
        if not np.array_equal(rf.point_id.astype(str).to_numpy(), ids):
            raise ValueError(f"Week04 RF IDs/order differ: {validation}")
        if not np.allclose(rf.actual_agbd.to_numpy(), actual, rtol=0, atol=1e-8):
            raise ValueError(f"Week04 RF targets differ: {validation}")
        dummy_value = float(data.iloc[train][TARGET].mean())
        predictions = {"RandomForest": rf.predicted_agbd.to_numpy(float),
                       "DummyMean": np.full(len(test), dummy_value)}
        for model, predicted in predictions.items():
            saved = metrics.loc[(metrics.validation == validation) & (metrics.model == model)].iloc[0]
            if (saved.train_n, saved.test_n, saved.buffer_excluded_n) != (len(train), len(test), len(excluded)):
                raise ValueError(f"Week04 split counts differ: {validation}, {model}")
            check_metrics(actual, predicted, saved)
            result.append(pd.DataFrame({"validation": validation, "model": model,
                                        "point_id": ids, "actual_agbd": actual,
                                        "predicted_agbd": predicted,
                                        "residual_actual_minus_predicted": actual - predicted}))
        distance = np.linalg.norm(coords[train, None, :] - coords[None, test, :], axis=2)
        if not np.isclose(distance.min(), metrics.loc[metrics.validation.eq(validation), "min_train_test_distance_m"].iloc[0]):
            raise ValueError("Minimum train/test distance differs from Week04")
        checks[validation] = {"train_n": len(train), "test_n": len(test),
                              "buffer_excluded_n": len(excluded),
                              "minimum_train_test_distance_m": float(distance.min()),
                              "train_test_id_overlap_n": len(set(data.iloc[train].point_id) & set(ids)),
                              "train_test_raster_pixel_overlap_n": len(set(zip(data.iloc[train].raster_row, data.iloc[train].raster_col)) & set(zip(data.iloc[test].raster_row, data.iloc[test].raster_col))),
                              "train_mean_agbd": dummy_value,
                              "test_mean_agbd": float(actual.mean())}
    return pd.concat(result, ignore_index=True), checks


def quality(data, existing, checks):
    stats = {}
    for column in (TARGET, *FEATURES):
        s = pd.to_numeric(data[column], errors="coerce")
        q1, q3 = s.quantile([.25, .75])
        lo, hi = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
        stats[column] = {"count": int(s.count()), "min": float(s.min()),
                         "q1": float(q1), "median": float(s.median()),
                         "mean": float(s.mean()), "q3": float(q3),
                         "max": float(s.max()), "sample_std": float(s.std()),
                         "iqr_lower_fence": float(lo), "iqr_upper_fence": float(hi),
                         "iqr_flagged_n": int(((s < lo) | (s > hi)).sum())}
    band_stats = pd.read_csv(WEEK04 / "outputs/reflectance_validation/band_statistics.csv")
    bands_2025 = band_stats.loc[band_stats.year.eq(2025)]
    sensitivity = pd.read_csv(WEEK04 / "outputs/reflectance_validation/index_comparison_statistics.csv")
    report = {"source": str((WEEK04.parent / "week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv").relative_to(ROOT)),
              "total_rows": len(data), "total_columns": len(data.columns),
              "missing_by_column": {k: int(v) for k, v in data.isna().sum().items()},
              "numeric_summary": stats,
              "agbd_ge_200_n": int((data.agbd >= 200).sum()),
              "agbd_ge_250_n": int((data.agbd >= 250).sum()),
              "invalid_index_range_n": {"ndvi": int((data.ndvi.abs() > 1).sum()),
                                          "ndmi": int((data.ndmi.abs() > 1).sum()),
                                          "evi_abs_gt_2": int((data.evi.abs() > 2).sum())},
              "duplicate_point_ids_n": int(data.point_id.duplicated().sum()),
              "duplicate_raster_pixels_n": int(data.duplicated(["raster_row", "raster_col"]).sum()),
              "date_gap_days": {"min": int(data.date_gap_days.min()), "max": int(data.date_gap_days.max())},
              "same_gedi_file_n": int(data.gedi_file.nunique()),
              "raw_65535_any_n": int(data.raw_65535_any.sum()),
              "raw_zero_any_n": int(data.raw_zero_any.sum()),
              "scl_boundary_disagreement_n": int(data.scl_boundary_disagreement.sum()),
              "split_checks": checks,
              "reflectance_2025_tiff_raw_band_summary": {
                  row.band: {"median_raw_uint16": float(row.median),
                             "median_if_divided_by_65535": float(row.median / 65535),
                             "valid_65535_pixels": int(row.valid_65535)}
                  for row in bands_2025.itertuples()},
              "evi_mean_by_scale_assumption": {
                  row.scenario: float(row.mean) for row in sensitivity.loc[sensitivity["index"].eq("EVI")].itertuples()},
              "week04_reflectance_assessment": existing}
    return report


def plot_residuals(predictions, data):
    rf = predictions.loc[predictions.model.eq("RandomForest")]
    specs = (("actual_vs_predicted.png", "actual"),
             ("residual_vs_predicted.png", "residual"),
             ("residual_histogram.png", "hist"))
    for filename, kind in specs:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
        for ax, validation, label in zip(axes, VALIDATIONS, LABELS):
            rows = rf.loc[rf.validation.eq(validation)]
            actual, pred = rows.actual_agbd, rows.predicted_agbd
            residual = rows.residual_actual_minus_predicted
            if kind == "actual":
                ax.scatter(actual, pred, alpha=.8)
                ax.plot([0, 310], [0, 310], "k--", linewidth=1)
                ax.set(xlabel="GEDI estimated AGBD (Mg/ha)", ylabel="RF predicted AGBD (Mg/ha)", xlim=(0, 310), ylim=(0, 310))
            elif kind == "residual":
                ax.scatter(pred, residual, alpha=.8)
                ax.axhline(0, color="black", linestyle="--", linewidth=1)
                ax.set(xlabel="RF predicted AGBD (Mg/ha)", ylabel="Actual − predicted (Mg/ha)")
            else:
                ax.hist(residual, bins=np.arange(-150, 201, 25), edgecolor="white")
                ax.axvline(0, color="black", linestyle="--", linewidth=1)
                ax.set(xlabel="Actual − predicted (Mg/ha)", ylabel="Count")
            ax.set_title(f"{label} (n={len(rows)})")
            ax.grid(alpha=.2)
        fig.savefig(OUT / filename, dpi=160)
        plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, values, label in zip(axes, (data.agbd, rf.loc[rf.validation.eq(VALIDATIONS[0]), "actual_agbd"], rf.loc[rf.validation.eq(VALIDATIONS[1]), "actual_agbd"]), ("All SCL-valid points", *LABELS)):
        ax.hist(values, bins=np.arange(0, 351, 25), edgecolor="white")
        ax.set(xlabel="GEDI estimated AGBD (Mg/ha)", ylabel="Count", title=f"{label} (n={len(values)})", xlim=(0, 350))
        ax.grid(axis="y", alpha=.2)
    fig.savefig(OUT / "agbd_distribution.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, feature in zip(axes, FEATURES):
        ax.hist(data[feature], bins=15, edgecolor="white")
        ax.set(xlabel=feature.upper(), ylabel="Count", title=feature.upper())
        ax.grid(axis="y", alpha=.2)
    fig.savefig(OUT / "feature_distributions.png", dpi=160)
    plt.close(fig)


def residual_summary(predictions):
    rows = []
    for (validation, model), group in predictions.groupby(["validation", "model"], sort=False):
        for label, subset in (("all", group), ("low_actual_le_50", group.loc[group.actual_agbd <= 50]),
                              ("high_actual_ge_150", group.loc[group.actual_agbd >= 150]),
                              ("very_high_actual_ge_200", group.loc[group.actual_agbd >= 200])):
            rows.append({"validation": validation, "model": model, "agbd_group": label,
                         "n": len(subset), "mean_actual": subset.actual_agbd.mean(),
                         "mean_predicted": subset.predicted_agbd.mean(),
                         "mean_residual_actual_minus_predicted": subset.residual_actual_minus_predicted.mean(),
                         "underpredicted_n": int((subset.residual_actual_minus_predicted > 0).sum()),
                         "overpredicted_n": int((subset.residual_actual_minus_predicted < 0).sum())})
    return pd.DataFrame(rows)


def report(metrics, summary, q):
    lines = ["# Week05 Step1: Week04 모델 성능 및 예측 오차 분석", "",
             "## 자료와 방법", "",
             "- Week03 SCL 통과 90점의 GEDI L4A 추정 AGBD(Mg/ha)를 목표값으로, Sentinel-2 NDVI·EVI·NDMI를 입력으로 사용했다. GEDI AGBD는 현장 실측값이 아니다.",
             "- Week04의 RF 지표·검증점 예측 CSV를 재사용했다. Dummy는 Week04의 분할 함수와 학습 집합 평균으로 예측값을 재현했다. 네 실험의 지표를 기존 CSV와 1e-8 이내로 대조했다.",
             "- 잔차 = 실제 GEDI 추정 AGBD − 모델 예측값. 양수는 과소 예측이다. 원본 데이터·Week04 결과는 변경하지 않았다.", "",
             "## Week04 실험 결과", "",
             "| 분할 | 모델 | 학습 n | 평가 n | R² | RMSE | MAE |", "|---|---|---:|---:|---:|---:|---:|"]
    for row in metrics.itertuples():
        lines.append(f"| {row.validation} | {row.model} | {row.train_n} | {row.test_n} | {row.R2:.3f} | {row.RMSE:.2f} | {row.MAE:.2f} |")
    lines += ["", "RMSE·MAE 단위는 Mg/ha이다. Random 분할에서 RF의 MAE는 Dummy보다 0.47 낮지만 RMSE는 1.39 높다. Spatial 분할에서 RF는 Dummy보다 MAE 9.64, RMSE 13.45 높다. 두 분할 모두 RF R²가 음수이며 단일 분할 결과다.", "",
              "## Random Split과 Spatial Split", "",
              f"- Random: 학습 72점, 평가 18점, 최소 학습·평가 거리 {q['split_checks'][VALIDATIONS[0]]['minimum_train_test_distance_m']:.1f}m.",
              f"- Spatial: 동쪽 끝 1km easting 블록 20점을 평가하고, 500m 버퍼의 15점을 제외해 55점으로 학습했다. 최소 거리 {q['split_checks'][VALIDATIONS[1]]['minimum_train_test_distance_m']:.1f}m. 두 분할의 평가점·학습 규모가 다르므로 점수 차이 전부를 공간 효과로 단정할 수 없다.",
              "- 공간 분할은 좌표만으로 정해지며 학습·평가 point ID 및 raster pixel 중복은 0개다. 검증점 AGBD를 분할 선택에 사용하지 않았다. 동일 GEDI 궤도와 공간 자기상관, 분석 과정의 반복적 모델 선택에 따른 평가 편향 가능성은 남아 있다. 이 코드는 500m 밖의 자기상관 소멸을 증명하지 않는다.", "",
              "## 예측 오차", ""]
    for validation in VALIDATIONS:
        for group_name in ("low_actual_le_50", "high_actual_ge_150", "very_high_actual_ge_200"):
            r = summary.loc[(summary.validation == validation) & (summary.model == "RandomForest") & (summary.agbd_group == group_name)].iloc[0]
            if r.n:
                lines.append(f"- {validation}, {group_name}: n={r.n}, 평균 잔차 {r.mean_residual_actual_minus_predicted:+.2f} Mg/ha, 과소 예측 {r.underpredicted_n}점, 과대 예측 {r.overpredicted_n}점.")
    lines += ["- 높은/낮은 AGBD 구간은 사전에 정한 기술적 기준(150 이상 / 50 이하 Mg/ha)이며, 표본이 작아 경향의 일반화는 어렵다. 잔차와 예측값의 관계는 별도 산점도에서 확인할 수 있다.", "",
              "## 데이터 품질과 결합 점검", "",
              f"- 최종 자료 {q['total_rows']}행·{q['total_columns']}열. 결측이 있는 컬럼: " + ", ".join(f"{k} {v}개" for k, v in q['missing_by_column'].items() if v) + ". 학습 입력·목표·좌표 결측은 0개다.",
              f"- AGBD: 최소 {q['numeric_summary']['agbd']['min']:.2f}, 최대 {q['numeric_summary']['agbd']['max']:.2f}, 평균 {q['numeric_summary']['agbd']['mean']:.2f}, 중앙값 {q['numeric_summary']['agbd']['median']:.2f}, 표본 표준편차 {q['numeric_summary']['agbd']['sample_std']:.2f} Mg/ha. IQR 1.5배 상한 {q['numeric_summary']['agbd']['iqr_upper_fence']:.2f} 초과 1점; 200 이상 {q['agbd_ge_200_n']}점, 250 이상 {q['agbd_ge_250_n']}점. 제거하지 않았다."]
    for feature in FEATURES:
        x = q["numeric_summary"][feature]
        lines.append(f"- {feature.upper()}: 범위 {x['min']:.3f}–{x['max']:.3f}, 중앙값 {x['median']:.3f}, IQR 이상치 표시 {x['iqr_flagged_n']}점. 분포는 `feature_distributions.png` 참조.")
    lines += [f"- NDVI·NDMI의 [-1, 1] 밖 값 0개, |EVI|>2 값 0개. 최종점의 raw 65535 표시 {q['raw_65535_any_n']}개, raw 0 표시 {q['raw_zero_any_n']}개, SCL 경계 판정 불일치 {q['scl_boundary_disagreement_n']}개. 단순 범위 검사만으로 반사율의 물리적 단위는 확인되지 않는다.",
              f"- 최종점의 GEDI–Sentinel 시차는 {q['date_gap_days']['min']}–{q['date_gap_days']['max']}일, GEDI 원본 파일 {q['same_gedi_file_n']}개, 중복 point ID {q['duplicate_point_ids_n']}개, 중복 raster pixel {q['duplicate_raster_pixels_n']}개다. 원시 결합 과정은 좌표를 UTM 52N으로 변환해 지수 raster 픽셀을 추출했고, SCL 필터로 최종 90점을 선택했다.",
              "- **반사율 단위 미확정:** Week03은 Copernicus Browser Analytical Raw uint16 TIFF를 `/65535`로 변환했다. Week04 TIFF 메타데이터에는 실질적인 radiometric scale/offset 근거가 없고, `/10000`은 민감도 비교 시나리오일 뿐이다. NDVI·NDMI는 동일 배율이면 값이 유지되지만 EVI는 상수항 때문에 달라진다. Week04 90점 비교에서 `/65535` 재계산값과 저장 지수는 거의 일치했으나 이는 물리적 단위의 검증이 아니다.",
              f"- 2025 TIFF 원시 중앙값은 B02 {q['reflectance_2025_tiff_raw_band_summary']['B02']['median_raw_uint16']:.0f}, B04 {q['reflectance_2025_tiff_raw_band_summary']['B04']['median_raw_uint16']:.0f}, B08 {q['reflectance_2025_tiff_raw_band_summary']['B08']['median_raw_uint16']:.0f}이다. `/65535` 가정 시 각각 {q['reflectance_2025_tiff_raw_band_summary']['B02']['median_if_divided_by_65535']:.3f}, {q['reflectance_2025_tiff_raw_band_summary']['B04']['median_if_divided_by_65535']:.3f}, {q['reflectance_2025_tiff_raw_band_summary']['B08']['median_if_divided_by_65535']:.3f}이지만 물리 반사율임이 입증되지는 않았다. 2025 TIFF 유효 픽셀의 65535도 존재한다. 90점 EVI 평균은 `/65535` 가정 {q['evi_mean_by_scale_assumption']['browser_65535_assumption']:.3f}, `/10000` 민감도 가정 {q['evi_mean_by_scale_assumption']['dn_10000_sensitivity']:.3f}으로 차이가 크다.", "",
              "## 성능 저하 원인: 확인 사실과 가설", "",
              "**확인 사실:** 표본은 단일 GEDI 파일의 90점이며 세 식생지수만 입력했다. Spatial 학습은 55점으로 줄고 평가 위치는 동쪽 끝에 집중된다. RF는 Spatial에서 Dummy보다 성능이 낮다. 반사율 변환은 미확정이다.", "",
              "**가설:** 제한된 입력이 수종·수관 구조·지형 차이를 설명하지 못하거나, 적은 표본과 위치별 분포 차이로 RF가 과적합했을 수 있다. EVI의 반사율 배율 가정도 영향을 줄 수 있다. 각 요인의 기여도는 현재 결과만으로 분리할 수 없다.", "",
              "## Week05 우선 작업", "",
              "1. 원본 TIFF 생성 설정/evalscript 또는 검증된 FLOAT32 반사율로 배율·offset을 확인한다. 확인 전에는 기존 지수를 수정하지 않는다.",
              "2. 버퍼를 둔 여러 공간 블록과 반복 분할로 성능 변동성을 평가하고, 각 fold에서 Dummy를 동일 학습점으로 비교한다.",
              "3. 평가점의 AGBD·입력 지수·공간 분포 차이와 극단값의 영향도를 점검한다. 극단값을 임의 제거하지 않는다.",
              "4. 검증 가능한 시공간 일치 자료와 수관 구조·지형 등 추가 변수를 검토하되, 독립 현장 자료가 없으면 현장 정확도로 해석하지 않는다.", "",
              "## 산출물", "",
              "`model_comparison.csv`, `validation_residuals.csv`, `residual_summary.csv`, `data_quality.json`과 PNG 그래프 5개. `data_quality.json`은 전체 컬럼 결측과 분할 검사 수치를 포함한다.", ""]
    return "\n".join(lines)


def main():
    data = load_data()
    coords, _ = validate_data(data)
    metrics = pd.read_csv(WEEK04 / "outputs/baseline_comparison/comparison_metrics.csv")
    original_rf = pd.read_csv(WEEK04 / "outputs/baseline_model/model_metrics.csv")
    rf_metrics = metrics.loc[metrics.model.eq("RandomForest")]
    for row in original_rf.itertuples():
        current = rf_metrics.loc[rf_metrics.validation.eq(row.validation)].iloc[0]
        if not np.allclose([row.R2, row.RMSE, row.MAE],
                           [current.R2, current.RMSE, current.MAE], rtol=0, atol=1e-8):
            raise ValueError(f"Two Week04 RF metric sources differ: {row.validation}")
    predictions, checks = make_predictions(data, coords, metrics)
    existing = json.loads((WEEK04 / "outputs/reflectance_validation/tiff_metadata.json").read_text(encoding="utf-8"))["interpretation"]
    q = quality(data, existing, checks)
    summary = residual_summary(predictions)
    OUT.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUT / "model_comparison.csv", index=False)
    predictions.to_csv(OUT / "validation_residuals.csv", index=False)
    summary.to_csv(OUT / "residual_summary.csv", index=False)
    (OUT / "data_quality.json").write_text(json.dumps(q, ensure_ascii=False, indent=2), encoding="utf-8")
    plot_residuals(predictions, data)
    (OUT / "report.md").write_text(report(metrics, summary, q), encoding="utf-8")
    print(metrics.to_string(index=False))
    print(summary.loc[summary.model.eq("RandomForest")].to_string(index=False))
    print(f"Saved to {OUT}")


if __name__ == "__main__":
    main()
