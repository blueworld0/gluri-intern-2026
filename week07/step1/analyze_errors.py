"""Diagnose frozen Week06 OOF predictions without fitting or changing inputs."""
from pathlib import Path
import argparse
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(frame):
    def fmt(v):
        return f'{v:.2f}' if isinstance(v, (float, np.floating)) else str(v)
    return '| ' + ' | '.join(frame.columns) + ' |\n| ' + ' | '.join(['---'] * len(frame.columns)) + ' |\n' + '\n'.join(
        '| ' + ' | '.join(fmt(v) for v in row) + ' |' for row in frame.itertuples(index=False, name=None))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).parent / 'outputs')
    out = parser.parse_args().output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use an empty output directory to preserve previous results.')
    pred_path = ROOT / 'week06/step8/outputs/predictions.csv'
    input_path = ROOT / 'week06/step7/outputs/new441_evaluation_input.csv'
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in (pred_path, input_path)}
    raw = pd.read_csv(pred_path)
    points = pd.read_csv(input_path, dtype={'shot_number': str})
    df = raw.loc[(raw['case'] == 'main_all441') & (raw['model'] == 'week05_spatial_B_fixed')].copy()
    assert len(df) == 441 and df.point_id.is_unique and points.point_id.is_unique
    metadata = ['point_id', 'agbd', 'latitude', 'longitude', 'spatial_group', 'easting_m', 'northing_m',
                'ndvi', 'ndmi', 'agbd_se', 'evi_review_required', 'forest_status', 'beam']
    joined = df.merge(points[metadata], on='point_id', validate='one_to_one', suffixes=('', '_input'), indicator=True)
    assert joined['_merge'].eq('both').all()
    for col in ['agbd', 'latitude', 'longitude']:
        assert np.allclose(joined[col], joined[col + '_input'], rtol=0, atol=1e-10)
    assert joined.spatial_group.eq(joined.spatial_group_input).all()
    assert np.allclose(df.residual, df.agbd - df.prediction)
    df = joined.drop(columns=['_merge', 'agbd_input', 'latitude_input', 'longitude_input', 'spatial_group_input'])
    assert np.isfinite(df[['agbd', 'prediction', 'residual', 'easting_m', 'northing_m']]).all().all()
    df['absolute_error'] = df.residual.abs()
    df['squared_error'] = df.residual ** 2
    df['high_agbd'] = df.agbd >= 100
    df = df.sort_values(['absolute_error', 'point_id'], ascending=[False, True]).reset_index(drop=True)
    df.insert(0, 'error_rank', np.arange(1, len(df) + 1))
    df['top20'] = df.error_rank <= 20
    top = df.head(20)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / 'all441_error_diagnostics.csv', index=False)
    top.to_csv(out / 'top20_absolute_errors.csv', index=False)
    grids = df.groupby(['spatial_group', 'fold']).agg(
        n=('point_id', 'size'), high_agbd_n=('high_agbd', 'sum'), top20_n=('top20', 'sum'),
        mean_agbd=('agbd', 'mean'), mean_prediction=('prediction', 'mean'),
        mean_residual=('residual', 'mean'), mae=('absolute_error', 'mean'), mse=('squared_error', 'mean')).reset_index()
    grids['rmse'] = np.sqrt(grids.pop('mse'))
    grids['top20_rate'] = grids.top20_n / grids.n
    grids = grids.sort_values(['top20_n', 'rmse'], ascending=False)
    grids.to_csv(out / 'spatial_group_metrics.csv', index=False)
    bins = df.groupby('agbd_bin', observed=True).agg(n=('point_id', 'size'),
        mean_agbd=('agbd', 'mean'), mean_prediction=('prediction', 'mean'),
        mean_residual=('residual', 'mean'), mae=('absolute_error', 'mean'), mse=('squared_error', 'mean')).reset_index()
    bins['rmse'] = np.sqrt(bins.pop('mse'))
    bins.to_csv(out / 'agbd_bin_metrics.csv', index=False)

    fig, ax = plt.subplots(figsize=(8, 7), layout='constrained')
    scatter = ax.scatter(df.agbd, df.prediction, c=df.fold, cmap='viridis', vmin=1, vmax=5, s=25, alpha=.65)
    ax.scatter(top.agbd, top.prediction, facecolors='none', edgecolors='red', s=100, label='Top 20 absolute errors')
    for row in top.itertuples():
        ax.annotate(str(row.error_rank), (row.agbd, row.prediction), xytext=(4, 4), textcoords='offset points', fontsize=7)
    limit = max(df.agbd.max(), df.prediction.max()) * 1.06
    ax.plot([0, limit], [0, limit], '--', color='gray', label='1:1')
    ax.set(xlim=(0, limit), ylim=(0, limit), xlabel='GEDI estimated AGBD (Mg/ha)',
           ylabel='OOF predicted AGBD (Mg/ha)', title='Week06 Spatial RF: observed vs predicted (n=441)')
    ax.set_aspect('equal')
    ax.legend(loc='upper left')
    fig.colorbar(scatter, ax=ax, ticks=range(1, 6), label='Spatial fold')
    fig.savefig(out / 'observed_vs_predicted.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 7), layout='constrained', sharex=True, sharey=True)
    bound = max(abs(df.residual.min()), abs(df.residual.max()))
    for ax in axes:
        for coord in ['easting_m', 'northing_m']:
            lo, hi = df[coord].min() / 1000, df[coord].max() / 1000
            for edge in np.arange(np.floor(lo / 5) * 5, hi + 5, 5):
                (ax.axvline if coord == 'easting_m' else ax.axhline)(edge, color='.9', linewidth=.6, zorder=0)
        ax.set_aspect('equal')
        ax.set(xlabel='UTM 52N easting (km)', ylabel='Northing (km)')
        ax.ticklabel_format(useOffset=False, style='plain')
    sc = axes[0].scatter(df.easting_m / 1000, df.northing_m / 1000, c=df.residual,
                         cmap='RdBu_r', vmin=-bound, vmax=bound, s=18)
    axes[0].scatter(top.easting_m / 1000, top.northing_m / 1000, facecolors='none', edgecolors='black', s=75)
    for row in top.itertuples():
        axes[0].annotate(str(row.error_rank), (row.easting_m / 1000, row.northing_m / 1000),
                         xytext=(4, 3), textcoords='offset points', fontsize=6)
    axes[0].set_title('Residuals and top 20 error ranks')
    fig.colorbar(sc, ax=axes[0], shrink=.7, label='GEDI minus prediction (Mg/ha)')
    axes[1].scatter(df.easting_m / 1000, df.northing_m / 1000, c='.8', s=14, label='All 441 points')
    high = df.loc[df.high_agbd]
    axes[1].scatter(high.easting_m / 1000, high.northing_m / 1000, c='darkorange', s=22, label='GEDI AGBD >= 100')
    axes[1].scatter(top.easting_m / 1000, top.northing_m / 1000, facecolors='none', edgecolors='black', s=75, label='Top 20 errors')
    axes[1].set_title('High AGBD and top 20 locations')
    axes[1].legend(fontsize=8)
    fig.suptitle('Week06 OOF error locations; 5 km grid (no basemap)')
    fig.savefig(out / 'error_locations.png', dpi=180)
    plt.close(fig)

    summary = dict(n=len(df), rmse=float(np.sqrt(df.squared_error.mean())), mae=float(df.absolute_error.mean()),
                   top20_high_agbd_n=int(top.high_agbd.sum()), top20_underprediction_n=int((top.residual > 0).sum()),
                   top20_squared_error_share=float(top.squared_error.sum() / df.squared_error.sum()),
                   top20_grid_n=int(top.spatial_group.nunique()), all_grid_n=int(df.spatial_group.nunique()),
                   high_agbd_n=int(df.high_agbd.sum()), high_agbd_underprediction_n=int((high.residual > 0).sum()),
                   high_agbd_grid_n=int(high.spatial_group.nunique()))
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    # Korean report keeps all interpretation tied to descriptive diagnostics.
    report = f'''# Week07 Step 1 — week06 예측 오차 진단

## 분석 조건

Week06 `main_all441` / `week05_spatial_B_fixed`의 저장된 공간 OOF 예측 441점만 분석했다. 각 점은 학습에 포함되지 않은 fold에서 얻은 예측이다. 기존 5km 격자·500m 버퍼 설정을 그대로 사용하며 재학습·튜닝·데이터 수정은 하지 않았다. 잔차 = GEDI AGBD − 예측값이며, 양수는 과소 예측이다. 절대오차 내림차순으로 상위 20점을 선택했고 동률은 point_id로 정렬했다.

## 주요 결과

- 전체 RMSE {summary['rmse']:.2f} Mg/ha, MAE {summary['mae']:.2f} Mg/ha로 week06 기본 결과와 일치한다.
- 절대오차 상위 20점 중 AGBD 100 Mg/ha 이상은 {summary['top20_high_agbd_n']}점, 과소 예측은 {summary['top20_underprediction_n']}점이다.
- 상위 20점은 전체 표본의 {20/441:.1%}이며 전체 제곱오차의 {summary['top20_squared_error_share']:.1%}를 차지한다. 이는 MSE 기여율이며 RMSE의 가산 기여율이 아니다.
- 상위 20점은 전체 {summary['all_grid_n']}개 격자 중 {summary['top20_grid_n']}개에 분포한다. 아래 표의 표본 수·고 AGBD 수·상위 20점 비율을 함께 보아야 한다.
- AGBD 100 이상 전체 {summary['high_agbd_n']}점 중 {summary['high_agbd_underprediction_n']}점이 과소 예측되며, 이 고 AGBD 표본은 {summary['high_agbd_grid_n']}개 격자에 있다. 100 기준은 기술적 진단 구간이며 모델 선택 기준이 아니다.
- 상위 20점이 가장 많은 격자는 `{grids.iloc[0].spatial_group}`로 {int(grids.iloc[0].top20_n)}점이다(전체 {int(grids.iloc[0].n)}점, 고 AGBD {int(grids.iloc[0].high_agbd_n)}점). 단일 격자에만 국한된 현상은 아니며, 지도에서도 여러 위치에서 고 AGBD 표본과 큰 과소 예측이 겹친다.

## 관측값과 예측값

![관측-예측 산점도](observed_vs_predicted.png)

빨간 원은 절대오차 상위 20점이며 숫자는 오차 순위다. 점선 아래는 과소 예측이다.

## 위치와 공간 격자별 오차

![오차 위치도](error_locations.png)

원은 상위 20점, 왼쪽 숫자는 오차 순위이며 오른쪽 주황색 점은 고 AGBD 표본이다. 배경지도 없이 EPSG:32652 좌표와 기존 5km 격자를 표시했다. 행정구역·산림 경계 판정은 하지 않았다.

{table(grids.drop(columns=['top20_rate']))}

격자별 `top20_rate`는 CSV에 0~1 비율로 저장했다. 표본이 적은 격자의 RMSE와 비율은 불안정할 수 있다. 고 AGBD 표본의 지역 분포와 오차 분포가 겹치는지만 기술하며, 지역 자체가 오차의 원인이라는 결론이나 공간 군집의 통계적 유의성을 주장하지 않는다.

## AGBD 구간별 오차

{table(bins)}

## 상위 20점

{table(top[['error_rank', 'spatial_group', 'fold', 'latitude', 'longitude', 'agbd', 'prediction', 'residual', 'absolute_error']])}

위 표의 좌표는 표시 자리수를 줄였다. `top20_absolute_errors.csv`에는 정확한 좌표, point_id, NDVI/EVI/NDMI, AGBD 표준오차, EVI 점검 플래그를 함께 보존했다.

## 해석과 다음 확인

고 AGBD에서의 과소 예측과 격자별 편차를 먼저 진단한다. 산점도와 격자 집계만으로 반사율 오류·혼합 피복·지수 포화·GEDI 추정 불확실성 중 어느 것이 원인인지 확정할 수 없다. 다음 단계에서는 상위 오차점의 원본 반사율·주변 픽셀 및 피복을 검토한다. 오차가 크다는 이유로 표본을 삭제하거나, 이미 확인한 OOF 점수에 맞춰 최종 검증 성능을 주장하지 않는다.

대상은 산림 미검증 SCL 식생 후보이며 GEDI L4A 모델 추정치를 비교 기준으로 사용한다. 같은 날짜·궤도 자료의 공간 CV로 독립 현장 검증이나 시간·궤도 일반화를 입증하지 않는다.

## 재현

저장소 루트에서 기존 Python 환경으로 실행한다. 출력 폴더가 비어 있지 않으면 중단하므로 재실행에는 새 폴더를 지정한다.

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\step1\\analyze_errors.py --output-dir .\\week07\\step1\\reproduction_outputs
```
'''
    (out / 'report.md').write_text(report, encoding='utf-8')
    assert all(sha(ROOT / p) == digest for p, digest in hashes.items())
    validation = dict(status='passed', n=441, unique_point_ids=True, metadata_join='one_to_one, complete',
                      coordinates_and_target_match=True, residual_recalculated=True,
                      input_sha256_unchanged=hashes, model_refit=False)
    (out / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
