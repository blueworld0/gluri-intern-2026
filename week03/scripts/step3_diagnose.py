"""Step 3: why is R² low? Diagnose label noise, footprint mismatch, and subsets."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import rasterio
from scipy.spatial import cKDTree
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / 'data/processed'
NODATA = -9999.0
POWER = {'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011'}
SEED = 42


def rf():
    return RandomForestRegressor(n_estimators=200, min_samples_leaf=5, max_features=0.33,
                                 random_state=SEED, n_jobs=4)


def metrics(y, p):
    y, p = np.asarray(y), np.asarray(p)
    ss = ((y - y.mean()) ** 2).sum()
    return {'n': len(y), 'R2': 1 - ((y - p) ** 2).sum() / ss if ss else np.nan,
            'RMSE': float(np.sqrt(((y - p) ** 2).mean())), 'Bias': float((p - y).mean())}


def spatial_oof(df, feats, target='agbd', log=False):
    """Pooled out-of-fold predictions with 2 km block GroupKFold."""
    X, y = df[feats].to_numpy(), df[target].to_numpy()
    yt = np.log1p(y) if log else y
    pred = np.empty(len(df))
    for tr, te in GroupKFold(5).split(X, y, df.block_id):
        m = rf().fit(X[tr], yt[tr])
        pred[te] = m.predict(X[te])
    return np.expm1(pred) if log else pred


def window_means(df, names, sizes):
    """Mean of each feature over size×size pixel windows around each shot.

    Reads one band at a time and smooths it whole; per-point window reads on the
    compressed GeoTIFF were far too slow.
    """
    from scipy.ndimage import uniform_filter
    outs = {s: np.full((len(df), len(names)), np.nan, dtype='float32') for s in sizes}
    with rasterio.open(PROC / 'features_10m.tif') as src:
        rows, cols = rasterio.transform.rowcol(src.transform, df.easting, df.northing)
        rows, cols = np.asarray(rows), np.asarray(cols)
        for b in range(len(names)):
            arr = src.read(b + 1).astype('float32')
            valid = (arr != NODATA).astype('float32')
            arr[arr == NODATA] = 0
            for s in sizes:
                total = uniform_filter(arr, s, mode='constant')
                count = uniform_filter(valid, s, mode='constant')
                with np.errstate(all='ignore'):
                    outs[s][:, b] = np.where(count[rows, cols] > 0.99, total[rows, cols] / count[rows, cols], np.nan)
            print(f'  window band {b + 1}/{len(names)}', flush=True)
    return {s: pd.DataFrame(o, columns=[f'{n}_w{s}' for n in names], index=df.index) for s, o in outs.items()}


def main():
    names = json.loads((PROC / 'feature_names.json').read_text(encoding='utf-8'))
    df = pd.read_csv(PROC / 'training_all_samples.csv', dtype={'shot_number': 'uint64'})
    df = df[df.valid_target & df.agbd.le(600) & df.complete_features].reset_index(drop=True)
    df['power'] = df.beam.isin(POWER)
    df['night'] = df.solar_elevation.lt(0)
    rows = []

    # 1. Label noise ceiling: GEDI shots within 25 m of each other on different tracks.
    tree = cKDTree(df[['easting', 'northing']].to_numpy())
    pairs = np.array(sorted(tree.query_pairs(25)))
    if len(pairs):
        pairs = pairs[df.source_file.to_numpy()[pairs[:, 0]] != df.source_file.to_numpy()[pairs[:, 1]]]
    if len(pairs):
        a, b = df.agbd.to_numpy()[pairs[:, 0]], df.agbd.to_numpy()[pairs[:, 1]]
        r = np.corrcoef(a, b)[0, 1]
        noise = {'pairs': len(pairs), 'r': r, 'r2': r ** 2,
                 'median_abs_diff': float(np.median(np.abs(a - b))),
                 'mean_abs_diff': float(np.mean(np.abs(a - b)))}
    else:
        noise = {'pairs': 0}
    print('GEDI-GEDI 25 m pairs:', noise, flush=True)

    # 2. Baseline and target transform.
    base = spatial_oof(df, names)
    df['pred_base'] = base
    rows.append({'experiment': 'E0 기준 (10m 픽셀, 31개 변수)', **metrics(df.agbd, base)})
    rows.append({'experiment': 'E1 로그 변환 목표값', **metrics(df.agbd, spatial_oof(df, names, log=True))})

    # 3. Footprint mismatch: 25 m footprint vs 10 m pixel, plus ~10 m geolocation error.
    print('Sampling 3x3 and 5x5 window means...', flush=True)
    windows = window_means(df, names, (3, 5))
    w3, w5 = windows[3], windows[5]
    d3 = pd.concat([df, w3], axis=1).dropna(subset=list(w3.columns))
    d5 = pd.concat([df, w5], axis=1).dropna(subset=list(w5.columns))
    rows.append({'experiment': 'E2 30m 창 평균 변수', **metrics(d3.agbd, spatial_oof(d3, list(w3.columns)))})
    rows.append({'experiment': 'E3 50m 창 평균 변수', **metrics(d5.agbd, spatial_oof(d5, list(w5.columns)))})
    both = list(w3.columns) + names
    rows.append({'experiment': 'E4 10m + 30m 변수 함께', **metrics(d3.agbd, spatial_oof(d3, both))})

    # 4. Cleaner labels: retrain on stricter subsets.
    subsets = {
        'E5 power 빔만': df.power,
        'E6 야간만': df.night,
        'E7 sensitivity ≥ 0.98': df.sensitivity.ge(0.98),
        'E8 power + 야간': df.power & df.night,
        'E9 경사 < 25°': df.slope.lt(25),
    }
    for label, mask in subsets.items():
        sub = df[mask].reset_index(drop=True)
        if len(sub) >= 200 and sub.block_id.nunique() >= 5:
            rows.append({'experiment': label, **metrics(sub.agbd, spatial_oof(sub, names))})
        print(label, int(mask.sum()), flush=True)

    result = pd.DataFrame(rows)
    result.to_csv(PROC / 'step3_experiments.csv', index=False, encoding='utf-8-sig')

    # 5. Where does the baseline fail? Error by AGBD level, slope, beam, time of day.
    df['agb_bin'] = pd.cut(df.agbd, [0, 50, 100, 150, 200, 300, 600])
    df['slope_bin'] = pd.cut(df.slope, [0, 15, 25, 35, 90])
    breakdown = []
    for col in ['agb_bin', 'slope_bin', 'power', 'night']:
        for key, g in df.groupby(col, observed=True):
            breakdown.append({'group': col, 'value': str(key), **metrics(g.agbd, g.pred_base),
                              'mean_obs': g.agbd.mean(), 'mean_pred': g.pred_base.mean()})
    breakdown = pd.DataFrame(breakdown)
    breakdown.to_csv(PROC / 'step3_error_breakdown.csv', index=False, encoding='utf-8-sig')
    (PROC / 'step3_noise.json').write_text(json.dumps(noise, ensure_ascii=False, indent=2), encoding='utf-8')

    pd.set_option('display.width', 160)
    print(result.round(3).to_string(index=False))
    print(breakdown.round(1).to_string(index=False))


if __name__ == '__main__':
    main()
