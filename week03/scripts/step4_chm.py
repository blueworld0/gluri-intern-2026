"""Step 4: add Meta 1 m canopy height (CHM) as features, compare with/without, redo map + carbon."""
from pathlib import Path
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
# step2_model fixes the GDAL CA bundle path (Korean venv path) on import.
from step2_model import (CRS, NODATA, PROC, RAW, FIG, aligned, grid, read_float,  # noqa: E402
                         windows)
from step3_diagnose import POWER, metrics, rf, spatial_oof  # noqa: E402

from contextlib import ExitStack  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402
from rasterio.features import geometry_mask  # noqa: E402
from rasterio.warp import reproject, transform_bounds  # noqa: E402
from rasterio.windows import from_bounds, Window  # noqa: E402
from scipy.ndimage import maximum_filter, uniform_filter  # noqa: E402

CHM_URL = '/vsicurl/https://dataforgood-fb-data.s3.amazonaws.com/forests/v1/alsgedi_global_v6_float/chm/{}.tif'
TILES = ('132110312', '132110330')  # zoom-9 quadkeys covering the Hongcheon box
CHM_PATH = PROC / 'chm_10m.tif'
CHM_NAMES = ['chm_mean', 'chm_max', 'chm_mean_w3', 'chm_max_w3', 'chm_cover2m']
ROWS = 512


def build_chm(profile):
    """Resample the 1 m CHM to the 10 m grid as mean, max and canopy cover (>2 m)."""
    if CHM_PATH.exists():
        print('CHM cached:', CHM_PATH.name, flush=True)
        return
    shape = (profile['height'], profile['width'])
    total, count = np.zeros(shape, 'float64'), np.zeros(shape, 'float64')
    cover = np.zeros(shape, 'float64')
    peak = np.full(shape, np.nan, 'float32')
    start = time.monotonic()
    env = dict(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', CPL_VSIL_CURL_CHUNK_SIZE=10485760,
               GDAL_HTTP_MAX_RETRY=5, GDAL_HTTP_RETRY_DELAY=2)
    with rasterio.Env(**env):
        for tile in TILES:
            with rasterio.open(CHM_URL.format(tile)) as src:
                b = transform_bounds(CRS, src.crs, *rasterio.transform.array_bounds(*shape, profile['transform']))
                win = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
                win = win.intersection(Window(0, 0, src.width, src.height))
                print(f'tile {tile}: rows {win.height}, cols {win.width}', flush=True)
                for off in range(0, int(win.height), ROWS):
                    w = Window(win.col_off, win.row_off + off, win.width, min(ROWS, win.height - off))
                    chunk = src.read(1, window=w).astype('float32')
                    t = src.window_transform(w)
                    for data, acc, how in ((chunk, total, Resampling.sum),
                                           (np.ones_like(chunk), count, Resampling.sum),
                                           ((chunk > 2).astype('float32'), cover, Resampling.sum)):
                        tmp = np.zeros(shape, 'float32')
                        reproject(data, tmp, src_transform=t, src_crs=src.crs, dst_transform=profile['transform'],
                                  dst_crs=CRS, resampling=how, src_nodata=None, dst_nodata=0)
                        acc += tmp
                    tmp = np.full(shape, np.nan, 'float32')
                    reproject(chunk, tmp, src_transform=t, src_crs=src.crs, dst_transform=profile['transform'],
                              dst_crs=CRS, resampling=Resampling.max, src_nodata=None, dst_nodata=np.nan)
                    np.fmax(peak, tmp, out=peak)
                    done = off + w.height
                    print(f'  {tile}: {done}/{win.height} rows, {(time.monotonic() - start) / 60:.1f} min', flush=True)
    with np.errstate(all='ignore'):
        mean = np.where(count > 0.5, total / count, np.nan).astype('float32')
        cov = np.where(count > 0.5, cover / count, np.nan).astype('float32')
    bands = [mean, peak]
    for arr, f in ((mean, 'mean'), (peak, 'max')):
        filled = np.nan_to_num(arr)
        valid = np.isfinite(arr).astype('float32')
        if f == 'mean':
            with np.errstate(all='ignore'):
                bands.append((uniform_filter(filled, 3) / uniform_filter(valid, 3)).astype('float32'))
        else:
            bands.append(maximum_filter(filled, 3).astype('float32'))
    bands.append(cov)
    prof = dict(profile, count=len(bands))
    with rasterio.open(CHM_PATH, 'w', **prof) as dst:
        for i, (arr, name) in enumerate(zip(bands, CHM_NAMES), 1):
            dst.write(np.where(np.isfinite(arr), arr, NODATA).astype('float32'), i)
            dst.set_band_description(i, name)
    print('CHM written', flush=True)


def load_training(names):
    df = pd.read_csv(PROC / 'training_all_samples.csv', dtype={'shot_number': 'uint64'})
    df = df[df.valid_target & df.agbd.le(600) & df.complete_features].reset_index(drop=True)
    with rasterio.open(CHM_PATH) as src:
        vals = np.array(list(src.sample(zip(df.easting, df.northing))), dtype='float32')
    vals[vals == NODATA] = np.nan
    for i, n in enumerate(CHM_NAMES):
        df[n] = vals[:, i]
    df['power'] = df.beam.isin(POWER)
    return df.dropna(subset=CHM_NAMES).reset_index(drop=True)


def experiments(df, names):
    terrain = ['elevation', 'slope', 'aspect']
    sets = {
        'A. 기존 (Sentinel-2 + 지형)': names,
        'B. 수관높이 + 지형만': CHM_NAMES + terrain,
        'C. 기존 + 수관높이': names + CHM_NAMES,
    }
    rows, oof = [], {}
    for label, feats in sets.items():
        p = spatial_oof(df, feats)
        oof[label] = p
        rows.append({'experiment': label, 'n_features': len(feats), **metrics(df.agbd, p)})
        print(rows[-1], flush=True)
    clean = df[df.power & df.slope.lt(25)].reset_index(drop=True)
    for label, feats in (('D. 기존, power빔+경사<25°', names), ('E. 기존+수관높이, power빔+경사<25°', names + CHM_NAMES)):
        rows.append({'experiment': label, 'n_features': len(feats), **metrics(clean.agbd, spatial_oof(clean, feats))})
        print(rows[-1], flush=True)
    # 250 m cell averages (less GEDI label noise), with and without CHM.
    df['cell'] = (np.floor(df.easting / 250).astype(int).astype(str) + '_'
                  + np.floor(df.northing / 250).astype(int).astype(str))
    allf = names + CHM_NAMES
    g = df.groupby('cell').agg(n=('agbd', 'size'), agbd=('agbd', 'mean'), easting=('easting', 'mean'),
                               northing=('northing', 'mean'), **{f: (f, 'mean') for f in allf}).reset_index()
    g = g[g.n >= 3].reset_index(drop=True)
    g['block_id'] = (np.floor(g.easting / 2000).astype(int).astype(str) + '_'
                     + np.floor(g.northing / 2000).astype(int).astype(str))
    for label, feats in (('F. 기존, 250m 칸 평균', names), ('G. 기존+수관높이, 250m 칸 평균', allf)):
        rows.append({'experiment': label, 'n_features': len(feats), **metrics(g.agbd, spatial_oof(g, feats))})
        print(rows[-1], flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(PROC / 'step4_chm_experiments.csv', index=False, encoding='utf-8-sig')
    # Error by AGBD level, with vs without CHM (does saturation ease?).
    df['agb_bin'] = pd.cut(df.agbd, [0, 50, 100, 150, 200, 300, 600])
    by_bin = []
    for key, idx in df.groupby('agb_bin', observed=True).groups.items():
        for label in ('A. 기존 (Sentinel-2 + 지형)', 'C. 기존 + 수관높이'):
            y, p = df.agbd[idx], oof[label][idx]
            by_bin.append({'agb_bin': str(key), 'model': label[:2].strip(), 'n': len(idx),
                           'mean_obs': y.mean(), 'mean_pred': p.mean(), 'bias': (p - y).mean()})
    pd.DataFrame(by_bin).to_csv(PROC / 'step4_bias_by_bin.csv', index=False, encoding='utf-8-sig')
    return res, oof


def fit_and_map(df, names, profile, polygon):
    """Fit RF on all points with CHM features, predict the 10 m map, total the carbon."""
    feats = names + CHM_NAMES
    model = rf().fit(df[feats].to_numpy('float32'), df.agbd.to_numpy())
    imp = pd.DataFrame({'feature': feats, 'importance': model.feature_importances_}).sort_values('importance', ascending=False)
    imp.to_csv(PROC / 'step4_feature_importance.csv', index=False, encoding='utf-8-sig')
    wc_path = RAW / 'worldcover/ESA_WorldCover_10m_2021_v200_N36E126_Map.tif'
    pix, agb_sum, tree_pix, miss = 0, 0.0, 0, 0
    with ExitStack() as stack:
        fs = stack.enter_context(rasterio.open(PROC / 'features_10m.tif'))
        cs = stack.enter_context(rasterio.open(CHM_PATH))
        wc = aligned(stack, wc_path, profile, True)
        dst = stack.enter_context(rasterio.open(PROC / 'agb_map_chm_10m.tif', 'w', **profile))
        for win in windows(profile):
            arr = np.concatenate([fs.read(window=win), cs.read(window=win)]).astype('float32')
            ok = (arr != NODATA).all(0) & np.isfinite(arr).all(0)
            inside = geometry_mask([polygon], out_shape=ok.shape, invert=True,
                                   transform=rasterio.windows.transform(win, profile['transform']))
            tree = (read_float(wc, win) == 10) & inside
            good = ok & tree
            tree_pix += int(tree.sum())
            miss += int((tree & ~ok).sum())
            out = np.full(ok.shape, NODATA, 'float32')
            if good.any():
                pred = np.maximum(model.predict(arr[:, good].T), 0).astype('float32')
                out[good] = pred
                pix += len(pred)
                agb_sum += float(pred.sum(dtype='float64')) * 0.01
            dst.write(out, 1, window=win)
    area = pix * 0.01
    tc = agb_sum * 1.25 * 0.47
    return {'model': 'RF + 수관높이', 'tree_cover_area_ha': tree_pix * 0.01, 'mapped_tree_area_ha': area,
            'coverage_percent': 100 * pix / tree_pix, 'mean_AGB_Mg_ha': agb_sum / area,
            'mean_carbon_tC_ha': tc / area, 'carbon_tC': tc, 'CO2_equivalent_t': tc * 44 / 12,
            'carbon_low_tC': agb_sum * 1.20 * 0.47, 'carbon_high_tC': agb_sum * 1.30 * 0.49}, imp


def figures(df, oof, res, profile):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'Malgun Gothic'
    plt.rcParams['axes.unicode_minus'] = False
    # 1. Observed vs predicted, before / after CHM.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, label in zip(axes, ('A. 기존 (Sentinel-2 + 지형)', 'C. 기존 + 수관높이')):
        p = oof[label]
        s = res.set_index('experiment').loc[label]
        ax.scatter(df.agbd, p, s=6, alpha=.3, edgecolors='none')
        ax.plot([0, 600], [0, 600], 'k--', lw=1)
        ax.set(xlim=(0, 600), ylim=(0, 600), xlabel='GEDI AGBD (Mg/ha)', ylabel='공간 교차검증 예측 (Mg/ha)',
               title=f'{label}\nR²={s.R2:.3f}, RMSE={s.RMSE:.1f}')
    fig.tight_layout()
    fig.savefig(FIG / 'chm_obs_vs_pred_compare.png', dpi=170)
    plt.close(fig)
    # 2. R² bar chart for all experiments.
    fig, ax = plt.subplots(figsize=(9, 4.5))
    colors = ['#9aa5b1' if '수관' not in e else '#2e7d32' for e in res.experiment]
    ax.barh(res.experiment[::-1], res.R2[::-1], color=colors[::-1])
    for i, v in enumerate(res.R2[::-1]):
        ax.text(v + .005, i, f'{v:.2f}', va='center')
    ax.set(xlabel='R² (2km 공간 블록 교차검증)', title='수관높이 추가 전(회색)·후(초록) 성능')
    fig.tight_layout()
    fig.savefig(FIG / 'chm_r2_compare.png', dpi=170)
    plt.close(fig)
    # 3. Maps side by side.
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, (path, title) in zip(axes, ((PROC / 'agb_map_RF_10m.tif', '기존 (Sentinel-2 + 지형)'),
                                        (PROC / 'agb_map_chm_10m.tif', '기존 + 수관높이'))):
        with rasterio.open(path) as src:
            arr = src.read(1, out_shape=(src.height // 2, src.width // 2), masked=True)
            b = src.bounds
        im = ax.imshow(np.ma.masked_less(arr, 0), extent=(b.left, b.right, b.bottom, b.top), cmap='YlGn', vmin=0, vmax=400)
        ax.set(title=title, xticks=[], yticks=[])
    fig.colorbar(im, ax=axes, label='AGBD (Mg/ha)', shrink=.8)
    fig.savefig(FIG / 'chm_map_compare.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    names = json.loads((PROC / 'feature_names.json').read_text(encoding='utf-8'))
    profile, polygon = grid()
    build_chm(profile)
    df = load_training(names)
    print('training points with CHM:', len(df), flush=True)
    res, oof = experiments(df, names)
    carbon, imp = fit_and_map(df, names, profile, polygon)
    old = pd.read_csv(PROC / 'carbon_totals.csv')
    old = old[old.model == 'RF'].iloc[0].to_dict()
    old['model'] = 'RF 기존'
    cmp = pd.DataFrame([{k: old.get(k) for k in carbon}, carbon])
    cmp.to_csv(PROC / 'step4_carbon_compare.csv', index=False, encoding='utf-8-sig')
    figures(df, oof, res, profile)
    pd.set_option('display.width', 200)
    print(res.round(3).to_string(index=False))
    print(cmp.round(1).to_string(index=False))
    print(imp.head(10).round(3).to_string(index=False))


if __name__ == '__main__':
    main()
