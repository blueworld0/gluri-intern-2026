"""Step 5: GEDI L4A V3 vs V2.1 over Hongcheon 2021.

Questions: does V3 change the values, is V3 more self-consistent at repeat
locations, and does a model trained on V3 labels explain more (R²)?
Needs Earthdata credentials (EARTHDATA_USERNAME / EARTHDATA_PASSWORD or _netrc).
"""
from pathlib import Path
import json
import sys
import urllib.parse

sys.path.insert(0, str(Path(__file__).resolve().parent))
# step2_model fixes the GDAL CA bundle path (Korean venv path) on import.
from step2_model import CRS, NODATA, PROC, RAW, FIG, aligned, grid, read_float, windows  # noqa: E402
from step3_diagnose import POWER, metrics, rf, spatial_oof  # noqa: E402
from step4_chm import CHM_NAMES, CHM_PATH  # noqa: E402
from run_analysis import BBOX, BEAMS, TEMPORAL, credential_strategy, AUTH_HELP  # noqa: E402

from contextlib import ExitStack  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
from rasterio.features import geometry_mask  # noqa: E402
from pyproj import Transformer  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

SHORT_NAME = 'GEDI_L4A_AGB_Density_V3_2508'
RAW_V3 = RAW / 'gedi_l4a_v3_2021'
BEAM_VARS = ('shot_number', 'lat_lowestmode', 'lon_lowestmode', 'agbd', 'agbd_se',
             'l4a_quality_flag_rel3', 'l2a_quality_flag_rel3', 'degrade_flag', 'degrade_include_flag',
             'elev_highestreturn_outlier_flag', 'predictor_limit_flag', 'response_limit_flag',
             'sensitivity', 'solar_elevation', 'delta_time')
LC_VARS = ('leaf_off_flag',)


def h5_ok(path):
    try:
        with h5py.File(path, 'r') as f:
            beams = [b for b in BEAMS if b in f]
            for b in beams:
                ds = f[b]['agbd']
                if ds.size:
                    ds[-1]
            return bool(beams)
    except (OSError, KeyError):
        return False


def download():
    import earthaccess
    strategy = credential_strategy()
    if strategy is None:
        print(AUTH_HELP)
        sys.exit(2)
    if not earthaccess.login(strategy=strategy, persist=False).authenticated:
        raise RuntimeError('Earthdata authentication failed')
    granules = earthaccess.search_data(short_name=SHORT_NAME, bounding_box=BBOX, temporal=TEMPORAL, count=-1)
    print(f'V3 granules: {len(granules)}', flush=True)
    RAW_V3.mkdir(parents=True, exist_ok=True)
    paths = []
    for i, g in enumerate(granules, 1):
        link = next(l for l in g.data_links() if urllib.parse.urlparse(l).path.endswith('.h5'))
        dest = RAW_V3 / Path(urllib.parse.urlparse(link).path).name
        if not h5_ok(dest):
            if dest.exists():
                dest.unlink()  # only this truncated download
            print(f'[{i}/{len(granules)}] downloading {dest.name}', flush=True)
            earthaccess.download([g], local_path=str(RAW_V3), threads=1)
            if not h5_ok(dest):
                raise RuntimeError(f'Invalid download: {dest.name}')
        else:
            print(f'[{i}/{len(granules)}] cached {dest.name}', flush=True)
        paths.append(dest)
    return paths


def extract(paths):
    west, south, east, north = BBOX
    frames, missing = [], set()
    for path in sorted(paths):
        with h5py.File(path, 'r') as f:
            for beam in BEAMS:
                if beam not in f:
                    continue
                g = f[beam]
                lat, lon = g['lat_lowestmode'][:], g['lon_lowestmode'][:]
                idx = np.flatnonzero((lat >= south) & (lat <= north) & (lon >= west) & (lon <= east))
                if not len(idx):
                    continue
                cols = {}
                for v in BEAM_VARS:
                    if v in g:
                        cols[v] = g[v][idx]
                    else:
                        missing.add(v)
                for v in LC_VARS:
                    if 'land_cover_data' in g and v in g['land_cover_data']:
                        cols[v] = g['land_cover_data'][v][idx]
                    else:
                        missing.add(v)
                frame = pd.DataFrame(cols)
                frame['beam'] = beam
                frame['source_file'] = path.name
                frames.append(frame)
    if missing:
        print('V3 variables not found (skipped):', sorted(missing), flush=True)
    df = pd.concat(frames, ignore_index=True)
    df['shot_number'] = df.shot_number.astype('uint64')
    t = pd.Timestamp('2018-01-01', tz='UTC') + pd.to_timedelta(df.delta_time, unit='s')
    return df[(t >= '2021-01-01') & (t < '2022-01-01')].reset_index(drop=True)


def funnel(df):
    steps = [('a. 전체 점', pd.Series(True, index=df.index))]
    valid = np.isfinite(df.agbd) & df.agbd.ge(0)
    steps.append(('b. AGBD 유효값', valid))
    steps.append(('c. l4a_quality_flag_rel3 = 1', steps[-1][1] & df.l4a_quality_flag_rel3.eq(1)))
    steps.append(('d. degrade_flag = 0 (V2.1과 같은 기준, 최종)', steps[-1][1] & df.degrade_flag.eq(0)))
    strict = steps[-1][1].copy()
    for col in ('elev_highestreturn_outlier_flag', 'predictor_limit_flag', 'response_limit_flag', 'leaf_off_flag'):
        if col in df:
            strict &= df[col].fillna(0).eq(0)
    steps.append(('e. + 구름·안개 이상치, 학습범위 밖, 낙엽기 제외 (엄격)', strict))
    rows = [{'stage': s, 'count': int(m.sum()), 'percent_of_all': 100 * m.mean()} for s, m in steps]
    return pd.DataFrame(rows), steps[3][1], strict


def repeat_pairs(df, label):
    out = []
    for r in (25, 50):
        t = cKDTree(df[['easting', 'northing']].to_numpy())
        p = np.array(sorted(t.query_pairs(r)))
        if len(p):
            p = p[df.source_file.to_numpy()[p[:, 0]] != df.source_file.to_numpy()[p[:, 1]]]
        if len(p) < 3:
            out.append({'set': label, 'radius_m': r, 'pairs': len(p)})
            continue
        a, b = df.agbd.to_numpy()[p[:, 0]], df.agbd.to_numpy()[p[:, 1]]
        out.append({'set': label, 'radius_m': r, 'pairs': len(p), 'r': np.corrcoef(a, b)[0, 1],
                    'median_abs_diff': float(np.median(np.abs(a - b)))})
    return out


def sample_features(df, names):
    with rasterio.open(PROC / 'features_10m.tif') as s:
        v = np.array(list(s.sample(zip(df.easting, df.northing))), dtype='float32')
    with rasterio.open(CHM_PATH) as s:
        c = np.array(list(s.sample(zip(df.easting, df.northing))), dtype='float32')
    v[v == NODATA], c[c == NODATA] = np.nan, np.nan
    for i, n in enumerate(names):
        df[n] = v[:, i]
    for i, n in enumerate(CHM_NAMES):
        df[n] = c[:, i]
    return df


def predict_total(model, feats, profile, polygon):
    wc_path = RAW / 'worldcover/ESA_WorldCover_10m_2021_v200_N36E126_Map.tif'
    pix, agb_sum = 0, 0.0
    with ExitStack() as stack:
        fs = stack.enter_context(rasterio.open(PROC / 'features_10m.tif'))
        cs = stack.enter_context(rasterio.open(CHM_PATH))
        wc = aligned(stack, wc_path, profile, True)
        dst = stack.enter_context(rasterio.open(PROC / 'agb_map_v3_10m.tif', 'w', **profile))
        for win in windows(profile):
            arr = np.concatenate([fs.read(window=win), cs.read(window=win)]).astype('float32')
            ok = (arr != NODATA).all(0) & np.isfinite(arr).all(0)
            inside = geometry_mask([polygon], out_shape=ok.shape, invert=True,
                                   transform=rasterio.windows.transform(win, profile['transform']))
            good = ok & (read_float(wc, win) == 10) & inside
            out = np.full(ok.shape, NODATA, 'float32')
            if good.any():
                pred = np.maximum(model.predict(arr[:, good].T), 0).astype('float32')
                out[good] = pred
                pix += len(pred)
                agb_sum += float(pred.sum(dtype='float64')) * 0.01
            dst.write(out, 1, window=win)
    area = pix * 0.01
    tc = agb_sum * 1.25 * 0.47
    return {'model': 'RF + 수관높이, V3 정답', 'mapped_tree_area_ha': area, 'mean_AGB_Mg_ha': agb_sum / area,
            'mean_carbon_tC_ha': tc / area, 'carbon_tC': tc, 'CO2_equivalent_t': tc * 44 / 12,
            'carbon_low_tC': agb_sum * 1.20 * 0.47, 'carbon_high_tC': agb_sum * 1.30 * 0.49}


def main():
    names = json.loads((PROC / 'feature_names.json').read_text(encoding='utf-8'))
    raw = extract(download())
    tx = Transformer.from_crs('EPSG:4326', CRS, always_xy=True)
    raw['easting'], raw['northing'] = tx.transform(raw.lon_lowestmode.to_numpy(), raw.lat_lowestmode.to_numpy())
    fun, final, strict = funnel(raw)
    raw['v3_final'], raw['v3_strict'] = final, strict
    raw.to_csv(PROC / 'gedi_l4a_v3_2021_all.csv', index=False, encoding='utf-8-sig')
    v2_fun = pd.read_csv(PROC / 'filter_summary.csv')
    print(fun.to_string(index=False), flush=True)

    # 1. Same shot, two versions: how much did the values change?
    v2 = pd.read_csv(PROC / 'gedi_l4a_2021_all.csv', dtype={'shot_number': 'uint64'})
    m = v2[['shot_number', 'agbd', 'l4_quality_flag', 'degrade_flag']].merge(
        raw[['shot_number', 'agbd', 'v3_final']], on='shot_number', suffixes=('_v2', '_v3'))
    both = m[m.l4_quality_flag.eq(1) & m.degrade_flag.eq(0) & m.v3_final
             & m.agbd_v2.between(0, 600) & m.agbd_v3.between(0, 600)]
    same_shot = {'matched_shots': len(m), 'both_final': len(both)}
    if len(both) > 2:
        same_shot.update({'r': float(np.corrcoef(both.agbd_v2, both.agbd_v3)[0, 1]),
                          'mean_v2': float(both.agbd_v2.mean()), 'mean_v3': float(both.agbd_v3.mean()),
                          'median_abs_diff': float((both.agbd_v3 - both.agbd_v2).abs().median())})

    # 2. Repeat consistency and 3. model R², per version.
    v2_train = pd.read_csv(PROC / 'training_all_samples.csv', dtype={'shot_number': 'uint64'})
    v2_train = v2_train[v2_train.valid_target & v2_train.agbd.le(600) & v2_train.complete_features]
    v2_train = sample_features(v2_train.reset_index(drop=True), names).dropna(subset=CHM_NAMES)
    sets = {'V2.1 최종': v2_train}
    for label, mask in (('V3 최종', final), ('V3 엄격', strict)):
        d = raw[mask & raw.agbd.le(600)].reset_index(drop=True)
        d = sample_features(d, names).dropna(subset=names + CHM_NAMES).reset_index(drop=True)
        d['block_id'] = (np.floor(d.easting / 2000).astype(int).astype(str) + '_'
                         + np.floor(d.northing / 2000).astype(int).astype(str))
        sets[label] = d
    pairs, scores, slope = [], [], []
    for label, d in sets.items():
        pairs += repeat_pairs(d, label)
        for fl, feats in (('기존 변수', names), ('기존 + 수관높이', names + CHM_NAMES)):
            scores.append({'labels': label, 'features': fl, 'mean_agbd': d.agbd.mean(),
                           **metrics(d.agbd, spatial_oof(d, feats))})
            print(scores[-1], flush=True)
        means = d.groupby(pd.cut(d.slope, [0, 15, 25, 35, 90]), observed=True).agbd.mean()
        slope.append({'labels': label, **{f'slope {k}': v for k, v in means.items()},
                      'power': d[d.beam.isin(POWER)].agbd.mean(), 'coverage': d[~d.beam.isin(POWER)].agbd.mean()})

    profile, polygon = grid()
    v3 = sets['V3 최종']
    model = rf().fit(v3[names + CHM_NAMES].to_numpy('float32'), v3.agbd.to_numpy())
    carbon = predict_total(model, names + CHM_NAMES, profile, polygon)

    result = {'funnel_v3': fun.to_dict('records'), 'funnel_v2': v2_fun.to_dict('records'),
              'same_shot': same_shot, 'repeat_pairs': pairs, 'model': scores, 'slope_beam_means': slope,
              'carbon_v3': carbon}
    (PROC / 'step5_v3_compare.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, default=float),
                                               encoding='utf-8')
    pd.DataFrame(scores).to_csv(PROC / 'step5_v3_model.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(pairs).to_csv(PROC / 'step5_v3_pairs.csv', index=False, encoding='utf-8-sig')
    pd.set_option('display.width', 200)
    print('\n== same shot V2.1 vs V3 ==', same_shot)
    print(pd.DataFrame(pairs).round(2).to_string(index=False))
    print(pd.DataFrame(scores).round(3).to_string(index=False))
    print(pd.DataFrame(slope).round(1).to_string(index=False))
    print('carbon V3:', {k: round(v, 1) if isinstance(v, float) else v for k, v in carbon.items()})


if __name__ == '__main__':
    main()
