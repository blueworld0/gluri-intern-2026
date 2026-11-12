"""Reusable inference/evaluation functions. Requires Colab scientific dependencies.
Call configure(config, output_directory, projected_crs) before inference.
DBH inputs: H, CD in metres. AGB input DBH in centimetres; output kilograms.
"""
from pathlib import Path
import inspect
import math
import json, importlib.metadata
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.mask import mask
from shapely.geometry import Point, box, Polygon
from shapely.ops import unary_union
from scipy.optimize import linear_sum_assignment
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, maximum_filter, label as ndi_label
from skimage.segmentation import watershed
from rasterio.features import shapes as raster_shapes
from shapely.geometry import shape
AUDIT=[]
def configure(config, output_directory, projected_crs):
    global CFG, OUT, WORK_CRS, default
    CFG=config
    OUT=Path(output_directory)
    OUT.mkdir(parents=True,exist_ok=True)
    WORK_CRS=projected_crs
    default=config['agb_default']
def audit(stage,item,reason):
    AUDIT.append(dict(stage=stage,item=str(item),reason=str(reason)))


def require_columns(df, cols, label):
    missing = set(cols) - set(df.columns)
    if missing:
        raise ValueError(f'{label} 필드 누락: {missing}; 실제 필드={list(df.columns)}')

def offset_xy(e, n, distance, azimuth_deg):
    rad = np.deg2rad(azimuth_deg)
    return (e + distance * np.sin(rad), n + distance * np.cos(rad))

def estimate_dbh(H_m, CD_m, coefficients=None):
    c = coefficients or CFG['jucker']
    H = np.asarray(H_m, dtype=float)
    CD = np.asarray(CD_m, dtype=float)
    return np.where((H > 0) & (CD > 0), c['a'] * (H * CD) ** c['b'] * np.exp(c.get('g', 0) ** 2 / 2), np.nan)

def estimate_agb(dbh_cm, beta0=None, beta1=None):
    b0 = default['beta0'] if beta0 is None else beta0
    b1 = default['beta1'] if beta1 is None else beta1
    D = np.asarray(dbh_cm, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(D > 0, np.exp(b0 + b1 * np.log(D)), np.nan)

def kg_to_t_ha(kg, area_m2):
    return np.asarray(kg) / 1000 * 10000 / np.asarray(area_m2)

def carbon_kg(agb_kg, fraction=0.5):
    return np.asarray(agb_kg) * fraction

def tile_coordinates(centers, buffer_m=50):
    tiles = set()
    for x, y in centers:
        for e in range(int((x - buffer_m) // 1000), int((x + buffer_m) // 1000) + 1):
            for n in range(int((y - buffer_m) // 1000), int((y + buffer_m) // 1000) + 1):
                tiles.add((e * 1000, n * 1000))
    return sorted(tiles)

def product_site_entry(product, site):
    records = product.get('data', product).get('siteCodes', [])
    return next((r for r in records if r.get('siteCode') == site), dict(siteCode=site, availableMonths=[], availableReleases=[]))

def select_aop_pair(rgb_entry, chm_entry, target_year, allow_adjacent=False, priority=(2022, 2020), include_provisional=True):
    common = sorted(set(rgb_entry.get('availableMonths', [])) & set(chm_entry.get('availableMonths', [])))
    years = [int(target_year)] + ([int(y) for y in priority] if allow_adjacent else [])
    selected = next((y for y in years if any((m.startswith(str(y) + '-') for m in common))), None)
    if selected is None:
        raise ValueError(f'RGB/CHM 공통 연월 없음: 요청 연도 {years}; 공통 월 {common}')
    months = [m for m in common if m.startswith(str(selected) + '-')]
    provisional = False
    for entry in [rgb_entry, chm_entry]:
        stable = {m for r in entry.get('availableReleases', []) if r.get('release') != 'PROVISIONAL' for m in r.get('availableMonths', [])}
        prov = {m for r in entry.get('availableReleases', []) if r.get('release') == 'PROVISIONAL' for m in r.get('availableMonths', [])}
        provisional |= any((m in prov and m not in stable for m in months))
    if provisional and (not include_provisional):
        raise ValueError('선택 연도에 PROVISIONAL 자료가 있습니다. aop_include_provisional=True 필요; 다른 해로 숨겨 대체하지 않습니다.')
    return dict(imagery_year=selected, common_months=months, include_provisional=bool(provisional), provisional=bool(provisional), adjacent_imagery=selected != int(target_year))

def aop_tiff_files(folder):
    return sorted((p for p in Path(folder).rglob('*') if p.is_file() and p.suffix.lower() in ['.tif', '.tiff']))

def raster_index(folder, expected_resolution):
    rows = []
    for path in sorted(Path(folder).rglob('*')):
        if path.suffix.lower() not in ['.tif', '.tiff']:
            continue
        with rasterio.open(path) as ds:
            if not ds.crs or not ds.crs.is_projected:
                raise ValueError(f'미터 투영 CRS 필요: {path}')
            if not np.allclose(ds.res, expected_resolution, atol=expected_resolution * 0.05):
                audit('resolution', path, f'예상 {expected_resolution}; 실제 {ds.res}')
                continue
            rows.append(dict(path=str(path), crs=ds.crs, geometry=box(*ds.bounds)))
    return rows

def crop_mosaic(index, region, region_crs, dest, bands):
    if dest.exists():
        return dest
    datasets = []
    coverage = []
    try:
        for r in index:
            geom = gpd.GeoSeries([region], crs=region_crs).to_crs(r['crs']).iloc[0]
            if r['geometry'].intersects(geom):
                datasets.append(rasterio.open(r['path']))
                coverage.append(r['geometry'])
        if not datasets:
            raise ValueError('해당 위치를 덮는 타일 없음')
        crs = datasets[0].crs
        if any((d.crs != crs for d in datasets)):
            raise ValueError('mosaic 입력 CRS 불일치')
        geom = gpd.GeoSeries([region], crs=region_crs).to_crs(crs).iloc[0]
        if not unary_union(coverage).buffer(0.2).covers(geom):
            raise ValueError('조사구 주변 영상 coverage 불완전')
        data, trans = merge(datasets, bounds=geom.bounds, indexes=bands)
        profile = datasets[0].profile.copy()
        profile.update(driver='GTiff', height=data.shape[1], width=data.shape[2], count=len(bands), transform=trans, compress='deflate', tiled=True, blockxsize=256, blockysize=256)
        with rasterio.open(dest, 'w', **profile) as out:
            out.write(data)
    finally:
        for d in datasets:
            d.close()
    return dest

def raster_covers_geometry(index, geometry, geometry_crs):
    if not index:
        return False
    bounds = [gpd.GeoSeries([r['geometry']], crs=r['crs']).to_crs(geometry_crs).iloc[0] for r in index]
    return unary_union(bounds).buffer(0.2).covers(geometry)

def zonal_max(path, geom, geom_crs):
    with rasterio.open(path) as ds:
        g = gpd.GeoSeries([geom], crs=geom_crs).to_crs(ds.crs).iloc[0]
        try:
            arr, _ = mask(ds, [g], crop=True, filled=False, indexes=1)
        except ValueError:
            return np.nan
        v = arr.compressed()
        v = v[np.isfinite(v) & (v >= 0)]
        return float(v.max()) if len(v) else np.nan

def predict_crowns(m, rgb, chm, plot):
    cache = OUT / f'{plot.plotID}_raw_predictions.csv'
    if cache.exists():
        raw = pd.read_csv(cache)
    else:
        params = inspect.signature(m.predict_tile).parameters
        path_key = 'path' if 'path' in params else 'raster_path'
        args = {path_key: str(rgb), 'patch_size': CFG['patch_size'], 'patch_overlap': CFG['patch_overlap']}
        if 'iou_threshold' in params:
            args['iou_threshold'] = CFG['nms_iou']
        raw = m.predict_tile(**args)
        if isinstance(raw, tuple):
            raw = raw[0]
        if raw is None or len(raw) == 0:
            raw = pd.DataFrame(columns=['xmin', 'ymin', 'xmax', 'ymax', 'score'])
        raw.to_csv(cache, index=False)
    require_columns(raw, ['xmin', 'ymin', 'xmax', 'ymax', 'score'], 'DeepForest predictions')
    rows = []
    with rasterio.open(rgb) as ds:
        trans = ds.transform
        crs = ds.crs
        for idx, r in raw[raw.score >= CFG['score_threshold']].iterrows():
            corners = [trans * (x, y) for x, y in [(r.xmin, r.ymin), (r.xmax, r.ymin), (r.xmax, r.ymax), (r.xmin, r.ymax)]]
            geom = Polygon(corners)
            geom = gpd.GeoSeries([geom], crs=crs).to_crs(WORK_CRS).iloc[0]
            if not plot.geometry.covers(geom.centroid):
                continue
            h = zonal_max(chm, geom, WORK_CRS)
            if not np.isfinite(h):
                audit('CHM_nodata', f'{plot.plotID}_{idx}', 'CHM 없음')
                continue
            if h < CFG['min_height_m']:
                continue
            width = (r.xmax - r.xmin) * abs(trans.a)
            length = (r.ymax - r.ymin) * abs(trans.e)
            rows.append(dict(det_id=f'{plot.plotID}_{idx}', plotID=plot.plotID, siteID=plot.siteID, field_year=plot.field_year, score=r.score, xmin_px=r.xmin, ymin_px=r.ymin, xmax_px=r.xmax, ymax_px=r.ymax, H_m=h, CD_m=(width + length) / 2, width_m=width, length_m=length, geometry=geom))
    return rows

def label_visibility(f, rasters):
    f = f.copy()
    local = []
    for _, r in f.iterrows():
        local.append(zonal_max(rasters[r.plotID][1], r.geometry.buffer(CFG['visible_radius_m']), WORK_CRS))
    f['local_CHM_m'] = local
    canopy = f.get('canopyPosition', pd.Series('', index=f.index)).fillna('').str.strip().str.lower()
    f['visible_by_position'] = canopy.isin(['full sun', 'open grown', 'partially shaded'])
    f['visible_by_height'] = f.field_H_m.ge(CFG['visible_height_ratio'] * f.local_CHM_m) & f.local_CHM_m.gt(0)
    f['visible'] = f.visible_by_position | f.visible_by_height
    return f

def match_one_to_one(det, truth):
    if det.empty or truth.empty:
        return pd.DataFrame(columns=['det_id', 'individualID', 'plotID', 'distance_m'])
    d = det.reset_index(drop=True)
    t = truth.reset_index(drop=True)
    n, m = (len(d), len(t))
    weights = np.zeros((n, m + n))
    valid = np.zeros((n, m), dtype=bool)
    sizes = pd.to_numeric(t[CFG['match_size']], errors='coerce').fillna(0).to_numpy()
    rank = pd.Series(sizes).rank(method='dense').to_numpy() / (m + 1)
    for i, r in d.iterrows():
        for j, s in t.iterrows():
            if r.geometry.covers(s.geometry):
                valid[i, j] = True
                distance = r.geometry.centroid.distance(s.geometry)
                weights[i, j] = (n + 1) * 10 + rank[j] + 1e-06 / (1 + distance)
    ri, ci = linear_sum_assignment(-weights)
    rows = []
    for i, j in zip(ri, ci):
        if j < m and valid[i, j]:
            rows.append(dict(det_id=d.iloc[i].det_id, individualID=t.iloc[j].individualID, plotID=d.iloc[i].plotID, distance_m=d.iloc[i].geometry.centroid.distance(t.iloc[j].geometry)))
    return pd.DataFrame(rows, columns=['det_id', 'individualID', 'plotID', 'distance_m'])

def grouped_dbh_regression(df, n_splits=5):
    result = pd.Series(np.nan, index=df.index)
    fold_id = pd.Series(np.nan, index=df.index)
    good = df[['H_m', 'CD_m', 'field_DBH_cm']].gt(0).all(axis=1)
    good &= np.isfinite(df[['H_m', 'CD_m', 'field_DBH_cm']]).all(axis=1)
    z = df.loc[good]
    groups = z.plotID
    if groups.nunique() < 2:
        audit('regression', 'GroupKFold', '독립 조사구 2개 미만')
        return (result, fold_id, None)
    X = np.log(z[['H_m', 'CD_m']].to_numpy())
    y = np.log(z.field_DBH_cm.to_numpy())
    folds = GroupKFold(n_splits=min(n_splits, groups.nunique()))
    for fold, (tr, te) in enumerate(folds.split(X, y, groups)):
        if len(tr) < 4:
            audit('regression', fold, '학습 나무 4개 미만')
            continue
        m = LinearRegression().fit(X[tr], y[tr])
        smear = float(np.mean(np.exp(y[tr] - m.predict(X[tr]))))
        result.loc[z.index[te]] = np.exp(m.predict(X[te])) * smear
        fold_id.loc[z.index[te]] = fold
    m = LinearRegression().fit(X, y)
    fit = dict(intercept=float(m.intercept_), log_H=float(m.coef_[0]), log_CD=float(m.coef_[1]), smearing=float(np.mean(np.exp(y - m.predict(X)))), n=len(z), n_plots=groups.nunique())
    return (result, fold_id, fit)

def errors(truth, pred):
    a = np.asarray(truth, dtype=float)
    b = np.asarray(pred, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a = a[ok]
    b = b[ok]
    if not len(a):
        return dict(n=0, MAE=np.nan, RMSE=np.nan, Bias=np.nan, r=np.nan)
    e = b - a
    r = np.corrcoef(a, b)[0, 1] if len(a) > 1 and np.std(a) > 0 and (np.std(b) > 0) else np.nan
    return dict(n=len(a), MAE=float(np.mean(abs(e))), RMSE=float(np.sqrt(np.mean(e ** 2))), Bias=float(np.mean(e)), r=r)

def plot_biomass_totals(detections_all, truth_all, truth_coord, pairs, plot, carbon_fraction=0.5, carbon_alternative=0.47):

    def number(value):
        try:
            return float(value)
        except (ValueError, TypeError):
            return float('nan')

    def sum_valid(rows, column):
        values = [number(r.get(column)) for r in rows]
        return sum(values) if all((math.isfinite(v) for v in values)) else float('nan')
    area = number(plot.get('reported_area_m2'))
    if not math.isfinite(area) or area <= 0:
        area = number(plot.get('area_m2'))
    if not math.isfinite(area) or area <= 0:
        raise ValueError('양수 조사면적 필요')
    scale = lambda value: value / 1000 * 10000 / area
    by_id = {r['individualID']: r for r in truth_all}
    if len(by_id) != len(truth_all):
        raise ValueError('조사구 참조 개체 중복')
    det_by_id = {r['det_id']: r for r in detections_all}
    paired_truth = {r['individualID'] for r in pairs}
    paired_det = {r['det_id'] for r in pairs}
    coord_ids = {r['individualID'] for r in truth_coord}
    if not paired_truth <= by_id.keys() or not paired_det <= det_by_id.keys():
        raise ValueError('매칭 ID가 조사구 모집단에 없음')
    reference = sum_valid(truth_all, 'ref_AGB_kg')
    pred = sum_valid(detections_all, 'AGB_default_kg')
    coord_reference = sum_valid(truth_coord, 'ref_AGB_kg')
    missing_agb = sum((not math.isfinite(number(r.get('ref_AGB_kg'))) for r in truth_all))
    matched_ref = sum_valid([by_id[k] for k in paired_truth], 'ref_AGB_kg')
    matched_pred = sum_valid([det_by_id[k] for k in paired_det], 'AGB_default_kg')
    unmatched_pred = sum_valid([r for r in detections_all if r['det_id'] not in paired_det], 'AGB_default_kg')
    missed_ref = sum_valid([r for r in truth_all if r['individualID'] not in paired_truth], 'ref_AGB_kg')
    no_coord = sum_valid([r for r in truth_all if r['individualID'] not in coord_ids], 'ref_AGB_kg')
    missed_coord = sum_valid([r for r in truth_coord if r['individualID'] not in paired_truth], 'ref_AGB_kg')
    field_year = int(plot['field_year'])
    imagery_year = int(plot.get('imagery_year', field_year))
    row = dict(plotID=plot['plotID'], siteID=plot['siteID'], field_year=field_year, imagery_year=imagery_year, same_year=field_year == imagery_year, imagery_provisional=plot.get('imagery_provisional', False), area_m2=area, reported_area_m2=number(plot.get('reported_area_m2', area)), polygon_area_m2=number(plot.get('polygon_area_m2', area)), area_difference_m2=number(plot.get('polygon_area_m2', area)) - area, n_truth=len(truth_all), n_truth_coord=len(truth_coord), n_detected=len(detections_all), n_matched=len(pairs), n_field_before_coordinate_filter=len(truth_all), n_without_matching_coordinate=len(truth_all) - len(truth_coord), n_missing_ref_AGB=missing_agb, AGB_reference_coverage=1 - missing_agb / len(truth_all) if truth_all else 1.0, reference_complete=missing_agb == 0, coordinates_complete=len(truth_all) == len(truth_coord), AGB_est_t_ha=scale(pred), AGB_ref_t_ha=scale(reference), AGB_ref_coordonly_t_ha=scale(coord_reference), matched_error_t_ha=scale(matched_pred - matched_ref), unmatched_detection_t_ha=scale(unmatched_pred), missed_truth_t_ha=scale(missed_ref), unmapped_truth_t_ha=scale(no_coord), missed_mapped_truth_t_ha=scale(missed_coord))
    row['C_est_tC_ha'] = row['AGB_est_t_ha'] * carbon_fraction
    row['C_ref_tC_ha'] = row['AGB_ref_t_ha'] * carbon_fraction
    row['C_ref_coordonly_tC_ha'] = row['AGB_ref_coordonly_t_ha'] * carbon_fraction
    row['C_est_047_tC_ha'] = row['AGB_est_t_ha'] * carbon_alternative
    row['C_ref_047_tC_ha'] = row['AGB_ref_t_ha'] * carbon_alternative
    row['C_bias_tC_ha'] = row['C_est_tC_ha'] - row['C_ref_tC_ha']
    row['C_missed_tC_ha'] = row['missed_truth_t_ha'] * carbon_fraction
    row['estimate_reference_ratio'] = pred / reference if reference > 0 else float('nan')
    if math.isfinite(reference) and math.isfinite(pred):
        assert math.isclose(pred - reference, matched_pred - matched_ref + unmatched_pred - missed_ref, rel_tol=1e-09, abs_tol=1e-06)
    return row

def comparison_statistics(truth, pred):
    result = errors(truth, pred)
    a = np.asarray(truth, dtype=float)
    b = np.asarray(pred, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a = a[ok]
    b = b[ok]
    mean_ref = float(a.mean()) if len(a) else np.nan
    mean_est = float(b.mean()) if len(b) else np.nan
    result.update(mean_reference=mean_ref, mean_estimate=mean_est, estimate_reference_ratio=mean_est / mean_ref if mean_ref > 0 else np.nan)
    return result

def finalize_detection_frame(records, method):
    cols = DET_COLUMNS if 'DET_COLUMNS' in globals() else ['det_id', 'plotID', 'siteID', 'field_year', 'score', 'H_m', 'CD_m', 'geometry']
    frame = gpd.GeoDataFrame(records, geometry='geometry', crs=WORK_CRS) if len(records) else gpd.GeoDataFrame(columns=cols, geometry='geometry', crs=WORK_CRS)
    frame['method'] = method
    frame['DBH_jucker_cm'] = estimate_dbh(frame.H_m, frame.CD_m)
    frame['AGB_default_kg'] = estimate_agb(frame.DBH_jucker_cm)
    frame['C_default_kg'] = carbon_kg(frame.AGB_default_kg, CFG['carbon_fraction'])
    return frame

def deepforest_patch_candidates(detector, rgb, chm, plot, patch_size, cache_dir):
    saved = []
    for obj in [detector.model, getattr(detector.model, 'model', None)]:
        if obj is not None and hasattr(obj, 'score_thresh'):
            saved.append((obj, 'score_thresh', obj.score_thresh))
            obj.score_thresh = 0.05
    configs = []
    for obj in [detector.config, getattr(detector.config, 'retinanet', None)]:
        if obj is None:
            continue
        try:
            old = obj['score_thresh']
            obj['score_thresh'] = 0.05
            configs.append((obj, old))
        except (KeyError, TypeError, AttributeError):
            pass
    if not saved and (not configs):
        raise RuntimeError('DeepForest 내부 score_thresh를 0.05로 설정하지 못했습니다. 낮은 점수 탐색을 숨겨 진행하지 않습니다.')
    oldcfg = CFG.copy()
    oldout = globals()['OUT']
    directory = Path(cache_dir) / f'patch_{patch_size}'
    directory.mkdir(parents=True, exist_ok=True)
    fingerprint = dict(rgb=str(rgb), rgb_bytes=Path(rgb).stat().st_size, patch_size=int(patch_size), patch_overlap=CFG['patch_overlap'], nms_iou=CFG['nms_iou'], min_score=0.1, model_version=importlib.metadata.version('deepforest'), model=CFG['model_name'], revision=CFG['model_revision'])
    stamp = directory / f'{plot.plotID}_input.json'
    rawpath = directory / f'{plot.plotID}_raw_predictions.csv'
    if rawpath.exists() and (not stamp.exists() or json.loads(stamp.read_text()) != fingerprint):
        rawpath.unlink()
    try:
        CFG['patch_size'] = int(patch_size)
        CFG['score_threshold'] = 0.1
        globals()['OUT'] = directory
        rows = predict_crowns(detector, rgb, chm, plot)
        stamp.write_text(json.dumps(fingerprint, indent=2))
    finally:
        CFG.clear()
        CFG.update(oldcfg)
        globals()['OUT'] = oldout
        for obj, key, value in saved:
            setattr(obj, key, value)
        for obj, value in configs:
            obj['score_thresh'] = value
    for row in rows:
        row['det_id'] = f"{plot.plotID}_p{patch_size}_{row['det_id'].split('_')[-1]}"
    return finalize_detection_frame(rows, f'DF_patch_{patch_size}')

def variable_window_peaks(height, valid, resolution_m, diameter_a=2.51, diameter_b=0.00901, min_height=3.0):
    mask = valid & (height >= min_height)
    cross = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
    if resolution_m > (diameter_a + diameter_b * min_height ** 2) / 2:
        raise ValueError('해상도가 최소 가변 창 반지름보다 큽니다. 더 세밀한 CHM 또는 현지 창 계수를 사용하세요.')
    maxima = mask & (height >= maximum_filter(np.where(valid, height, -np.inf), footprint=cross, mode='constant', cval=-np.inf))
    components, count = ndi_label(maxima, structure=cross.astype(int))
    candidates = []
    for k in range(1, count + 1):
        coords = np.argwhere(components == k)
        order = np.lexsort((coords[:, 1], coords[:, 0], -height[coords[:, 0], coords[:, 1]]))
        r, c = coords[order[0]]
        H = float(height[r, c])
        radius = (diameter_a + diameter_b * H * H) / 2
        pixels = int(np.ceil(radius / resolution_m))
        r0 = max(0, r - pixels)
        r1 = min(height.shape[0], r + pixels + 1)
        c0 = max(0, c - pixels)
        c1 = min(height.shape[1], c + pixels + 1)
        rr, cc = np.ogrid[r0:r1, c0:c1]
        circle = ((rr - r) * resolution_m) ** 2 + ((cc - c) * resolution_m) ** 2 <= radius ** 2
        neighborhood = height[r0:r1, c0:c1][circle & valid[r0:r1, c0:c1]]
        if len(neighborhood) and H >= float(neighborhood.max()) - 1e-07:
            candidates.append((int(r), int(c), H, float(radius * 2)))
    return sorted(candidates, key=lambda x: (-x[2], x[0], x[1]))

def chm_watershed_crowns(chm_path, plot, sigma_m=1.0):
    with rasterio.open(chm_path) as ds:
        raw = ds.read(1, masked=True)
        transform = ds.transform
        source_crs = ds.crs
        resolution = float(abs(transform.a))
        if not np.isclose(abs(transform.a), abs(transform.e)):
            raise ValueError('CHM 정사각형 픽셀 필요')
        rawdata = np.asarray(raw.astype(float).filled(np.nan), dtype=float)
        valid = np.isfinite(rawdata) & (rawdata >= 0)
        weight = gaussian_filter(valid.astype(float), sigma=sigma_m / resolution)
        smooth = gaussian_filter(np.where(valid, rawdata, 0), sigma=sigma_m / resolution) / np.maximum(weight, 1e-12)
        peaks = variable_window_peaks(smooth, valid, resolution, CFG['chm_window_a'], CFG['chm_window_b'], CFG['min_height_m'])
        markers = np.zeros(rawdata.shape, dtype=np.int32)
        accepted = []
        for r, c, H, diameter in peaks:
            if rawdata[r, c] < CFG['min_height_m']:
                continue
            k = len(accepted) + 1
            markers[r, c] = k
            accepted.append((k, r, c, float(rawdata[r, c]), diameter))
        if not accepted:
            return finalize_detection_frame([], 'B_CHM')
        segments = watershed(-smooth, markers, mask=valid & (rawdata >= CFG['min_height_m']))
        crown_parts = {}
        for geom, value in raster_shapes(segments.astype('int32'), mask=segments > 0, transform=transform):
            crown_parts.setdefault(int(value), []).append(shape(geom))
    rows = []
    for k, r, c, H, diameter in accepted:
        xy = transform * (c + 0.5, r + 0.5)
        point = gpd.GeoSeries([Point(*xy)], crs=source_crs).to_crs(WORK_CRS).iloc[0]
        if not plot.geometry.covers(point):
            continue
        polygon = gpd.GeoSeries([unary_union(crown_parts[k])], crs=source_crs).to_crs(WORK_CRS).iloc[0]
        area = float(polygon.area)
        if area <= 0:
            continue
        rows.append(dict(det_id=f'{plot.plotID}_CHM_{k}', plotID=plot.plotID, siteID=plot.siteID, field_year=plot.field_year, imagery_year=plot.imagery_year, score=1.0, H_m=H, CD_m=2 * np.sqrt(area / np.pi), crown_area_m2=area, crown_geometry_wkt=polygon.wkt, peak_x=point.x, peak_y=point.y, window_diameter_m=diameter, geometry=polygon.envelope))
    return finalize_detection_frame(rows, 'B_CHM')

def combine_boxes_and_peaks(boxes, chm, dedup_m=1.0):
    rows = []
    used = set()
    for _, boxrow in boxes.sort_values(['score', 'det_id'], ascending=[False, True]).iterrows():
        candidates = chm[chm.apply(lambda p: boxrow.geometry.covers(Point(p.peak_x, p.peak_y)), axis=1)] if len(chm) else chm
        had_seed = len(candidates) > 0
        candidates = candidates[~candidates.det_id.isin(used)]
        if had_seed and candidates.empty:
            continue
        row = boxrow.to_dict()
        row['det_id'] = 'C_' + str(row['det_id'])
        if len(candidates):
            peak = candidates.sort_values(['H_m', 'det_id'], ascending=[False, True]).iloc[0]
            used.add(peak.det_id)
            row.update(H_m=peak.H_m, peak_x=peak.peak_x, peak_y=peak.peak_y)
        else:
            row.update(peak_x=boxrow.geometry.centroid.x, peak_y=boxrow.geometry.centroid.y)
        point = Point(row['peak_x'], row['peak_y'])
        if any((point.distance(Point(x['peak_x'], x['peak_y'])) < dedup_m for x in rows)):
            continue
        row['fusion_source'] = 'DeepForest_box'
        rows.append(row)
    for _, p in chm.iterrows():
        point = Point(p.peak_x, p.peak_y)
        if any((g.covers(point) for g in boxes.geometry)):
            continue
        if any((point.distance(Point(x['peak_x'], x['peak_y'])) < dedup_m for x in rows)):
            continue
        row = p.to_dict()
        row['det_id'] = 'C_' + row['det_id']
        row['fusion_source'] = 'outside_CHM_peak'
        rows.append(row)
    return finalize_detection_frame(rows, 'C_fusion')

def detection_counts(detections, truth, plot_ids):
    tp = fp = fn = 0
    for pid in plot_ids:
        d = detections[detections.plotID == pid]
        f = truth[truth.plotID == pid]
        hits = len(match_one_to_one(d, f))
        tp += hits
        fp += len(d) - hits
        fn += len(f) - hits
    return dict(TP=tp, FP=fp, FN=fn, precision=tp / (tp + fp) if tp + fp else np.nan, recall=tp / (tp + fn) if tp + fn else np.nan, F1=2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else np.nan)

def nested_group_parameter_selection(candidates, plot_table, truth, n_splits=5):
    ids = np.array(sorted(plot_table.plotID.unique()))
    outer = GroupKFold(n_splits=min(n_splits, len(ids)))
    predictions = []
    selections = []
    inner_rows = []
    split_rows = []
    keys = sorted(candidates)
    for fold, (train, test) in enumerate(outer.split(ids, groups=ids)):
        trainids = ids[train]
        testids = ids[test]
        assert set(trainids).isdisjoint(testids)
        for pid in trainids:
            split_rows.append(dict(outer_fold=fold, plotID=pid, role='train'))
        for pid in testids:
            split_rows.append(dict(outer_fold=fold, plotID=pid, role='test'))
        inner = GroupKFold(n_splits=min(CFG['search_inner_folds'], len(trainids)))
        scores = {key: dict(TP=0, FP=0, FN=0) for key in keys}
        for inner_fold, (itr, ival) in enumerate(inner.split(trainids, groups=trainids)):
            assert set(trainids[itr]).isdisjoint(trainids[ival])
            assert set(trainids[ival]).isdisjoint(testids)
            for key in keys:
                result = detection_counts(candidates[key], truth, trainids[ival])
                inner_rows.append(dict(outer_fold=fold, inner_fold=inner_fold, patch_size=key[0], score_threshold=key[1], **result))
                for k in ['TP', 'FP', 'FN']:
                    scores[key][k] += result[k]

        def ranking(key):
            counts = scores[key]
            den = 2 * counts['TP'] + counts['FP'] + counts['FN']
            f1 = 2 * counts['TP'] / den if den else -1
            return (-f1, key[0], key[1])
        best = min(keys, key=ranking)
        chosen = candidates[best][candidates[best].plotID.isin(testids)].copy()
        chosen['outer_fold'] = fold
        chosen['selected_patch_size'] = best[0]
        chosen['selected_score_threshold'] = best[1]
        predictions.append(chosen)
        selections.append(dict(outer_fold=fold, patch_size=best[0], score_threshold=best[1], train_plots=';'.join(trainids), test_plots=';'.join(testids), **detection_counts(chosen, truth, testids)))
    oof = gpd.GeoDataFrame(pd.concat(predictions, ignore_index=True), geometry='geometry', crs=WORK_CRS)
    deployment_rows = []
    for key in keys:
        deployment_rows.append(dict(patch_size=key[0], score_threshold=key[1], **detection_counts(candidates[key], truth, ids)))
    deployment = pd.DataFrame(deployment_rows).sort_values(['F1', 'patch_size', 'score_threshold'], ascending=[False, True, True])
    return (oof, pd.DataFrame(selections), pd.DataFrame(inner_rows), pd.DataFrame(split_rows), deployment)

def evaluate_detection_method(det, plot_table, all_truth, coordinate_truth):
    metric_rows = []
    pair_tables = {}
    joined = None
    for population in ['visible', 'all']:
        pairs = []
        for _, p in plot_table.iterrows():
            d = det[det.plotID == p.plotID]
            f = coordinate_truth[coordinate_truth.plotID == p.plotID]
            if population == 'visible':
                f = f[f.visible]
            match = match_one_to_one(d, f)
            pairs.append(match)
            metric_rows.append(dict(plotID=p.plotID, siteID=p.siteID, field_year=p.field_year, population=population, **detection_counts(d, f, [p.plotID])))
        pair_tables[population] = pd.concat(pairs, ignore_index=True)
    visible = pair_tables['visible']
    joined = visible.merge(det.drop(columns='geometry'), on=['det_id', 'plotID'], validate='one_to_one')
    joined = joined.merge(coordinate_truth.drop(columns=['geometry', 'plotID', 'siteID', 'field_year']), on='individualID', validate='one_to_one')
    joined['DBH_regression_oof_cm'], joined['cv_fold'], _ = grouped_dbh_regression(joined, CFG['cv_folds'])
    joined['AGB_regression_oof_kg'] = estimate_agb(joined.DBH_regression_oof_cm)
    joined['AGB_known_species_kg'] = estimate_agb(joined.DBH_jucker_cm, joined.species_beta0, joined.species_beta1)
    joined['species_effect_kg'] = joined.AGB_known_species_kg - joined.AGB_default_kg
    joined['species_effect_pct'] = 100 * (joined.AGB_known_species_kg / joined.AGB_default_kg - 1)
    plotrows = []
    for _, p in plot_table.iterrows():
        ds = det[det.plotID == p.plotID]
        fa = all_truth[all_truth.plotID == p.plotID]
        fc = coordinate_truth[coordinate_truth.plotID == p.plotID]
        pairs = pair_tables['all']
        pairs = pairs[pairs.plotID == p.plotID]
        plotrows.append(plot_biomass_totals(ds.to_dict('records'), fa.to_dict('records'), fc.to_dict('records'), pairs.to_dict('records'), p.to_dict(), CFG['carbon_fraction'], CFG['carbon_alternative']))
    return dict(detections=det, matched=joined, metrics=pd.DataFrame(metric_rows), pairs=pair_tables, plots=pd.DataFrame(plotrows))

def method_comparison_rows(method_results, requested_sites):
    rows = []
    for method, result in method_results.items():
        for site in ['ALL'] + list(requested_sites):
            pc = result['plots']
            mt = result['matched']
            dm = result['metrics']
            if site != 'ALL':
                pc = pc[pc.siteID == site]
                mt = mt[mt.siteID == site]
                dm = dm[dm.siteID == site]
            stats = comparison_statistics(pc.C_ref_tC_ha, pc.C_est_tC_ha)
            row = dict(method=method, siteID=site, n_plots=len(pc), n_field=int(pc.n_truth.sum()), n_detected=int(pc.n_detected.sum()), n_carbon_plots=stats['n'], C_est_mean_tC_ha=stats['mean_estimate'], C_ref_mean_tC_ha=stats['mean_reference'], estimate_reference_ratio=stats['estimate_reference_ratio'], carbon_r=stats['r'], carbon_RMSE_tC_ha=stats['RMSE'], carbon_MAE_tC_ha=stats['MAE'])
            for pop in ['visible', 'all']:
                q = dm[dm.population == pop]
                tp = q.TP.sum()
                fp = q.FP.sum()
                fn = q.FN.sum()
                row['F1_' + pop] = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else np.nan
                row['precision_' + pop] = tp / (tp + fp) if tp + fp else np.nan
                row['recall_' + pop] = tp / (tp + fn) if tp + fn else np.nan
            for name, a, b in [('height_Bias_m', 'field_H_m', 'H_m'), ('crown_Bias_m', 'field_CD_m', 'CD_m'), ('DBH_Bias_cm', 'field_DBH_cm', 'DBH_jucker_cm'), ('individual_AGB_Bias_kg', 'ref_AGB_kg', 'AGB_default_kg')]:
                row[name] = errors(mt[a], mt[b])['Bias']
            rows.append(row)
    return pd.DataFrame(rows)

def detection_algorithm_smoke_checks():
    from types import SimpleNamespace
    yy, xx = np.mgrid[:64, :64]
    canopy = 20 * np.exp(-((xx - 18) ** 2 + (yy - 32) ** 2) / 8) + 30 * np.exp(-((xx - 46) ** 2 + (yy - 32) ** 2) / 8)
    path = OUT / 'synthetic_CHM_smoke.tif'
    with rasterio.open(path, 'w', driver='GTiff', width=64, height=64, count=1, dtype='float32', crs=WORK_CRS, transform=rasterio.transform.from_origin(0, 64, 1, 1), nodata=-9999) as ds:
        ds.write(canopy.astype('float32'), 1)
    p = SimpleNamespace(plotID='SMOKE', siteID='SMOKE', field_year=2021, imagery_year=2021, geometry=box(0, 0, 64, 64))
    result = chm_watershed_crowns(path, p, 1.0)
    assert len(result) == 2, '합성 두 나무의 꼭대기 2개를 찾아야 합니다'
    assert np.allclose(sorted(result.H_m), [20, 30], atol=0.01)
    assert np.allclose(result.CD_m, 2 * np.sqrt(result.crown_area_m2 / np.pi))
    assert np.isclose(result.crown_area_m2.sum(), np.count_nonzero(canopy >= CFG['min_height_m']))
    first_peak = result.iloc[0].to_dict()
    duplicate = dict(first_peak)
    duplicate['det_id'] = 'duplicate_box'
    duplicate['score'] = 0.8
    first_peak['score'] = 0.9
    boxes = gpd.GeoDataFrame([first_peak, duplicate], geometry='geometry', crs=WORK_CRS)
    fused = combine_boxes_and_peaks(boxes, result, 1.0)
    assert len(fused) == 2, '동일 꼭대기 중복 상자는 제거하고 상자 밖 두 번째 꼭대기는 추가해야 합니다'
    ids = ['P0', 'P1', 'P2', 'P3']
    drows = []
    frows = []
    for i, pid in enumerate(ids):
        x = i * 20
        drows.append(dict(det_id='d' + pid, plotID=pid, score=0.8, H_m=20, CD_m=5, geometry=box(x - 3, -3, x + 3, 3)))
        frows.append(dict(individualID='f' + pid, plotID=pid, stemDiameter=30, geometry=Point(x, 0)))
    ds = gpd.GeoDataFrame(drows, geometry='geometry', crs=WORK_CRS)
    fs = gpd.GeoDataFrame(frows, geometry='geometry', crs=WORK_CRS)
    plotdf = pd.DataFrame({'plotID': ids})
    grid = {(300, 0.1): ds, (400, 0.2): ds.iloc[:2].copy()}
    _, first, _, _, _ = nested_group_parameter_selection(grid, plotdf, fs, n_splits=2)
    changed = fs.copy()
    changed.loc[changed.plotID == 'P0', 'geometry'] = Point(9999, 9999)
    _, second, _, _, _ = nested_group_parameter_selection(grid, plotdf, changed, n_splits=2)
    fold = first[first.test_plots.str.split(';').apply(lambda z: 'P0' in z)].index[0]
    assert first.loc[fold, ['patch_size', 'score_threshold']].equals(second.loc[fold, ['patch_size', 'score_threshold']]), '바깥 평가 정답이 선택에 새어 들어갔습니다'
    for _, r in first.iterrows():
        assert set(r.train_plots.split(';')).isdisjoint(r.test_plots.split(';'))
    path.unlink(missing_ok=True)
    print('합성 CHM 꼭대기·watershed 면적·CD, GroupKFold 정답 격리 검증 통과')

def plot_carbon_scatter(plot_table, output_folder):
    for (year, image_year), cohort in plot_table.groupby(['field_year', 'imagery_year']):
        z = cohort.replace([np.inf, -np.inf], np.nan).dropna(subset=['C_ref_tC_ha', 'C_est_tC_ha'])
        fig, ax = plt.subplots(figsize=(6, 6))
        for site, q in z.groupby('siteID'):
            ax.scatter(q.C_ref_tC_ha, q.C_est_tC_ha, s=35, alpha=0.8, label=site, color={'GRSM': 'tab:blue', 'ORNL': 'tab:orange'}.get(site))
        if len(z):
            hi = max(z[['C_ref_tC_ha', 'C_est_tC_ha']].max()) * 1.05
            hi = max(hi, 1.0)
            ax.plot([0, hi], [0, hi], 'k--', label='1:1')
            ax.set_xlim(0, hi)
            ax.set_ylim(0, hi)
        ax.set(xlabel='Field-based reference carbon (tC/ha)', ylabel='Estimated carbon (tC/ha)', title=f'Plot carbon: field {year} / imagery {image_year} (n={len(z)})')
        ax.legend()
        ax.set_aspect('equal', adjustable='box')
        filename = f'plot_carbon_{year}.png' if year == image_year else f'plot_carbon_field{year}_imagery{image_year}.png'
        fig.savefig(Path(output_folder) / filename, dpi=170, bbox_inches='tight')
        plt.close(fig)

def build_site_summary(plot_table, detection_metrics, matched_table, requested_sites):
    rows = []
    for site in requested_sites:
        site_plots = plot_table[plot_table.siteID == site]
        if site_plots.empty:
            rows.append(dict(siteID=site, n_plots=0, comparison_note='분석 결과 없음'))
            continue
        for (field_year, image_year), q in site_plots.groupby(['field_year', 'imagery_year']):
            ids = set(q.plotID)
            dm = detection_metrics[detection_metrics.plotID.isin(ids)]
            mt = matched_table[matched_table.plotID.isin(ids)]
            row = dict(siteID=site, field_year=field_year, imagery_year=image_year, comparison_note='영상 연도 = 실측 연도' if field_year == image_year else f'{site}은 영상 연도 ≠ 실측 연도', n_plots=len(q), n_field_trees=int(q.n_truth.sum()), n_coord_trees=int(q.n_truth_coord.sum()), n_detections=int(q.n_detected.sum()))
            for population, label in [('visible', 'F1_visible'), ('all', 'F1_all')]:
                z = dm[dm.population == population]
                tp = z.TP.sum()
                fp = z.FP.sum()
                fn = z.FN.sum()
                row[label] = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else np.nan
            for label, refcol, predcol in [('height_Bias_m', 'field_H_m', 'H_m'), ('DBH_Bias_cm', 'field_DBH_cm', 'DBH_jucker_cm'), ('individual_AGB_Bias_kg', 'ref_AGB_kg', 'AGB_default_kg')]:
                row[label] = errors(mt[refcol], mt[predcol])['Bias']
            stats = comparison_statistics(q.C_ref_tC_ha, q.C_est_tC_ha)
            row.update(n_carbon_plots=stats['n'], C_est_mean_tC_ha=stats['mean_estimate'], C_ref_mean_tC_ha=stats['mean_reference'], estimate_reference_ratio=stats['estimate_reference_ratio'])
            rows.append(row)
    return pd.DataFrame(rows)
