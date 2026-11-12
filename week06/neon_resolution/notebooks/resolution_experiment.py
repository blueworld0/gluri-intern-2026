"""Colab-only extension; embedded in the notebook by resolution_extension.py.

Reuses the completed third run. Never downloads imagery or reruns old methods.
Dependencies/functions/globals come from the original notebook definition cells.
"""
import shutil
import zipfile
from contextlib import contextmanager
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
from rasterio.io import MemoryFile
from shapely import wkt

RESOLUTION_SCHEMA = 'resolution-v1-average-native-bilinear-10cm'
RESOLUTION_GSDS = (10, 25, 51)
RESOLUTION_THRESHOLDS = (0.1, 0.2, 0.3, 0.4, 0.5)
RESOLUTION_PATCH = 600
RESOLUTION_FIXED_SCORE = 0.2
RESOLUTION_OVERLAYS = ('GRSM_016', 'ORNL_006', 'GRSM_051')


def resolution_sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def resolution_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.part')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    tmp.replace(path)


def resolution_read_geo(path):
    frame = pd.read_csv(path)
    if 'geometry' not in frame:
        raise ValueError(f'WKT geometry column missing: {path}')
    frame['geometry'] = frame.geometry.map(wkt.loads)
    return gpd.GeoDataFrame(frame, geometry='geometry', crs=WORK_CRS)


def resolution_bool(series):
    return series.map(lambda v: str(v).strip().lower() in ('true', '1'))


def resolution_required_files():
    files = ['config.json', 'versions.json', 'selected_detection_method.json',
             'comparison_cohort.csv', 'sampled_footprints.gpkg', 'raster_manifest.csv',
             'field_trees_all.csv', 'field_trees.csv', 'method_comparison.csv',
             'A_outer_groupkfold.csv', 'A_inner_grid_scores.csv',
             'A_fold_plot_assignments.csv', 'A_deployment_grid.csv', 'report.md']
    for method in ('baseline_015', 'A_DF_nested_OOF', 'B_CHM', 'C_fusion_nested_OOF'):
        files += [f'methods/{method}/{f}' for f in ('detected_trees.csv', 'matched_trees_visible.csv',
                  'detection_metrics.csv', 'plot_comparison.csv', 'matches_visible.csv', 'matches_all.csv')]
    return files


def resolution_find_base(root, explicit=None):
    required = resolution_required_files()
    expected = set(CFG['comparison_plot_ids'])

    def valid(path):
        if not all((path / f).is_file() for f in required):
            return False
        selection = json.loads((path / 'selected_detection_method.json').read_text())
        cohort = pd.read_csv(path / 'comparison_cohort.csv')
        return (selection['default_method'] == 'A_DF_nested_OOF'
                and selection['deployment_DF_patch_size'] == 600
                and np.isclose(selection['deployment_DF_score_threshold'], 0.2)
                and len(cohort) == 23 and set(cohort.plotID) == expected)

    if explicit:
        chosen = Path(explicit)
        if not valid(chosen):
            raise FileNotFoundError(f'3차 완료 결과/필수 파일/23개 cohort 확인 필요: {chosen}')
        return chosen
    folders = [root / 'outputs', root / 'outputs_v3'] + sorted(p for p in (root / 'outputs').glob('*') if p.is_dir())
    # An original result ZIP may have been unpacked on an earlier run. Keep using
    # that immutable snapshot even after the final ZIP is replaced by this run.
    folders += [p.parent for p in (root / 'data/cache/resolution_base').rglob('selected_detection_method.json')]
    candidates = [p for p in folders if valid(p)]
    # Prefer original result folders, not copies produced by this extension.
    originals = [p for p in candidates if not (p / 'resolution_provenance.json').exists()]
    candidates = originals or candidates
    if not candidates and (root / 'outputs/neon_results.zip').is_file():
        archive = root / 'outputs/neon_results.zip'
        dest = root / 'data/cache/resolution_base' / resolution_sha(archive)[:16]
        dest.mkdir(parents=True, exist_ok=True)
        if not (dest / '.complete').exists():
            with zipfile.ZipFile(archive) as z:
                for member in z.infolist():
                    target = (dest / member.filename).resolve()
                    if not target.is_relative_to(dest.resolve()):
                        raise ValueError('Unsafe member in local result ZIP')
                z.extractall(dest)
            (dest / '.complete').touch()
        candidates = [p for p in [dest] + [f.parent for f in dest.rglob('selected_detection_method.json')] if valid(p)]
    if not candidates:
        raise FileNotFoundError('Drive의 3차 결과 폴더/outputs/neon_results.zip을 찾지 못했습니다. '
                                'RESOLUTION_BASE_DIR에 완료된 결과 폴더를 지정하세요. 다운로드/기존 방법 재계산은 하지 않습니다.')
    signatures = {tuple(resolution_sha(p / f) for f in
                       ('method_comparison.csv', 'field_trees_all.csv', 'field_trees.csv', 'A_outer_groupkfold.csv')) for p in candidates}
    if len(signatures) != 1:
        raise RuntimeError('서로 다른 완료 결과가 여러 개입니다. RESOLUTION_BASE_DIR로 지난주 3차 폴더를 지정하세요.')
    return sorted(set(candidates), key=lambda p: (len(p.parts), str(p)))[0]


def resolution_crop_path(recorded, plot_id, kind, year):
    recorded = Path(recorded)
    if recorded.is_file():
        return recorded
    relocated = ROOT / 'data/crops' / recorded.name
    if relocated.is_file():
        return relocated
    # Never substitute a different crop with a similar name or download a tile.
    raise FileNotFoundError(f'{plot_id} {year} {kind} 기존 crop 없음: {recorded}; '
                            f'Drive data/crops/{recorded.name}를 확인하세요.')


def resolution_load_state():
    global CFG, OUT, WORK_CRS, default, plots, field_all, field_geo, RASTERS
    global RESOLUTION_BASE, RESOLUTION_CACHE, RESOLUTION_BASE_SIGNATURE, RESOLUTION_INPUTS
    global RESOLUTION_BASE_METHODS, RESOLUTION_BASE_SUMMARY, RESOLUTION_SELECTION, RUN_ID
    global RESOLUTION_DETECTOR
    RESOLUTION_DETECTOR = None
    RESOLUTION_BASE = resolution_find_base(ROOT, RESOLUTION_BASE_DIR)
    original_cfg = json.loads((RESOLUTION_BASE / 'config.json').read_text())
    required_cfg = dict(jucker=dict(a=0.367, b=0.887, g=0.0, biome=11), carbon_fraction=0.5,
                        min_height_m=3.0, match_size='stemDiameter', cv_folds=5, search_inner_folds=4)
    for key, value in required_cfg.items():
        if original_cfg[key] != value:
            raise ValueError(f'기존 설정 불일치: {key}: {original_cfg[key]} != {value}')
    assert np.isclose(original_cfg['agb_default']['beta0'], -2.2118)
    assert np.isclose(original_cfg['agb_default']['beta1'], 2.4133)
    assert original_cfg['aop_include_provisional'] is True
    versions = json.loads((RESOLUTION_BASE / 'versions.json').read_text())
    if importlib.metadata.version('deepforest') != versions['deepforest']:
        raise RuntimeError(f'DeepForest {versions["deepforest"]} 필요. 첫 설치 셀/런타임 재시작을 확인하세요.')
    CFG = original_cfg
    default = CFG['agb_default']
    WORK_CRS = CRS.from_epsg(32617)
    cohort = pd.read_csv(RESOLUTION_BASE / 'comparison_cohort.csv')
    if cohort.siteID.value_counts().to_dict() != {'GRSM': 15, 'ORNL': 8}:
        raise ValueError('GRSM15 + ORNL8 이외의 cohort')
    assert set(cohort.plotID) == set(CFG['comparison_plot_ids'])
    plots = gpd.read_file(RESOLUTION_BASE / 'sampled_footprints.gpkg').to_crs(WORK_CRS)
    plots = plots[plots.plotID.isin(cohort.plotID)].copy()
    for col in ('siteID', 'field_year', 'imagery_year'):
        plots[col] = plots.plotID.map(cohort.set_index('plotID')[col])
    assert plots.plotID.is_unique and len(plots) == 23
    field_all = pd.read_csv(RESOLUTION_BASE / 'field_trees_all.csv')
    field_geo = resolution_read_geo(RESOLUTION_BASE / 'field_trees.csv')
    field_geo['visible'] = resolution_bool(field_geo.visible)
    for name, frame in [('field_all', field_all), ('field_geo', field_geo)]:
        if not set(frame.plotID).issubset(set(cohort.plotID)):
            raise ValueError(f'{name}에 다른 조사구 포함')
        assert not frame.individualID.duplicated().any()
    manifest = pd.read_csv(RESOLUTION_BASE / 'raster_manifest.csv').set_index('plotID')
    assert set(cohort.plotID).issubset(set(manifest.index))
    plots['imagery_provisional'] = plots.plotID.map(resolution_bool(manifest.imagery_provisional))
    assert plots.loc[plots.siteID == 'ORNL', 'imagery_provisional'].all()
    RASTERS = {}; inputs = []
    for _, p in plots.iterrows():
        row = manifest.loc[p.plotID]
        paths = tuple(resolution_crop_path(row[k], p.plotID, k, p.imagery_year) for k in ('rgb', 'chm'))
        RASTERS[p.plotID] = paths
        for kind, path, gsd in zip(('rgb', 'chm'), paths, (0.1, 1.0)):
            with rasterio.open(path) as ds:
                if not ds.crs or not ds.crs.is_projected or not np.allclose(ds.res, gsd, atol=1e-6):
                    raise ValueError(f'{kind} 원본 CRS/GSD 확인: {path}: {ds.crs}, {ds.res}')
                if kind == 'rgb' and (ds.count < 3 or ds.dtypes[0] != 'uint8'):
                    raise ValueError('원본 RGB uint8 3밴드 필요')
            inputs.append(dict(plotID=p.plotID, kind=kind, path=str(path), sha256=resolution_sha(path)))
    RESOLUTION_INPUTS = pd.DataFrame(inputs)
    original_hashes = {f: resolution_sha(RESOLUTION_BASE / f) for f in resolution_required_files()}
    spec = dict(schema=RESOLUTION_SCHEMA, code_sha=RESOLUTION_CODE_SHA,
                source_files=original_hashes, rasters=inputs,
                resolutions=list(RESOLUTION_GSDS), patch=600, thresholds=list(RESOLUTION_THRESHOLDS))
    RESOLUTION_BASE_SIGNATURE = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
    RUN_ID = 'resolution_' + RESOLUTION_BASE_SIGNATURE[:12]
    OUT = ROOT / 'outputs' / RUN_ID
    RESOLUTION_CACHE = ROOT / 'data/cache/resolution' / RESOLUTION_BASE_SIGNATURE[:16]
    OUT.mkdir(parents=True, exist_ok=True); RESOLUTION_CACHE.mkdir(parents=True, exist_ok=True)
    if OUT.resolve() == RESOLUTION_BASE.resolve():
        raise ValueError('기존 결과와 해상도 실험 결과 폴더를 분리해야 합니다.')
    # Preserve all old method files, reports and figures verbatim in the new ZIP.
    for source in list(RESOLUTION_BASE.rglob('*')):
        rel = source.relative_to(RESOLUTION_BASE)
        if (not source.is_file() or source.is_relative_to(OUT) or source.suffix.lower() == '.zip'
                or (len(rel.parts) > 1 and rel.parts[0] not in ('methods', 'figures'))):
            continue
        dest = OUT / rel
        if not dest.exists():
            dest.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, dest)
    RESOLUTION_INPUTS.to_csv(OUT / 'resolution_input_manifest.csv', index=False)
    resolution_json(OUT / 'resolution_provenance.json', dict(base_directory=str(RESOLUTION_BASE),
                    base_signature=RESOLUTION_BASE_SIGNATURE, schema=RESOLUTION_SCHEMA,
                    source_files=original_hashes, settings=spec, no_imagery_download=True,
                    old_methods_recomputed=False))
    RESOLUTION_BASE_METHODS = {}
    for name in ('baseline_015', 'A_DF_nested_OOF', 'B_CHM', 'C_fusion_nested_OOF'):
        folder = RESOLUTION_BASE / 'methods' / name
        RESOLUTION_BASE_METHODS[name] = dict(detections=resolution_read_geo(folder / 'detected_trees.csv'),
            matched=pd.read_csv(folder / 'matched_trees_visible.csv'),
            metrics=pd.read_csv(folder / 'detection_metrics.csv'), plots=pd.read_csv(folder / 'plot_comparison.csv'),
            pairs={pop: pd.read_csv(folder / f'matches_{pop}.csv') for pop in ('visible', 'all')})
    RESOLUTION_BASE_SUMMARY = pd.read_csv(RESOLUTION_BASE / 'method_comparison.csv')
    RESOLUTION_SELECTION = json.loads((RESOLUTION_BASE / 'selected_detection_method.json').read_text())
    if (RESOLUTION_BASE / 'excluded_taxon_coefficients.csv').exists():
        conflict = pd.read_csv(RESOLUTION_BASE / 'excluded_taxon_coefficients.csv')
        if conflict.astype(str).apply(lambda s: s.str.contains('QUDO')).any().any():
            warnings.warn('기존 QUDO 계수 충돌: 해당 수종 계수 제외 상태를 유지합니다. 전체 참조 나무는 제거하지 않습니다.')
    print('기존 3차 결과:', RESOLUTION_BASE)
    print('해상도 실험 결과:', OUT)
    print('캐시:', RESOLUTION_CACHE)
    print('조사구:', plots.siteID.value_counts().to_dict(), '전체 참조 나무:', len(field_all), '좌표 나무:', len(field_geo))


def resolution_grid(width, height, xres, yres, requested_m):
    """Integer dimensions preserve the original extent; record actual GSD."""
    w = max(1, int(round(width * xres / requested_m)))
    h = max(1, int(round(height * yres / requested_m)))
    return w, h


def resolution_rgb(source, requested_cm, cache_directory):
    source = Path(source)
    with rasterio.open(source) as ds:
        profile = ds.profile.copy(); shape0 = (ds.height, ds.width)
        transform0 = ds.transform; crs0 = ds.crs; bounds0 = tuple(ds.bounds)
        if not np.isclose(transform0.b, 0) or not np.isclose(transform0.d, 0):
            raise ValueError('북향 정렬 원본 RGB 필요')
        w, h = resolution_grid(ds.width, ds.height, *ds.res, requested_cm / 100)
        low_transform = transform0 * rasterio.Affine.scale(ds.width / w, ds.height / h)
        actual = (abs(low_transform.a), abs(low_transform.e))
    metadata = dict(schema=RESOLUTION_SCHEMA, source_sha256=resolution_sha(source),
        nominal_cm=requested_cm, coarse_width=w, coarse_height=h,
        actual_gsd_x_m=actual[0], actual_gsd_y_m=actual[1], bounds=list(bounds0),
        coarse_transform=list(low_transform), input_transform=list(transform0),
        input_width=shape0[1], input_height=shape0[0], crs=crs0.to_string(),
        downsampling='rasterio.Resampling.average', upsampling='rasterio.Resampling.bilinear',
        rasterio_version=rasterio.__version__, gdal_version=rasterio.__gdal_version__)
    cache_directory = Path(cache_directory); cache_directory.mkdir(parents=True, exist_ok=True)
    coarse = cache_directory / 'coarse_native.tif'; expanded = cache_directory / 'model_input_10cm.tif'
    stamp = cache_directory / 'resampling.json'
    if not (coarse.is_file() and expanded.is_file() and stamp.is_file() and json.loads(stamp.read_text()) == metadata):
        with rasterio.open(source) as ds:
            # Float average avoids a premature uint8 rounding before bilinear expansion.
            low = ds.read([1, 2, 3], out_shape=(3, h, w), out_dtype='float32',
                          masked=True, resampling=Resampling.average)
        low_profile = profile.copy()
        low_profile.update(driver='GTiff', count=3, dtype='float32', nodata=np.nan,
                           height=h, width=w, transform=low_transform, compress='deflate')
        tmp = coarse.with_name(coarse.stem + '.part.tif')
        with rasterio.open(tmp, 'w', **low_profile) as ds:
            ds.write(low.filled(np.nan)); ds.update_tags(nominal_gsd_cm=str(requested_cm), experiment=RESOLUTION_SCHEMA)
        tmp.replace(coarse)
        with rasterio.open(coarse) as ds:
            expanded_data = ds.read(out_shape=(3, *shape0), resampling=Resampling.bilinear, masked=True)
        good = ~np.ma.getmaskarray(expanded_data).any(axis=0)
        expanded_data = np.clip(np.rint(expanded_data.filled(0)), 0, 255).astype('uint8')
        model_profile = profile.copy(); model_profile.update(driver='GTiff', count=3, dtype='uint8', compress='deflate')
        tmp = expanded.with_name(expanded.stem + '.part.tif')
        with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
            with rasterio.open(tmp, 'w', **model_profile) as ds:
                ds.write(expanded_data); ds.write_mask(good.astype('uint8') * 255)
                ds.update_tags(effective_nominal_gsd_cm=str(requested_cm), downsample='average', upsample='bilinear')
        tmp.replace(expanded)
        resolution_json(stamp, metadata)
    # Always read the files themselves, including on a cache hit.
    with rasterio.open(coarse) as low, rasterio.open(expanded) as model_input:
        assert low.crs == model_input.crs == crs0
        assert np.allclose(low.bounds, bounds0, rtol=0, atol=1e-8)
        assert np.allclose(low.res, actual, rtol=0, atol=1e-10)
        assert model_input.transform == transform0 and model_input.shape == shape0
        assert tuple(model_input.bounds) == bounds0
    return expanded, coarse, metadata


@contextmanager
def resolution_context(**changes):
    previous_cfg = CFG.copy(); previous_out = globals()['OUT']
    try:
        if 'output' in changes:
            globals()['OUT'] = Path(changes.pop('output'))
        CFG.update(changes)
        yield
    finally:
        CFG.clear(); CFG.update(previous_cfg); globals()['OUT'] = previous_out


def resolution_model():
    global RESOLUTION_DETECTOR
    if globals().get('RESOLUTION_DETECTOR') is not None:
        return RESOLUTION_DETECTOR
    location = RESOLUTION_CACHE / 'model'; location.mkdir(exist_ok=True)
    stamp = location / 'resolved_model.json'
    if stamp.exists():
        info = json.loads(stamp.read_text())
    else:
        # Resolve once; partial reruns retain the same immutable HF commit.
        from huggingface_hub import HfApi
        model_info = HfApi().model_info(CFG['model_name'], revision=CFG['model_revision'])
        info = dict(model_name=CFG['model_name'], source_revision=CFG['model_revision'],
                    resolved_revision=model_info.sha, deepforest_version=importlib.metadata.version('deepforest'),
                    limitation='3차 결과는 main으로만 기록되어 당시 가중치 SHA 동일 여부는 직접 검증할 수 없음')
        resolution_json(stamp, info)
    with resolution_context(output=location, model_revision=info['resolved_revision']):
        RESOLUTION_DETECTOR = load_detector()
    resolution_json(OUT / 'resolution_model_provenance.json', info)
    return RESOLUTION_DETECTOR


def resolution_candidates(plot, cm, rgb, chm):
    directory = RESOLUTION_CACHE / 'detections' / f'{cm}cm' / plot.plotID
    directory.mkdir(parents=True, exist_ok=True)
    processed = directory / 'candidates.csv'; stamp = directory / 'candidates_input.json'
    fingerprint = dict(schema=RESOLUTION_SCHEMA, code_sha=RESOLUTION_CODE_SHA,
        rgb_sha256=resolution_sha(rgb), chm_sha256=resolution_sha(chm), plot_wkt=plot.geometry.wkt,
        config=CFG, deepforest_version=importlib.metadata.version('deepforest'),
        resampling=json.loads((Path(rgb).parent / 'resampling.json').read_text()))
    # Model revision is stored once per experiment and in every prediction stamp.
    model_stamp = RESOLUTION_CACHE / 'model/resolved_model.json'
    fingerprint['resolved_model'] = json.loads(model_stamp.read_text()) if model_stamp.exists() else None
    if processed.exists() and stamp.exists() and json.loads(stamp.read_text()) == fingerprint:
        return resolution_read_geo(processed)
    detector = resolution_model()
    fingerprint['resolved_model'] = json.loads(model_stamp.read_text())
    with resolution_context(model_revision=fingerprint['resolved_model']['resolved_revision']):
        frame = deepforest_patch_candidates(detector, rgb, chm, plot, 600, directory / 'raw')
    frame['det_id'] = frame.det_id.map(lambda s: f'r{cm}_' + s)
    write_geo_csv(frame, processed.with_suffix('.part.csv'))
    processed.with_suffix('.part.csv').replace(processed)
    resolution_json(stamp, fingerprint)
    return frame


def resolution_verify_cached_winners(outer_records, inner_records):
    """Pure Python proof that the cached full-grid winner wins the subset."""
    winners = {}; totals = {}
    for row in inner_records:
        if int(row['patch_size']) != 600:
            continue
        key = (int(row['outer_fold']), float(row['score_threshold']))
        counts = totals.setdefault(key, [0, 0, 0])
        for i, column in enumerate(('TP', 'FP', 'FN')):
            counts[i] += int(row[column])
    assert {t for fold, t in totals} == set(RESOLUTION_THRESHOLDS)
    for row in outer_records:
        fold = int(row['outer_fold']); threshold = float(row['score_threshold'])
        assert int(row['patch_size']) == 600 and math.isclose(threshold, .2)
        def rank(t):
            tp, fp, fn = totals[fold, t]; den = 2 * tp + fp + fn
            return (-(2 * tp / den if den else -1), t)
        best = min(RESOLUTION_THRESHOLDS, key=rank)
        assert math.isclose(best, threshold), f'Cached subset winner differs in fold{fold}'
        winners[fold] = best
    assert set(winners) == set(range(5))
    return winners


def resolution_cached_10cm():
    outer = pd.read_csv(RESOLUTION_BASE / 'A_outer_groupkfold.csv')
    inner = pd.read_csv(RESOLUTION_BASE / 'A_inner_grid_scores.csv')
    splits = pd.read_csv(RESOLUTION_BASE / 'A_fold_plot_assignments.csv')
    deployment = pd.read_csv(RESOLUTION_BASE / 'A_deployment_grid.csv')
    # The original unrestricted optimum is also optimal in the patch600 subset.
    # Verify this explicitly rather than relabel a mixed-patch OOF as fixed600.
    restricted = inner[inner.patch_size == 600].copy()
    resolution_verify_cached_winners(outer.to_dict('records'), inner.to_dict('records'))
    result = RESOLUTION_BASE_METHODS['A_DF_nested_OOF']
    return result, outer, restricted, splits, deployment[deployment.patch_size == 600].copy()


def resolution_runtime_checks():
    # Real Rasterio average/bilinear test on a small synthetic checkerboard.
    with tempfile.TemporaryDirectory(prefix='resolution_check_') as td:
        td = Path(td); source = td / 'synthetic.tif'
        yy, xx = np.indices((100, 100)); tile = ((xx + yy) % 2 * 200 + 20).astype('uint8')
        with rasterio.open(source, 'w', driver='GTiff', width=100, height=100, count=3,
                           dtype='uint8', crs='EPSG:32617', transform=rasterio.Affine(.1, 0, 100, 0, -.1, 200)) as ds:
            ds.write(np.stack([tile] * 3))
        for cm in (25, 51):
            expanded, coarse, meta = resolution_rgb(source, cm, td / str(cm))
            with rasterio.open(coarse) as ds:
                assert np.var(ds.read(1)) < np.var(tile)
                assert abs(float(ds.read(1).mean()) - float(tile.mean())) < 1
            with rasterio.open(expanded) as ds:
                assert ds.width == ds.height == 100 and np.allclose(ds.res, .1)
            # The second request must use the finished cache without a rewrite.
            previous = expanded.stat().st_mtime_ns
            resolution_rgb(source, cm, td / str(cm))
            assert expanded.stat().st_mtime_ns == previous
    detection_algorithm_smoke_checks()
    original = RESOLUTION_BASE_METHODS['A_DF_nested_OOF']['plots']
    for _, p in original.iterrows():
        fa = field_all[field_all.plotID == p.plotID]
        assert len(fa) == p.n_truth
        if fa.ref_AGB_kg.notna().all():
            expected = fa.ref_AGB_kg.sum() / 1000 * 10000 / p.area_m2 * .5
            assert np.isclose(expected, p.C_ref_tC_ha)
    print('해상도/격자/실제 GSD/평균 보존/캐시 재사용/매칭·CV 격리/좌표 없는 참조 포함 검증 통과')


def resolution_save_result(result, folder):
    folder.mkdir(parents=True, exist_ok=True)
    for key, filename in [('detections', 'detected_trees.csv'), ('matched', 'matched_trees_visible.csv'),
                          ('metrics', 'detection_metrics.csv'), ('plots', 'plot_comparison.csv')]:
        write_geo_csv(result[key], folder / filename)
    for pop, frame in result['pairs'].items():
        frame.to_csv(folder / f'matches_{pop}.csv', index=False)


def resolution_score_text(outer, plot_ids):
    ids = set(plot_ids)
    selected = {int(r.outer_fold): float(r.score_threshold) for _, r in outer.iterrows()
                if ids.intersection(str(r.test_plots).split(';'))}
    values = set(selected.values())
    return next(iter(values)) if len(values) == 1 else json.dumps(selected, sort_keys=True)


def resolution_summaries(results, selections):
    rows = []; plotrows = []
    for (cm, setting), result in results.items():
        if cm == 10:
            # Reuse saved 10cm numbers, including the individual-error population.
            table = RESOLUTION_BASE_SUMMARY[RESOLUTION_BASE_SUMMARY.method == 'A_DF_nested_OOF'].copy()
        else:
            table = method_comparison_rows({f'{cm}_{setting}': result}, CFG['sites'])
        outer = selections[cm]['outer']
        for _, row in table.iterrows():
            r = row.to_dict(); r.pop('method', None)
            ids = plots.plotID if row.siteID == 'ALL' else plots.loc[plots.siteID == row.siteID, 'plotID']
            r.update(resolution_cm=cm, setting=setting, patch_size=600,
                     selected_score_threshold=.2 if setting == 'fixed' else resolution_score_text(outer, ids),
                     deployment_score_threshold=float(selections[cm]['deployment'].iloc[0].score_threshold),
                     source='cached_3rd_run_A' if cm == 10 else 'new_RGB_resolution_experiment')
            for old, new in [('height_Bias_m', 'height_Bias'), ('crown_Bias_m', 'crown_Bias'), ('DBH_Bias_cm', 'DBH_Bias'),
                             ('individual_AGB_Bias_kg', 'individual_AGB_Bias'), ('C_est_mean_tC_ha', 'C_est_mean'),
                             ('C_ref_mean_tC_ha', 'C_ref_mean'), ('carbon_RMSE_tC_ha', 'carbon_RMSE'), ('carbon_MAE_tC_ha', 'carbon_MAE')]:
                r[new] = r[old]
            rows.append(r)
        frame = result['plots'].copy(); frame['resolution_cm'] = cm; frame['setting'] = setting
        frame['selected_score_threshold'] = .2
        frame['outer_fold'] = -1
        if setting == 'retuned':
            for _, selection in outer.iterrows():
                mask = frame.plotID.isin(str(selection.test_plots).split(';'))
                frame.loc[mask, 'selected_score_threshold'] = selection.score_threshold
                frame.loc[mask, 'outer_fold'] = selection.outer_fold
        plotrows.append(frame)
    summary = pd.DataFrame(rows).sort_values(['resolution_cm', 'setting', 'siteID']).reset_index(drop=True)
    per_plot = pd.concat(plotrows, ignore_index=True).sort_values(['plotID', 'resolution_cm', 'setting'])
    assert len(summary) == 3 * 2 * 3 and len(per_plot) == 23 * 3 * 2
    assert per_plot.groupby('plotID').C_ref_tC_ha.nunique().eq(1).all()
    summary.to_csv(OUT / 'resolution_comparison.csv', index=False)
    per_plot.to_csv(OUT / 'resolution_plot_comparison.csv', index=False)
    return summary, per_plot


def resolution_validate_10cm(summary):
    base = RESOLUTION_BASE_SUMMARY.query("method == 'A_DF_nested_OOF' and siteID == 'ALL'").iloc[0]
    current = summary.query("resolution_cm == 10 and setting == 'fixed' and siteID == 'ALL'").iloc[0]
    rows = []
    for key, target, tolerance in [('F1_visible', .55, .02), ('estimate_reference_ratio', .68, .02)]:
        delta = float(current[key] - base[key]); passed = abs(delta) <= 1e-8
        approximate = abs(float(current[key]) - target) <= tolerance
        rows.append(dict(metric=key, cached_A=float(base[key]), fixed_10cm=float(current[key]),
                         difference=delta, matches_cached_A=passed, approximate_expected=target,
                         within_approximate_tolerance=approximate))
        if not passed or not approximate:
            warnings.warn(f'10cm 고정 결과와 기존 A/예상값 차이: {rows[-1]}')
    pd.DataFrame(rows).to_csv(OUT / 'resolution_10cm_validation.csv', index=False)
    print('10cm 고정 vs 기존 A 자동 비교 (각 outer fold600/.2 동일 확인):')
    display(pd.DataFrame(rows))
    return rows


def resolution_figures(summary, results, images):
    figures = OUT / 'figures'; figures.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for setting, label, color in [('fixed', 'Fixed: score .2 / patch600', '#3877a8'),
                                  ('retuned', 'Retuned: nested OOF / patch600', '#df8639')]:
        q = summary[(summary.siteID == 'ALL') & (summary.setting == setting)].sort_values('resolution_cm')
        for ax, y in zip(axes, ('F1_visible', 'estimate_reference_ratio')):
            ax.plot(q.resolution_cm, q[y], 'o-' if setting == 'fixed' else 's--', label=label, color=color)
            ax.set_xticks([10, 25, 51]); ax.set_xlabel('Nominal RGB GSD (cm)'); ax.grid(alpha=.2)
    axes[0].set_ylabel('Visible-tree F1'); axes[0].set_ylim(0, 1)
    axes[1].set_ylabel('Mean estimated / mean reference carbon'); axes[1].axhline(1, ls=':', c='gray')
    axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
    fig.suptitle('Same 23 plots; CHM held at 1m; RGB expanded to original 10cm grid')
    fig.tight_layout(); fig.savefig(figures / 'resolution_trend.png', dpi=170, bbox_inches='tight'); plt.close(fig)
    # Three columns; rows distinguish fixed settings from honest outer-OOF retuning.
    for pid in RESOLUTION_OVERLAYS:
        fig, axes = plt.subplots(2, 3, figsize=(17, 11), sharex=True, sharey=True)
        p = plots.set_index('plotID').loc[pid]
        for column, cm in enumerate(RESOLUTION_GSDS):
            rgb = images[cm][pid]
            for row, setting in enumerate(('fixed', 'retuned')):
                ax = axes[row, column]
                with rasterio.open(rgb) as ds:
                    show(ds, ax=ax); crs = ds.crs
                f = field_geo[field_geo.plotID == pid].to_crs(crs)
                d = results[cm, setting]['detections']; d = d[d.plotID == pid].to_crs(crs)
                if len(d): d.boundary.plot(ax=ax, color='cyan', linewidth=.8)
                for visible, color in [(True, 'lime'), (False, 'magenta')]:
                    z = f[f.visible == visible]
                    ax.scatter(z.geometry.x, z.geometry.y, c=color, s=19, edgecolors='black', linewidths=.35, zorder=4)
                footprint = gpd.GeoSeries([p.geometry], crs=WORK_CRS).to_crs(crs)
                footprint.boundary.plot(ax=ax, color='yellow', linewidth=1.5)
                bounds = footprint.total_bounds; padding = 8
                ax.set_xlim(bounds[0] - padding, bounds[2] + padding); ax.set_ylim(bounds[1] - padding, bounds[3] + padding)
                score = .2 if setting == 'fixed' else float(d.selected_score_threshold.iloc[0]) if len(d) and 'selected_score_threshold' in d else None
                ax.set_title(f'{cm}cm | {setting} | n={len(d)}' + (f' | score={score:g}' if score is not None else ''))
                ax.ticklabel_format(useOffset=False, style='plain'); ax.tick_params(labelsize=7)
        fig.suptitle(pid + ' | RGB resolution only; same CHM and footprints')
        fig.text(.5, .012, 'Cyan: detections | Green: visible field stems | Pink: occluded field stems | Yellow: plot boundary | Rows: fixed / nested OOF', ha='center', fontsize=10)
        fig.tight_layout(rect=(0, .035, 1, .96))
        fig.savefig(figures / f'{pid}_resolution_side_by_side.png', dpi=160, bbox_inches='tight'); plt.close(fig)


def resolution_report(summary, validation):
    columns = ['resolution_cm', 'setting', 'siteID', 'n_detected', 'F1_visible', 'F1_all',
               'C_est_mean', 'C_ref_mean', 'estimate_reference_ratio', 'carbon_RMSE', 'selected_score_threshold']
    table = summary[columns].to_markdown(index=False, floatfmt='.3f')
    all_sites = summary[summary.siteID == 'ALL'].set_index(['resolution_cm', 'setting'])
    changes = []
    for setting in ('fixed', 'retuned'):
        base = all_sites.loc[10, setting]
        for cm in (25, 51):
            q = all_sites.loc[cm, setting]
            changes.append(f'- {setting} {cm}cm: 10cm 대비 F1_visible {q.F1_visible-base.F1_visible:+.3f}, '
                           f'탄소 추정/참조 {q.estimate_reference_ratio-base.estimate_reference_ratio:+.3f}, '
                           f'탄소 RMSE {q.carbon_RMSE-base.carbon_RMSE:+.3f}tC/ha.')
    intro = '\n'.join(changes)
    model_file = RESOLUTION_CACHE / 'model/resolved_model.json'
    if model_file.exists(): shutil.copy2(model_file, OUT / 'resolution_model_provenance.json')
    section = f'''<!-- resolution-experiment:start -->
# 해상도 실험

원본10cm RGB를 면적 평균(Resampling.average)으로25/51cm로 줄이고, DeepForest 입력은
bilinear로 원본10cm 격자에 확대했다. **확대해도 원래 선명도와 세부 정보는 늘지 않는다.**
CHM은 동일한 원본1m 파일을 사용했다. H=상자 내부 CHM최대, CD=상자 가로·세로 평균,
CHM≥3m, Jucker biome11, Chojnacky 일반활엽수 대체식, CF0.5, point-in-box 헝가리안 매칭을 유지했다.

## 결과 요약

탄소 평균·RMSE 단위는 tC/ha, 높이·수관 Bias는m, DBH Bias는cm, 개체AGB Bias는kg이다.
탄소 비율은 유효 조사구의 평균 추정/평균 참조다. 상세 Bias·precision/recall은 resolution_comparison.csv 참고.

{table}

{intro}

## 실행·비교 검증

- GRSM15+ORNL8, 총23개 조사구 고정. 조사구 참조 총량에는 좌표 없는 나무를 포함했다.
  F1_all의 모집단은 기존과 동일하게 매칭 좌표가 있는 전체 나무이며, 좌표 없는 나무는 F1 분모로 만들지 않았다.
- 10cm baseline_015/A/B/C는 기존 결과를 재사용했다. 기존 A의5개 outer fold가 모두600/.2이고,
  patch600 안쪽 탐색에서도 같은 score가 최적임을 확인해10cm 고정/재선택 행에 그대로 재사용했다.
- 10cm 고정 결과 vs 기존A: {json.dumps(validation, ensure_ascii=False)}.
- 25/51cm 고정은 score.2, 재선택은 patch600에서.1/.2/.3/.4/.5를 안쪽 보이는 나무F1로 선택했다.
  바깥 평가 조사구의 정답은 선택에 쓰지 않았다. fold별 임계값과 train/test목록을 별도CSV로 저장했다.
  deployment_score_threshold는 전 조사구 선택값으로, OOF 평가에 역으로 적용하지 않았다.
- 실제 파일 CRS·bounds·GSD와 원본/확대본transform·shape의 일치를 검사했다.
  coarse 격자는 원본 범위와 좌상단을 보존하고 정수 크기로 만들므로 실제 GSD는25/51cm와 소폭 다를 수 있다.
  nominal과 실제x/y GSD는 resolution_raster_metadata.csv에 모두 기록했다.
- Drive crop만 읽고 새 RGB·CHM 다운로드는 수행하지 않았다. 줄인RGB/10cm 확대본/탐지 후보는Drive캐시에 보존했다.
- QUDO 충돌 수종계수 제외, gzip 응답 크기검사 예외, ORNL2021 provisional, HF_HOME=/content/hf_cache를 유지했다.

## 해석상 한계

- 인위적 평균 축소 영상은 실제25/51cm 촬영의 센서·광학특성·압축·촬영각·노출·계절과 다를 수 있다.
- CHM을 고정했으므로 RGB 해상도 변화만 분리한 실험이며, 국내 CHM의 오차/해상도 영향은 평가하지 않는다.
  상자가 달라지면 같은 CHM에서도 추출 최대 수고와 CD가 달라진다.
- 이23개 조사구에서 배포 설정이 이미 선택되었으므로 고정 설정의10cm 점수는 독립적인 외부 검증이 아니다.
  NEON 사전학습 모델·지역·수종·산림구조가 국내 영상과 달라 국내 성능을 보장하지 않는다.
- 기존3차 모델은 revision=main으로만 기록되었다. 새 추론은 동일DeepForest버전과 실험 첫 실행 때 고정한
  HF커밋을 사용하나, 과거 가중치 SHA가 없어 당시와 동일한 가중치인지 직접 검증할 수 없다.
  resolution_model_provenance.json을 보존한다. 재실행 시 새해상도들 사이의 커밋은 바꾸지 않는다.
- 실측 기반NEONForestAGB 역시 상대생장식 참조 추정값이다. 파괴조사 무게 정답이나 수관IoU 검증이 아니다.

[Rasterio 리샘플링 문서](https://rasterio.readthedocs.io/en/stable/topics/resampling.html)

그림: figures/resolution_trend.png 및 세 *_resolution_side_by_side.png.
기존 결과 출처: {RESOLUTION_BASE}; 해상도 실험ID: {RUN_ID}.
<!-- resolution-experiment:end -->

'''
    previous = (RESOLUTION_BASE / 'report.md').read_text(encoding='utf-8')
    if '<!-- resolution-experiment:end -->' in previous:
        previous = previous.split('<!-- resolution-experiment:end -->', 1)[1].lstrip()
    (OUT / 'report.md').write_text(section + previous, encoding='utf-8')
    pointer = ROOT / 'outputs' / ('resolution_latest.md' if RESOLUTION_BASE.resolve() == (ROOT/'outputs').resolve() else 'report.md')
    pointer.write_text(f'최근 완료 결과: [{RUN_ID}/report.md]({RUN_ID}/report.md)\n', encoding='utf-8')


def resolution_completed_outputs():
    outputs = ['resolution_comparison.csv', 'resolution_plot_comparison.csv',
            'resolution_raster_metadata.csv', 'resolution_10cm_validation.csv', 'report.md',
            'figures/resolution_trend.png'] + [f'figures/{pid}_resolution_side_by_side.png' for pid in RESOLUTION_OVERLAYS]
    for cm in RESOLUTION_GSDS:
        outputs += [f'resolution/{cm}cm/{name}' for name in ('outer_groupkfold.csv', 'inner_grid_scores.csv',
                    'fold_plot_assignments.csv', 'deployment_grid.csv')]
        for setting in ('fixed', 'retuned'):
            outputs += [f'resolution/{cm}cm/{setting}/{name}' for name in ('detected_trees.csv',
                        'matched_trees_visible.csv', 'detection_metrics.csv', 'plot_comparison.csv',
                        'matches_visible.csv', 'matches_all.csv')]
    return outputs


def run_resolution_experiment():
    global RESOLUTION_SUMMARY, RESOLUTION_PLOT_COMPARISON, RESOLUTION_RESULTS
    resolution_load_state()
    completed = OUT / 'resolution_completed.json'
    if completed.is_file():
        record = json.loads(completed.read_text())
        if record['base_signature'] == RESOLUTION_BASE_SIGNATURE and all(
            (OUT / name).is_file() and resolution_sha(OUT / name) == digest
            for name, digest in record['outputs'].items()
        ):
            RESOLUTION_SUMMARY = pd.read_csv(OUT / 'resolution_comparison.csv')
            RESOLUTION_PLOT_COMPARISON = pd.read_csv(OUT / 'resolution_plot_comparison.csv')
            print('완료된 동일 실험 재사용: RGB 재생성·GPU 추론·CV/평가 재계산 없이 저장 결과를 불러왔습니다.')
            display(pd.read_csv(OUT / 'resolution_10cm_validation.csv'))
            return RESOLUTION_SUMMARY
    resolution_runtime_checks()
    results = {}; selections = {}; metadata = []
    images = {10: {pid: paths[0] for pid, paths in RASTERS.items()}}
    cached, outer, inner, splits, deployment = resolution_cached_10cm()
    for setting in ('fixed', 'retuned'): results[10, setting] = cached
    selections[10] = dict(outer=outer, inner=inner, splits=splits, deployment=deployment)
    for cm in (25, 51):
        frames = []; images[cm] = {}
        for _, p in plots.iterrows():
            rgb, chm = RASTERS[p.plotID]
            model_input, coarse, meta = resolution_rgb(rgb, cm, RESOLUTION_CACHE / 'rgb' / f'{cm}cm' / p.plotID)
            images[cm][p.plotID] = model_input
            meta.update(plotID=p.plotID, resolution_cm=cm, coarse_path=str(coarse), model_input_path=str(model_input),
                        chm_path=str(chm), chm_sha256=resolution_sha(chm), chm_gsd_m=1.0)
            metadata.append(meta)
            print(f'{p.plotID} {cm}cm: actual GSD={meta["actual_gsd_x_m"]:.6f} × '
                  f'{meta["actual_gsd_y_m"]:.6f}m; input=10cm original grid; CHM=1m unchanged')
            frames.append(resolution_candidates(p, cm, model_input, chm))
        candidates = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), geometry='geometry', crs=WORK_CRS)
        grid = {(600, t): candidates[candidates.score >= t].copy() for t in RESOLUTION_THRESHOLDS}
        oof, outer, inner, splits, deployment = nested_group_parameter_selection(
            grid, plots, field_geo[field_geo.visible], CFG['cv_folds'])
        reference_splits = selections[10]['splits'].sort_values(['outer_fold', 'plotID', 'role']).reset_index(drop=True)
        pd.testing.assert_frame_equal(splits.sort_values(['outer_fold', 'plotID', 'role']).reset_index(drop=True), reference_splits)
        selections[cm] = dict(outer=outer, inner=inner, splits=splits, deployment=deployment)
        results[cm, 'fixed'] = evaluate_detection_method(grid[600, .2], plots, field_all, field_geo)
        results[cm, 'retuned'] = evaluate_detection_method(oof, plots, field_all, field_geo)
    for pid, rgb in images[10].items():
        with rasterio.open(rgb) as ds:
            metadata.append(dict(plotID=pid, resolution_cm=10, nominal_cm=10,
                actual_gsd_x_m=ds.res[0], actual_gsd_y_m=ds.res[1], bounds=list(ds.bounds),
                input_transform=list(ds.transform), input_width=ds.width, input_height=ds.height,
                coarse_path=str(rgb), model_input_path=str(rgb), crs=ds.crs.to_string(),
                chm_path=str(RASTERS[pid][1]), chm_sha256=resolution_sha(RASTERS[pid][1]), chm_gsd_m=1.0,
                downsampling='none: original 10cm', upsampling='none'))
    pd.DataFrame(metadata).to_csv(OUT / 'resolution_raster_metadata.csv', index=False)
    for cm, selection in selections.items():
        folder = OUT / 'resolution' / f'{cm}cm'; folder.mkdir(parents=True, exist_ok=True)
        for key, filename in [('outer', 'outer_groupkfold.csv'), ('inner', 'inner_grid_scores.csv'),
                              ('splits', 'fold_plot_assignments.csv'), ('deployment', 'deployment_grid.csv')]:
            selection[key].to_csv(folder / filename, index=False)
        for setting in ('fixed', 'retuned'):
            if cm == 10:
                # Byte-for-byte copies; no old inference, evaluation or carbon recomputation.
                shutil.copytree(RESOLUTION_BASE / 'methods/A_DF_nested_OOF', folder / setting, dirs_exist_ok=True)
            else:
                resolution_save_result(results[cm, setting], folder / setting)
    summary, per_plot = resolution_summaries(results, selections)
    validation = resolution_validate_10cm(summary)
    resolution_figures(summary, results, images)
    resolution_report(summary, validation)
    # Confirm old CSVs have not been mutated by evaluation/plotting or copying.
    provenance = json.loads((OUT / 'resolution_provenance.json').read_text())
    for name, digest in provenance['source_files'].items():
        assert resolution_sha(RESOLUTION_BASE / name) == digest, f'기존 결과 변경 감지: {name}'
    RESOLUTION_SUMMARY = summary; RESOLUTION_PLOT_COMPARISON = per_plot; RESOLUTION_RESULTS = results
    resolution_json(completed, dict(base_signature=RESOLUTION_BASE_SIGNATURE,
                    outputs={name: resolution_sha(OUT / name) for name in resolution_completed_outputs()}))
    return summary
