"""Find GEDI L4A V3 footprints over Geumjeongsan and the PNU Busan campus.

Run with the gedi_hongcheon venv:
  ..\\gedi_hongcheon\\.venv\\Scripts\\python.exe scripts\\find_gedi_sites.py
Needs environment variables EARTHDATA_USERNAME / EARTHDATA_PASSWORD.

Granules are streamed (only the needed parts of each HDF5 file are read), so
matching granules are not downloaded in full. Failed streams alone fall back
to a local download. --limit 2 --collect-only benchmarks two new granules;
running again without a limit resumes from the same atomic per-granule cache.

Uses only numpy/h5py/scipy/matplotlib: Windows application control on this
PC currently blocks the pandas and rasterio DLLs. Tree cover comes from the
GEDI file itself (land_cover_data/landsat_treecover, pft_class).
"""
from pathlib import Path
import argparse
import csv
import json
import math
import os
import sys
import time
import urllib.parse
import urllib.request
from xml.sax.saxutils import escape

import numpy as np

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs'
RAW = ROOT / 'data'
SHORT_NAME = 'GEDI_L4A_AGB_Density_V3_2508'
BBOX = (129.00, 35.20, 129.12, 35.34)
CACHE = RAW / 'gedi_v3_geumjeongsan'
TEMPORAL = ('2019-04-01T00:00:00Z', '2025-12-31T23:59:59Z')
BEAMS = ('BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011')
POWER = {'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011'}
VARS = ('shot_number', 'lat_lowestmode', 'lon_lowestmode', 'agbd', 'agbd_se', 'l4a_quality_flag_rel3',
        'degrade_flag', 'sensitivity', 'solar_elevation', 'delta_time')
LC_VARS = ('landsat_treecover', 'pft_class')
PFT = {0: '수역', 1: '상록침엽수', 2: '상록활엽수', 3: '낙엽침엽수', 4: '낙엽활엽수',
       5: '관목', 6: '초지', 7: '곡물 경작지', 8: '활엽 작물 경작지', 9: '시가지',
       10: '만년설·빙하', 11: '나지·희박식생', 255: '결측'}
COLUMNS = ['shot_number', 'lat', 'lon', 'agbd', 'agbd_se', 'l4a_quality_flag_rel3', 'degrade_flag', 'sensitivity',
           'solar_elevation', 'delta_time', 'landsat_treecover', 'pft_class', 'beam', 'granule']

# Ordered, exclusive zones: campus first, then the remainder of BBOX.
ZONES = {
    '부산대 부산캠퍼스': (129.0738, 35.2299, 129.0847, 35.2386),
    '금정산 일대': BBOX,
}
TRAIL_INFO = {}
TRAIL_LINES = []
LANDMARKS = {'부산대': (35.2343, 129.0784), '고당봉': (35.2801, 129.0506), '범어사': (35.2841, 129.0687)}


def zone_of(lat, lon):
    for name, (w, s, e, n) in ZONES.items():
        if w <= lon <= e and s <= lat <= n:
            return name
    return '그 외'


def credentials_ok():
    if os.getenv('EARTHDATA_USERNAME') and os.getenv('EARTHDATA_PASSWORD'):
        return 'environment'
    return None


def read_granule(fobj, name, stats=None):
    """Rows (dicts) for every shot inside BBOX in one granule."""
    import h5py
    w, s, e, n = BBOX
    rows = []
    def read(dataset, selection=slice(None)):
        values = dataset[selection]
        if stats is not None:
            stats['dataset_bytes'] += values.nbytes
        return values
    with h5py.File(fobj, 'r') as f:
        for beam in BEAMS:
            if beam not in f:
                continue
            g = f[beam]
            lat, lon = read(g['lat_lowestmode']), read(g['lon_lowestmode'])
            idx = np.flatnonzero((lat >= s) & (lat <= n) & (lon >= w) & (lon <= e))
            if not len(idx):
                continue
            lo, hi = idx.min(), idx.max() + 1  # one contiguous slice = few remote reads
            # Reuse coordinates already read instead of fetching them again.
            cols = {v: read(g[v], slice(lo, hi))[idx - lo]
                    for v in VARS if v not in ('lat_lowestmode', 'lon_lowestmode')}
            cols['lat_lowestmode'], cols['lon_lowestmode'] = lat[idx], lon[idx]
            for v in LC_VARS:
                if 'land_cover_data' in g and v in g['land_cover_data']:
                    cols[v] = read(g['land_cover_data'][v], slice(lo, hi))[idx - lo]
            cols['lat'], cols['lon'] = cols.pop('lat_lowestmode'), cols.pop('lon_lowestmode')
            for i in range(len(idx)):
                row = {k: (int(v[i]) if k == 'shot_number' else v[i].item()) for k, v in cols.items()}
                row['beam'], row['granule'] = beam, name
                rows.append(row)
    return rows


def load_cache(path):
    if not path.exists():
        return []
    with path.open(encoding='utf-8-sig', newline='') as fh:
        out = []
        for r in csv.DictReader(fh):
            for k in r:
                if k == 'shot_number':
                    r[k] = int(r[k])
                elif k not in ('beam', 'granule') and r[k] != '':
                    r[k] = float(r[k])
            out.append(r)
        return out


def save_rows(rows, path, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(path.suffix + '.tmp')
    with pending.open('w', encoding='utf-8-sig', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=columns, extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)
    pending.replace(path)


def save_json(data, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(path.suffix + '.tmp')
    pending.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
    pending.replace(path)


def worker_login():
    import earthaccess
    if not earthaccess.login(strategy='environment', persist=False).authenticated:
        raise RuntimeError('Earthdata worker authentication failed')


def fetch_granule(g):
    """Read one independent granule; subprocesses never share an HDF5 handle."""
    import earthaccess
    name = g['meta']['native-id']
    t0 = time.monotonic()
    stats = {'granule': name, 'dataset_bytes': 0, 'remote_bytes': 0, 'range_requests': 0,
             'download_bytes': 0, 'stream_error_type': ''}
    remote = None
    try:
        remote = earthaccess.open([g], show_progress=False,
                                  open_kwargs={'block_size': 2**20, 'cache_type': 'blockcache'})[0]
        stats['file_bytes'] = remote.size
        # Response body bytes, including block overhead but excluding HTTP/TLS.
        fetch = remote.cache.fetcher
        def counted_fetch(a, b):
            body = fetch(a, b)
            stats['remote_bytes'] += len(body)
            stats['range_requests'] += 1
            return body
        remote.cache.fetcher = counted_fetch
        new = read_granule(remote, name, stats)
        how = 'stream'
    except Exception as exc:
        stats['stream_error_type'] = type(exc).__name__
        print(f'  stream failed ({type(exc).__name__}); downloading this granule only: {name}', flush=True)
        fallback = CACHE / 'fallback_downloads'
        fallback.mkdir(exist_ok=True)
        paths = earthaccess.download([g], local_path=str(fallback), threads=1, show_progress=False)
        if not paths:
            raise RuntimeError(f'Fallback download failed: {name}') from exc
        path = Path(paths[0])
        stats['download_bytes'] = path.stat().st_size
        stats['file_bytes'] = path.stat().st_size
        new = read_granule(path, name, stats)
        how = 'download'
    finally:
        if remote is not None:
            remote.close()
    stats.update(method=how, seconds=round(time.monotonic() - t0, 3), shots=len(new))
    return new, stats


def granule_results(todo, workers):
    if workers == 1:
        for granule in todo:
            yield fetch_granule(granule)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        # HDF5 file-object callbacks have process-local locks. Use separate
        # processes, never shared HDF5 readers or threads.
        with ProcessPoolExecutor(max_workers=workers, initializer=worker_login) as pool:
            pending = [pool.submit(fetch_granule, granule) for granule in todo]
            for future in as_completed(pending):
                yield future.result()


def collect(limit=None, workers=1):
    strategy = credentials_ok()
    if strategy is None:
        print('Earthdata 인증 정보가 없습니다. EARTHDATA_USERNAME / EARTHDATA_PASSWORD 를 설정한 뒤 다시 실행하세요.')
        sys.exit(2)
    import earthaccess
    CACHE.mkdir(parents=True, exist_ok=True)
    scope = {'bbox': list(BBOX), 'temporal': list(TEMPORAL), 'short_name': SHORT_NAME}
    scope_path = CACHE / 'scope.json'
    if scope_path.exists() and json.loads(scope_path.read_text(encoding='utf-8')) != scope:
        raise RuntimeError('Cache scope differs; use a separate cache directory')
    save_json(scope, scope_path)
    # One atomic CSV per granule is its completion marker, including empty files.
    # A crash cannot mark incomplete rows as finished.
    metrics_path = CACHE / 'granule_metrics.json'
    metrics = json.loads(metrics_path.read_text(encoding='utf-8')) if metrics_path.exists() else []
    if not earthaccess.login(strategy=strategy, persist=False).authenticated:
        raise RuntimeError('Earthdata authentication failed')
    granules = sorted(earthaccess.search_data(short_name=SHORT_NAME, bounding_box=BBOX,
                                            temporal=TEMPORAL, count=-1),
                      key=lambda g: g['meta']['native-id'])
    names = [g['meta']['native-id'] for g in granules]
    if any(Path(name).name != name for name in names):
        raise RuntimeError('Unexpected granule name')
    save_json({'query': scope, 'count': len(granules), 'granules': names}, CACHE / 'manifest.json')
    done = {name for name in names if (CACHE / 'granules' / (name + '.csv')).exists()}
    todo = [g for g in granules if g['meta']['native-id'] not in done]
    print(f'V3 granules: {len(granules)} (남은 것 {len(todo)})', flush=True)
    if limit is not None:
        todo = todo[:limit]
    start = time.monotonic()
    for i, (new, stats) in enumerate(granule_results(todo, workers), 1):
        name, how = stats['granule'], stats['method']
        metrics.append(stats)
        save_json(metrics, metrics_path)
        save_rows(new, CACHE / 'granules' / (name + '.csv'), COLUMNS)
        done.add(name)
        save_json(sorted(done), CACHE / 'granules_done.json')
        el = time.monotonic() - start
        print(f'[{i}/{len(todo)}] {name}: {len(new)} shots, {how}, {stats["seconds"]:.1f}s, '
              f'원격 {stats["remote_bytes"] / 2**20:.2f} MiB, 다운로드 {stats["download_bytes"] / 2**20:.2f} MiB '
              f'(경과 {el / 60:.1f}분, 남은 예상 {el / i * (len(todo) - i) / 60:.0f}분)', flush=True)
    rows = [row for name in names if name in done
            for row in load_cache(CACHE / 'granules' / (name + '.csv'))]
    run_info = {'query': scope, 'total_granules': len(granules), 'completed_granules': len(done),
                'complete': len(done) == len(granules), 'metrics': metrics,
                'workers': workers, 'collection_wall_seconds': round(time.monotonic() - start, 3)}
    return rows, run_info


def osm_paths():
    """Hiking paths and tracks from OpenStreetMap, for 'how far to walk' per point."""
    cache = RAW / 'osm_paths.json'
    if not cache.exists():
        w, s, e, n = BBOX
        q = f'[out:json][timeout:90];way["highway"~"path|footway|track|steps"]({s},{w},{n},{e});out geom;'
        body = urllib.parse.urlencode({'data': q}).encode()
        for server in ('https://overpass-api.de/api/interpreter', 'https://overpass.kumi.systems/api/interpreter',
                       'https://overpass.private.coffee/api/interpreter'):
            try:
                req = urllib.request.Request(server, data=body, headers={'User-Agent': 'gluri-intern-research/0.1'})
                cache.write_bytes(urllib.request.urlopen(req, timeout=120).read())
                break
            except OSError as exc:
                print(f'  Overpass {server} 실패 ({exc})', flush=True)
        else:
            raise RuntimeError('All Overpass servers failed')
    data = json.loads(cache.read_text(encoding='utf-8'))
    from pyproj import Transformer
    from shapely.geometry import LineString
    transform = Transformer.from_crs(4326, 5179, always_xy=True)
    lines = []
    for element in data.get('elements', []):
        nodes = element.get('geometry', [])
        if len(nodes) >= 2:
            x, y = transform.transform([p['lon'] for p in nodes], [p['lat'] for p in nodes])
            lines.append(LineString(zip(x, y)))
    return lines


def forest_paths():
    """Read line geometry only: no DBF decoding or GDAL dependencies."""
    import shapefile
    from pyproj import CRS
    from shapely.geometry import LineString
    lines, sources = [], []
    for district in ('금정구', '북구', '양산시'):
        files = sorted((RAW / 'forest_map' / '등산로' / district).glob('*.shp'))
        if not files:
            raise FileNotFoundError(f'Missing trail shapefile: {district}')
        for path in files:
            prj = path.with_suffix('.prj')
            if prj.exists() and CRS.from_wkt(prj.read_text(encoding='utf-8-sig')).to_epsg() != 5179:
                raise ValueError(f'Trail CRS is not EPSG:5179: {path.name}')
            before = len(lines)
            with shapefile.Reader(shp=str(path), shx=str(path.with_suffix('.shx'))) as reader:
                for shape in reader.iterShapes():
                    if shape.shapeType == 0:
                        continue
                    if shape.shapeType not in (3, 13, 23):
                        raise ValueError(f'Expected polyline, got shape type {shape.shapeType}')
                    parts = list(shape.parts) + [len(shape.points)]
                    for a, b in zip(parts, parts[1:]):
                        coords = np.asarray(shape.points[a:b], dtype=float)[:, :2]
                        if len(coords) >= 2 and np.all(np.isfinite(coords)):
                            line = LineString(coords)
                            if line.length > 0:
                                lines.append(line)
            if len(lines) == before:
                raise ValueError(f'Empty trail file: {district}/{path.name}')
            sources.append({'path': str(path.relative_to(ROOT)), 'lines': len(lines) - before})
    return lines, sources


def load_paths():
    global TRAIL_INFO, TRAIL_LINES
    try:
        TRAIL_LINES, sources = forest_paths()
        TRAIL_INFO = {'source': '산림청 등산로', 'crs': 'EPSG:5179', 'files': sources,
                      'distance_method': 'point to nearest polyline segment'}
    except Exception as exc:
        print(f'산림청 등산로 읽기 실패 ({type(exc).__name__}); OpenStreetMap으로 대체', flush=True)
        TRAIL_LINES = osm_paths()
        TRAIL_INFO = {'source': 'OpenStreetMap', 'crs': 'EPSG:5179',
                      'fallback_reason': str(exc), 'distance_method': 'point to nearest polyline segment'}
    if not TRAIL_LINES:
        raise RuntimeError('No readable trail geometry')
    print(f'등산로: {TRAIL_INFO["source"]}, {len(TRAIL_LINES):,}개 선형', flush=True)
    return TRAIL_LINES


def path_distance(lat, lon, lines):
    from pyproj import Transformer
    import shapely
    if not lines or not len(lat):
        return np.full(len(lat), np.nan)
    x, y = Transformer.from_crs(4326, 5179, always_xy=True).transform(lon, lat)
    points = shapely.points(x, y)
    tree = shapely.STRtree(lines)
    return shapely.distance(points, tree.geometries[tree.nearest(points)])


def enrich(rows):
    seen, out = set(), []
    for r in rows:
        if zone_of(r['lat'], r['lon']) == '그 외':
            continue
        if r['shot_number'] in seen:
            continue
        seen.add(r['shot_number'])
        r['date'] = time.strftime('%Y-%m-%d', time.gmtime(1514764800 + float(r['delta_time'])))  # GEDI epoch 2018-01-01
        r['zone'] = zone_of(r['lat'], r['lon'])
        r['power'] = r['beam'] in POWER
        agbd = r['agbd']
        r['quality'] = (agbd != '' and math.isfinite(agbd) and agbd >= 0
                        and r['l4a_quality_flag_rel3'] == 1 and r['degrade_flag'] == 0)
        pft = r.get('pft_class')
        r['pft_name'] = PFT.get(int(pft), '기타') if pft != '' and pft is not None else ''
        tc = r.get('landsat_treecover')
        threshold = 20 if r['zone'] == '부산대 부산캠퍼스' else 30
        r['forest'] = tc not in ('', None) and math.isfinite(tc) and threshold <= tc <= 100
        out.append(r)
    try:
        d = path_distance(np.array([r['lat'] for r in out]), np.array([r['lon'] for r in out]), load_paths())
    except Exception as exc:
        print(f'등산로 거리 계산 실패 ({type(exc).__name__}); 금정산 후보는 선정하지 않음', flush=True)
        TRAIL_INFO.update(error=type(exc).__name__)
        d = np.full(len(out), np.nan)
    for r, x in zip(out, d):
        r['path_dist_m'] = float(x)  # Keep full precision for the 150 m threshold.
        r['path_source'] = TRAIL_INFO.get('source', 'unavailable')
    return out


def pick_sites(good):
    """15 mountain sites >=300 m apart, <=150 m from a trail; 10 campus sites.
    Campus has no trail-distance or minimum-spacing requirement.
    Preserve power-beam/night preference without an extra AGBD cutoff.
    """
    from pyproj import Geod
    geod = Geod(ellps='WGS84')
    cand = [r for r in good if r['quality'] and r['forest']
            and (r['zone'] == '부산대 부산캠퍼스'
                 or (math.isfinite(r['path_dist_m']) and r['path_dist_m'] <= 150))]
    for r in cand:
        r['score'] = 2 * r['power'] + (r['solar_elevation'] < 0)
        if r['zone'] == '금정산 일대':
            r['score'] += r['path_dist_m'] <= 50
    sites = []
    for zone in ZONES:
        chosen = []
        campus = zone == '부산대 부산캠퍼스'
        for r in sorted((r for r in cand if r['zone'] == zone),
                        key=lambda r: (-r['score'], r['agbd_se'] if campus else r['path_dist_m'], r['shot_number'])):
            if campus or all(geod.inv(r['lon'], r['lat'], c['lon'], c['lat'])[2] >= 300 for c in chosen):
                chosen.append(r)
            if len(chosen) >= (10 if campus else 15):
                break
        sites += chosen
    for i, r in enumerate(sites, 1):
        r['site_id'] = f'S{i:02d}'
        r['google_maps'] = f"https://maps.google.com/?q={r['lat']:.6f},{r['lon']:.6f}"
    return sites


def write_kml(sites, path):
    marks = []
    for r in sites:
        desc = (f"{r['zone']} | {r['pft_name']} | GEDI V3 AGBD {r['agbd']:.0f} Mg/ha (SE {r['agbd_se']:.0f}) | "
                f"{r['date']} | {r['beam']} | 등산로까지 약 {r['path_dist_m']:.0f}m | 약 25m 지름 GEDI footprint의 중심점")
        marks.append(f"<Placemark><name>{r['site_id']}</name><description>{escape(desc)}</description>"
                     f"<Point><coordinates>{r['lon']:.6f},{r['lat']:.6f},0</coordinates></Point></Placemark>")
    path.write_text('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
                    '<name>GEDI 현장 방문 후보</name>' + ''.join(marks) + '</Document></kml>', encoding='utf-8')


def write_gpx(sites, path):
    pts = [f"<wpt lat=\"{r['lat']:.6f}\" lon=\"{r['lon']:.6f}\"><name>{r['site_id']}</name>"
           f"<desc>{escape(r['zone'])} AGBD {r['agbd']:.0f}</desc></wpt>" for r in sites]
    path.write_text('<?xml version="1.0" encoding="UTF-8"?>\n<gpx version="1.1" creator="gluri" '
                    'xmlns="http://www.topografix.com/GPX/1/1">' + ''.join(pts) + '</gpx>', encoding='utf-8')


def figure(rows, good, sites):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    plt.rcParams['font.family'] = 'Malgun Gothic'
    plt.rcParams['axes.unicode_minus'] = False
    fig, ax = plt.subplots(figsize=(10, 10))
    if TRAIL_LINES:
        from matplotlib.collections import LineCollection
        from pyproj import Transformer
        transform = Transformer.from_crs(5179, 4326, always_xy=True)
        segments = []
        for line in TRAIL_LINES:
            coords = np.asarray(line.coords)
            lon, lat = transform.transform(coords[:, 0], coords[:, 1])
            segments.append(np.c_[lon, lat])
        ax.add_collection(LineCollection(segments, colors='#a87e57', linewidths=.45, alpha=.6,
                                        label=TRAIL_INFO.get('source', '등산로')))
    ax.scatter([r['lon'] for r in rows], [r['lat'] for r in rows], s=1, c='0.8', label=f'전체 GEDI 점 ({len(rows):,})')
    sc = ax.scatter([r['lon'] for r in good], [r['lat'] for r in good], s=4,
                    c=np.clip([r['agbd'] for r in good], 0, 400), cmap='YlGn', label=f'품질 통과 ({len(good):,})')
    if sites:
        ax.scatter([r['lon'] for r in sites], [r['lat'] for r in sites], s=70, marker='*', c='red',
                   label=f'현장 방문 후보 ({len(sites)})')
        for r in sites:
            if r['zone'] != '부산대 부산캠퍼스':
                dx, align = (-4, 'right') if r['lon'] > BBOX[2] - .004 else (4, 'left')
                ax.annotate(r['site_id'], (r['lon'], r['lat']), fontsize=7, xytext=(dx, 3),
                            ha=align, textcoords='offset points')
    for name, (w, s, e, n) in ZONES.items():
        ax.add_patch(Rectangle((w, s), e - w, n - s, fill=False, ls='--', lw=1, ec='tab:blue'))
        ax.text(w, n, name, fontsize=9, color='tab:blue', va='bottom')
    for name, (la, lo) in LANDMARKS.items():
        ax.plot(lo, la, 'k^', ms=6)
        ax.annotate(name, (lo, la), fontsize=8, xytext=(4, -10), textcoords='offset points')
    ax.set_aspect(1 / math.cos(math.radians(35.27)))
    ax.set(xlim=(BBOX[0], BBOX[2]), ylim=(BBOX[1], BBOX[3]), xlabel='경도', ylabel='위도',
           title='금정산·부산대 부산캠퍼스 GEDI L4A V3 (2019–2025)')
    fig.colorbar(sc, ax=ax, label='AGBD (Mg/ha, 품질 통과)', shrink=.7)
    ax.legend(loc='lower right', fontsize=8)
    fig.savefig(OUT / 'gedi_sites_map.png', dpi=170, bbox_inches='tight')
    plt.close(fig)
    # A separate campus map keeps the ten nearby labels readable.
    w, s, e, n = ZONES['부산대 부산캠퍼스']
    campus_rows = [r for r in rows if r['zone'] == '부산대 부산캠퍼스']
    campus_good = [r for r in good if r['zone'] == '부산대 부산캠퍼스']
    campus_sites = [r for r in sites if r['zone'] == '부산대 부산캠퍼스']
    fig, ax = plt.subplots(figsize=(9, 8))
    if TRAIL_LINES:
        ax.add_collection(LineCollection(segments, colors='#a87e57', linewidths=1, alpha=.7,
                                        label=TRAIL_INFO.get('source', '등산로')))
    ax.scatter([r['lon'] for r in campus_rows], [r['lat'] for r in campus_rows],
               s=13, c='0.78', label=f'전체 ({len(campus_rows)})')
    sc = ax.scatter([r['lon'] for r in campus_good], [r['lat'] for r in campus_good],
                    s=25, c=np.clip([r['agbd'] for r in campus_good], 0, 400),
                    cmap='YlGn', vmin=0, vmax=400, label=f'품질 통과 ({len(campus_good)})')
    ax.scatter([r['lon'] for r in campus_sites], [r['lat'] for r in campus_sites],
               s=110, marker='*', c='red', label=f'현장 후보 ({len(campus_sites)})', zorder=5)
    for i, r in enumerate(campus_sites):
        offset = (8, 7) if i % 2 == 0 else (8, -12)
        ax.annotate(r['site_id'], (r['lon'], r['lat']), fontsize=9, xytext=offset,
                    textcoords='offset points', zorder=6,
                    bbox=dict(boxstyle='round,pad=.12', fc='white', ec='none', alpha=.85))
    ax.add_patch(Rectangle((w, s), e-w, n-s, fill=False, ls='--', ec='tab:blue'))
    ax.set(xlim=(w-.0003, e+.0003), ylim=(s-.0003, n+.0003), xlabel='경도', ylabel='위도',
           title='부산대 부산캠퍼스 GEDI 관측점·현장 후보 확대도')
    ax.set_aspect(1 / math.cos(math.radians(35.234)))
    ax.ticklabel_format(axis='both', useOffset=False)
    ax.grid(alpha=.15)
    ax.legend(loc='upper right', fontsize=8)
    fig.colorbar(sc, ax=ax, label='AGBD (Mg/ha, 품질 통과)', shrink=.75)
    fig.savefig(OUT / 'gedi_campus_detail.png', dpi=180, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, help='Process at most N uncached granules')
    parser.add_argument('--collect-only', action='store_true', help='Benchmark/cache only; no final site outputs')
    parser.add_argument('--workers', type=int, default=1, help='Independent reader processes (1 to 8)')
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error('--limit must be positive')
    if not 1 <= args.workers <= 8:
        parser.error('--workers must be between 1 and 8')
    raw_rows, run_info = collect(args.limit, args.workers)
    report_dir = OUT if run_info['complete'] else OUT / 'trial'
    report_dir.mkdir(parents=True, exist_ok=True)
    save_json(run_info, report_dir / 'run_metadata.json')
    metrics = run_info['metrics']
    if metrics:
        save_rows(metrics, report_dir / 'granule_metrics.csv', list(metrics[0]))
        mean_s = float(np.mean([r['seconds'] for r in metrics]))
        print(f'평균 {mean_s:.1f}초/파일; 83개 예상 {83 * mean_s / 60:.1f}분; '
              f'현재 검색 {run_info["total_granules"]}개, 완료 {run_info["completed_granules"]}개', flush=True)
    if args.collect_only:
        return
    if not run_info['complete']:
        raise RuntimeError('Collection incomplete. Resume without --limit before generating final outputs.')
    export_results(raw_rows, run_info)


def export_results(raw_rows, run_info):
    """Regenerate outputs from an already complete cache without remote reads."""
    if not run_info['complete']:
        raise RuntimeError('Cannot export a partial collection')
    run_info['generated_at_kst'] = time.strftime('%Y-%m-%d %H:%M:%S KST', time.gmtime(time.time() + 9 * 3600))
    rows = enrich(raw_rows)
    OUT.mkdir(parents=True, exist_ok=True)
    cols = COLUMNS + ['date', 'zone', 'power', 'quality', 'pft_name', 'forest', 'path_dist_m', 'path_source']
    save_rows(rows, OUT / 'gedi_v3_all_points.csv', cols)
    good = [r for r in rows if r['quality']]
    save_rows(good, OUT / 'gedi_v3_quality_points.csv', cols)
    sites = pick_sites(good)
    save_rows(sites, OUT / 'field_sites.csv',
              ['site_id', 'zone', 'lat', 'lon', 'pft_name', 'agbd', 'agbd_se', 'landsat_treecover', 'date', 'beam',
               'path_dist_m', 'path_source', 'google_maps', 'shot_number'])
    write_kml(sites, OUT / 'field_sites.kml')
    write_gpx(sites, OUT / 'field_sites.gpx')
    figure(rows, good, sites)
    summary = []
    for zone in ZONES:
        z = [r for r in rows if r['zone'] == zone]
        q = [r for r in z if r['quality']]
        f = [r for r in q if r['forest']]
        summary.append({'zone': zone, '전체': len(z), '품질통과': len(q), '품질통과_숲': len(f),
                        '평균AGBD_품질통과': round(float(np.mean([r['agbd'] for r in q])), 1) if q else '',
                        '평균AGBD_숲': round(float(np.mean([r['agbd'] for r in f])), 1) if f else '',
                        '관측연도': ', '.join(sorted({r['date'][:4] for r in z})),
                        '품질통과_관측연도': ', '.join(sorted({r['date'][:4] for r in q})),
                        '수관피복_기준': 20 if zone == '부산대 부산캠퍼스' else 30,
                        '현장후보': sum(r['zone'] == zone for r in sites)})
    save_rows(summary, OUT / 'zone_summary.csv', list(summary[0]))
    run_info['trails'] = TRAIL_INFO
    save_json(run_info, OUT / 'run_metadata.json')
    write_report(summary, sites, run_info)
    for s in summary:
        print(s)
    if summary[0]['품질통과'] == 0:
        print('부산대 부산캠퍼스에는 품질 통과 GEDI L4A V3 관측점이 없습니다.')
    print(f'현장 방문 후보: {len(sites)}곳 → outputs/field_sites.csv, .kml, .gpx, gedi_sites_map.png')


def write_report(summary, sites, run_info):
    lines = ['# 금정산·부산대 부산캠퍼스 GEDI L4A V3 조사', '',
             f'검색·처리: {run_info["completed_granules"]}/{run_info["total_granules"]}개 궤도 파일.',
             f'범위: {BBOX}, 기간: 2019–2025년. 구역은 캠퍼스 BBOX와 나머지 범위로 상호 배타적으로 분류.',
             '품질: l4a_quality_flag_rel3 = 1, degrade_flag = 0, 유한한 AGBD ≥ 0.',
             '숲 위 점: 품질 통과 중 Landsat 수관피복이 금정산 ≥30%, 캠퍼스 ≥20%인 점.',
             '평균 AGBD 단위는 Mg/ha. 아래 표는 품질 통과 점 평균과 그중 숲 위 점 평균을 구분한다. '
             '관측점의 단순 평균이며 구역 전체의 면적 가중 평균이 아니다.',
             f'등산로: {TRAIL_INFO.get("source", "없음")}. EPSG:5179에서 선분까지의 최단 직선거리.', '',
             '| 구역 | 전체 점 | 품질 통과 | 숲 위 점 | 평균 AGBD(품질 통과) | 평균 AGBD(숲) | 관측 연도 |',
             '|---|---:|---:|---:|---:|---:|---|']
    for s in summary:
        lines.append(f'| {s["zone"]} | {s["전체"]} | {s["품질통과"]} | {s["품질통과_숲"]} | '
                     f'{s["평균AGBD_품질통과"]} | {s["평균AGBD_숲"]} | {s["관측연도"]} |')
    if summary[0]['품질통과'] == 0:
        lines += ['', '**부산대 부산캠퍼스에는 품질 통과 관측점이 없습니다.**']
    lines += ['', '금정산 후보: 등산로 ≤150 m, 후보 간 WGS84 측지거리 ≥300 m, 최대 15곳.',
              '캠퍼스 후보: 등산로·최소 간격 조건 없이 최대 10곳. 파워 빔과 야간 관측 우선.', '',
              '| 번호 | 구역 | 위도 | 경도 | 지도 | AGBD | 등산로 거리(m) |',
              '|---|---|---:|---:|---|---:|---:|']
    for r in sites:
        distance = f'{r["path_dist_m"]:.1f}' if math.isfinite(r['path_dist_m']) else '미확인'
        lines.append(f'| {r["site_id"]} | {r["zone"]} | {r["lat"]:.6f} | {r["lon"]:.6f} | '
                     f'[구글맵]({r["google_maps"]}) | {r["agbd"]:.1f} | {distance} |')
    lines += ['', 'AGBD는 각 관측 시점의 약 25 m 지름 footprint 추정값이며 현재 현장 수목량과 다를 수 있다.',
              'Landsat 수관피복은 2010년 자료를 바탕으로 한 GEDI 보조자료이며 현장 방문 시점의 산림 경계가 아니다.',
              '등산로 거리는 직선거리로, 실제 진입 가능 여부·보행거리·경사는 현장에서 확인해야 한다.', '',
              '자료 설명: [ORNL DAAC GEDI L4A V3](https://daac.ornl.gov/GEDI/guides/GEDI_L4A_AGB_Density_V3.html)', '',
              '## 처리량', '',
              'granule_metrics.csv의 remote_bytes는 range 응답 본문 바이트 수(블록 여유분 포함), '
              'dataset_bytes는 압축 해제 후 읽은 배열 바이트 수, download_bytes는 대체 다운로드 파일 크기다. '
              'HTTP/TLS 헤더·인증·검색 트래픽은 포함하지 않는다.', '']
    metrics = run_info['metrics']
    lines += [f'- 원격 range 응답: {sum(r["remote_bytes"] for r in metrics) / 2**30:.3f} GiB',
              f'- 대체 전체 다운로드: {sum(r["download_bytes"] for r in metrics) / 2**30:.3f} GiB',
              f'- 전체 파일 크기 합계: {sum(r["file_bytes"] for r in metrics) / 2**30:.3f} GiB',
              f'- 이번 이어받기 실행 경과: {run_info.get("collection_wall_seconds", 0) / 60:.1f}분 '
              f'({run_info.get("workers", 1)}개 독립 프로세스)', '']
    trial_path = OUT / 'trial' / 'run_metadata.json'
    if trial_path.exists():
        trial_metrics = json.loads(trial_path.read_text(encoding='utf-8'))['metrics']
        lines += ['## 먼저 수행한 2개 파일 시험', '',
                  '| 파일 | 초 | 원격 수신량(MiB) | 압축 해제 배열(MiB) | 범위 내 점 |',
                  '|---|---:|---:|---:|---:|']
        for r in trial_metrics:
            lines.append(f'| {r["granule"]} | {r["seconds"]:.1f} | {r["remote_bytes"] / 2**20:.2f} | '
                         f'{r["dataset_bytes"] / 2**20:.2f} | {r["shots"]} |')
        if trial_metrics:
            avg = sum(r['seconds'] for r in trial_metrics) / len(trial_metrics)
            lines += ['', f'시험 평균 {avg:.1f}초/파일 기준 83개 순차 처리 예상은 {83 * avg / 60:.1f}분이었다.', '']
    lines += ['## 생성 파일', '']
    lines += [f'- {p.name}' for p in sorted(OUT.iterdir()) if p.is_file() and p.name != 'report.md']
    lines += ['- report.md']
    (OUT / 'report.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
