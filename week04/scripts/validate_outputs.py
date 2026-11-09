"""Verify completed scientific outputs without rereading the remote files."""
import csv
import json
import math
import sys
import xml.etree.ElementTree as ET
from itertools import combinations

from pyproj import Geod
import find_gedi_sites as m

def read(name):
    with (m.OUT / name).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

meta = json.loads((m.OUT / 'run_metadata.json').read_text(encoding='utf-8'))
assert meta['complete'] and meta['completed_granules'] == meta['total_granules'] == 83
assert meta['trails']['source'] == '산림청 등산로'
manifest = json.loads((m.CACHE / 'manifest.json').read_text(encoding='utf-8'))
assert len(meta['metrics']) == 83
assert {r['granule'] for r in meta['metrics']} == set(manifest['granules'])
assert all((m.CACHE / 'granules' / (name + '.csv')).exists() for name in manifest['granules'])
rows = read('gedi_v3_all_points.csv')
quality = read('gedi_v3_quality_points.csv')
sites = read('field_sites.csv')
summary = read('zone_summary.csv')
assert len({r['shot_number'] for r in rows}) == len(rows)
assert sum(int(s['전체']) for s in summary) == len(rows)
assert sum(int(s['품질통과']) for s in summary) == len(quality)
for row in rows:
    lat, lon = float(row['lat']), float(row['lon'])
    assert m.BBOX[0] <= lon <= m.BBOX[2] and m.BBOX[1] <= lat <= m.BBOX[3]
    assert row['zone'] == m.zone_of(lat, lon)
for row in quality:
    assert float(row['l4a_quality_flag_rel3']) == 1 and float(row['degrade_flag']) == 0
    assert math.isfinite(float(row['agbd'])) and float(row['agbd']) >= 0
for s in summary:
    z = [r for r in rows if r['zone'] == s['zone']]
    q = [r for r in quality if r['zone'] == s['zone']]
    threshold = 20 if s['zone'] == '부산대 부산캠퍼스' else 30
    forest = [r for r in q if r['landsat_treecover'] and threshold <= float(r['landsat_treecover']) <= 100]
    assert (int(s['전체']), int(s['품질통과']), int(s['품질통과_숲'])) == (len(z), len(q), len(forest))
    if q:
        assert round(sum(float(r['agbd']) for r in q)/len(q), 1) == float(s['평균AGBD_품질통과'])
    if forest:
        assert round(sum(float(r['agbd']) for r in forest)/len(forest), 1) == float(s['평균AGBD_숲'])
    assert int(s['현장후보']) == sum(r['zone'] == s['zone'] for r in sites)
    assert s['관측연도'] == ', '.join(sorted({r['date'][:4] for r in z}))
good_ids = {r['shot_number'] for r in quality}
mountain, campus = [], []
for r in sites:
    assert r['shot_number'] in good_ids
    assert r['google_maps'].endswith(f'{float(r["lat"]):.6f},{float(r["lon"]):.6f}')
    if r['zone'] == '금정산 일대':
        mountain.append(r)
        assert 30 <= float(r['landsat_treecover']) <= 100
        assert 0 <= float(r['path_dist_m']) <= 150
    else:
        campus.append(r)
        assert r['zone'] == '부산대 부산캠퍼스' and 20 <= float(r['landsat_treecover']) <= 100
assert len(mountain) <= 15 and len(campus) <= 10
geod = Geod(ellps='WGS84')
distances = [geod.inv(float(a['lon']), float(a['lat']), float(b['lon']), float(b['lat']))[2]
             for a, b in combinations(mountain, 2)]
assert not distances or min(distances) >= 300
kml = ET.parse(m.OUT / 'field_sites.kml')
gpx = ET.parse(m.OUT / 'field_sites.gpx')
assert len(kml.findall('.//{http://www.opengis.net/kml/2.2}Placemark')) == len(sites)
assert len(gpx.findall('.//{http://www.topografix.com/GPX/1/1}wpt')) == len(sites)
assert not {'pandas', 'geopandas', 'rasterio', 'pyogrio'} & set(sys.modules)
print(json.dumps({'validation': 'PASS', 'all_points': len(rows), 'quality_points': len(quality),
                  'mountain_sites': len(mountain), 'campus_sites': len(campus),
                  'min_mountain_spacing_m': round(min(distances), 3) if distances else None,
                  'remote_MiB': sum(r['remote_bytes'] for r in meta['metrics'])/2**20,
                  'downloaded_MiB': sum(r['download_bytes'] for r in meta['metrics'])/2**20}, ensure_ascii=False))
