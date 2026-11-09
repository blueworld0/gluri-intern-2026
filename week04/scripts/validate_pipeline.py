"""Offline checks for spatial rules, precision, local HDF5 reading and trails."""
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
from pyproj import Geod, Transformer
from shapely.geometry import LineString

import find_gedi_sites as m

assert m.zone_of(35.234, 129.078) == '부산대 부산캠퍼스'
assert m.zone_of(35.205, 129.005) == '금정산 일대'
assert m.zone_of(35.27, 129.15) == '그 외'
trails = m.load_paths()
assert m.TRAIL_INFO['source'] == '산림청 등산로', m.TRAIL_INFO
assert len(m.TRAIL_INFO['files']) == 3
transform = Transformer.from_crs(4326, 5179, always_xy=True)
x, y = transform.transform(129.05, 35.28)
d = m.path_distance(np.array([35.28]), np.array([129.05]),
                    [LineString([(x - 1000, y + 100), (x + 1000, y + 100)])])
assert abs(d[0] - 100) < 1e-6, d

with tempfile.TemporaryDirectory(dir=m.ROOT) as tmp:
    path = Path(tmp) / 'synthetic.h5'
    shot = 12345678901234567  # must never round through a float
    with h5py.File(path, 'w') as f:
        beam = f.create_group('BEAM0101')
        for key in m.VARS:
            values = {'shot_number': np.array([shot, shot+1], dtype='uint64'),
                      'lat_lowestmode': [35.28, 35.5], 'lon_lowestmode': [129.05, 129.05],
                      'l4a_quality_flag_rel3': [1, 1], 'degrade_flag': [0, 0]}.get(key, [50., 50.])
            beam.create_dataset(key, data=values)
        lc = beam.create_group('land_cover_data')
        lc.create_dataset('landsat_treecover', data=[30., 0.])
        lc.create_dataset('pft_class', data=[4, 4])
    stats = {'dataset_bytes': 0}
    rows = m.read_granule(path, 'synthetic', stats)
    assert len(rows) == 1 and rows[0]['shot_number'] == shot
    m.save_rows(rows, Path(tmp) / 'roundtrip.csv', m.COLUMNS)
    assert m.load_cache(Path(tmp) / 'roundtrip.csv')[0]['shot_number'] == shot
    assert stats['dataset_bytes'] > 0

def candidate(i, zone, cover, distance, lon=129.05, quality=True):
    return dict(shot_number=i, zone=zone, forest=cover >= (20 if zone == '부산대 부산캠퍼스' else 30),
                quality=quality, path_dist_m=distance, power=True, solar_elevation=-1,
                lat=35.28, lon=lon, agbd=5., agbd_se=2.)

mountain = '금정산 일대'
campus = '부산대 부산캠퍼스'
inputs = [candidate(1, mountain, 30, 150), candidate(2, mountain, 30, 150.0001, 129.1),
          candidate(3, mountain, 30, float('nan'), 129.1), candidate(4, mountain, 29, 0, 129.1),
          candidate(5, mountain, 30, 140, 129.0501), candidate(6, mountain, 30, 0, 129.1, False),
          candidate(7, campus, 20, float('nan')), candidate(8, campus, 19, 0)]
sites = m.pick_sites(inputs)
assert {r['shot_number'] for r in sites} == {5, 7}, sites
geod = Geod(ellps='WGS84')
separated = m.pick_sites([candidate(i, mountain, 30, 50, 129.01 + i * .01) for i in range(20)])
assert len(separated) == 15
assert min(geod.inv(a['lon'], a['lat'], b['lon'], b['lat'])[2]
           for i, a in enumerate(separated) for b in separated[i+1:]) >= 300
assert len(m.pick_sites([candidate(i, campus, 20, 9999) for i in range(20)])) == 10
assert not {'pandas', 'geopandas', 'rasterio', 'pyogrio'} & set(sys.modules)
print('PASS: trail files, segment distance, exclusive zones, HDF5/CSV shot precision, candidate rules and limits')
