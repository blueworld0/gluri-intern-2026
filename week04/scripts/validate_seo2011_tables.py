"""Check final CSV relations, source-cell traceability, units and KML values."""
from collections import Counter
import csv
import json
import math
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import extract_seo2011_fieldnotes as m

sys.stdout.reconfigure(encoding='utf-8')
def read(name):
    with (m.FIELD / name).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

plots = read('seo2011_plots.csv')
trees = read('seo2011_trees.csv')
summaries = read('seo2011_plot_summary.csv')
sources = read('seo2011_source_index.csv')
provenance = read('seo2011_tree_provenance.csv')
layers = read('seo2011_layers.csv')
abnormal = read('seo2011_abnormal_trees.csv')
expected = {f'{z}-{i:02}' for z in m.ZONES for i in range(1, 11)}
assert len(plots) == len(summaries) == len(sources) == 50
assert {r['plot_id'] for r in plots} == expected
assert len(trees) == len(provenance) == 1750
assert sum(int(p['나무 수']) for p in plots) == len(trees)
assert len(layers) == 200
assert not abnormal
assert set(trees[0]) == set(m.TREE_COLUMNS)
assert set(plots[0]) == set(m.PLOT_COLUMNS)
for tree, origin in zip(trees, provenance):
    assert origin['tree_csv_row'] == str(provenance.index(origin) + 2)
    assert tree['plot_id'] == origin['plot_id']
    cell = m.TRIPLE.fullmatch(origin['raw_cell'])
    assert cell and (tree['수종'], float(tree['흉고직경_cm']), float(tree['수고_m'])) == (
        cell.group(1), float(cell.group(2)), float(cell.group(3)))
    assert 0 < float(tree['흉고직경_cm']) <= 100 and 0 < float(tree['수고_m']) <= 40
for plot, summary, source in zip(plots, summaries, sources):
    assert plot['plot_id'] == summary['plot_id'] == source['plot_id']
    subset = [t for t in trees if t['plot_id'] == plot['plot_id']]
    assert len(subset) == int(plot['나무 수']) == int(source['매목_본수_plain']) == int(source['매목_본수_layout'])
    assert int(source['pdf_page']) - int(source['printed_page']) == 12
    for field, raw in [('위도', '위도_원문'), ('경도', '경도_원문')]:
        d, minutes, seconds = map(float, re.findall(r'\d+(?:\.\d+)?', source[raw]))
        assert math.isclose(float(plot[field]), (d*3600 + minutes*60 + seconds)/3600, abs_tol=1e-12)
    diameter = np.array([float(r['흉고직경_cm']) for r in subset])
    heights = np.array([float(r['수고_m']) for r in subset])
    assert float(summary['1ha당_본수']) == len(subset) * 25
    assert math.isclose(float(summary['평균_흉고직경_cm']), diameter.mean(), abs_tol=1e-10)
    assert math.isclose(float(summary['평균_수고_m']), heights.mean(), abs_tol=1e-10)
    assert math.isclose(float(summary['흉고단면적_합계_m2_ha']), np.sum(np.pi*diameter**2/40000)*25, abs_tol=1e-10)
    group = [r for r in layers if r['plot_id'] == plot['plot_id']]
    assert len(group) == 4 and len({r['층'] for r in group}) == 4
    assert all(0 <= float(r['식피율_pct']) <= 100 for r in group)

root = ET.parse(m.OUT / 'seo2011_plots.kml')
ns = {'k': 'http://www.opengis.net/kml/2.2'}
marks = root.findall('.//k:Placemark', ns)
assert len(marks) == 50
lookup = {r['plot_id']: r for r in plots}
summary_lookup = {r['plot_id']: r for r in summaries}
assert {mark.findtext('k:name', namespaces=ns) for mark in marks} == expected
for mark in marks:
    name = mark.findtext('k:name', namespaces=ns)
    lon, lat, _ = map(float, mark.findtext('k:Point/k:coordinates', namespaces=ns).split(','))
    assert abs(lon-float(lookup[name]['경도'])) < 1e-9 and abs(lat-float(lookup[name]['위도'])) < 1e-9
    description = mark.findtext('k:description', namespaces=ns)
    assert summary_lookup[name]['우점_수종'] in description
    assert f'{float(summary_lookup[name]["평균_흉고직경_cm"]):.2f} cm' in description
assert not {'pandas', 'geopandas', 'rasterio', 'pyogrio'} & set(sys.modules)
metadata_path = m.FIELD / 'seo2011_validation.json'
meta = json.loads(metadata_path.read_text(encoding='utf-8'))
meta['output_checks'] = 'PASS'
meta['kml_placemarks'] = len(marks)
metadata_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'validation': 'PASS', 'plots': len(plots), 'GPS': len(sources), 'trees': len(trees),
                  'source_cell_matches': len(provenance), 'layers': len(layers), 'abnormal_trees': len(abnormal),
                  'KML_placemarks': len(marks), 'zone_tree_counts': meta['zone_tree_counts']}, ensure_ascii=False))
