"""Round-trip SHP, area-conservation, point-attribute and visit-list checks."""
import sys
import json
import math
import xml.etree.ElementTree as ET
import shapefile
import shapely
from shapely.geometry import shape,Point
from pyproj import CRS,Geod
from forest_map_common import *

sys.stdout.reconfigure(encoding='utf-8')
meta=json.loads((OUT/'forest_analysis_metadata.json').read_text(encoding='utf-8'))
geometries={}; attrs_by_year={}; errors={}
for key,stem in [('2013','forest2013_geumjeong'),('2025','forest2025_geumjeong'),('trails','trails_geumjeong')]:
    path=MERGED/(stem+'.shp')
    assert all(path.with_suffix(suffix).exists() for suffix in ('.shp','.dbf','.shx','.prj','.cpg'))
    assert CRS.from_wkt(path.with_suffix('.prj').read_text()).to_epsg()==5179
    with shapefile.Reader(str(path),encoding='cp949') as reader:
        fields=[f[0] for f in reader.fields[1:]]
        features=[(shape(sr.shape.__geo_interface__),sr.record.as_dict()) for sr in reader.iterShapeRecords()]
    geometries[key]=[g for g,a in features]; attrs_by_year[key]=[a for g,a in features]
    assert len(features)==(meta['trail_parts'] if key=='trails' else meta['layer_counts'][key])
    assert len({a['FID_U'] for g,a in features})==len(features)
    assert all(g.is_valid and not g.is_empty for g,a in features)
    assert all(g.difference(CLIP.buffer(.00001)).is_empty for g,a in features)
    if key!='trails':
        assert all(name in fields for name in ('FRTP_NM','KOFTR_NM','DMCLS_NM','AGCLS_NM','HEIGHT','갱신년도'))
        assert all(math.isclose(g.area,a['Shape_Area'],abs_tol=.001) for g,a in features)
        if key=='2013':
            assert all(not str(a['HEIGHT']).strip() for a in attrs_by_year[key])
        union=shapely.union_all(geometries[key])
        duplicate_area=sum(g.area for g in geometries[key])-union.area
        assert abs(duplicate_area)<.01,(key,duplicate_area)
        errors[key+'_duplicate_area_m2']=duplicate_area

areas=read_csv(OUT/'forest_area_summary.csv'); coverage=read_csv(OUT/'forest_coverage_summary.csv')
for c in coverage:
    year=c['year']; region=CAMPUS_GEOM if c['scope']=='부산대 캠퍼스' else CLIP
    actual=shapely.union_all(geometries[year]).intersection(region).area/10000
    assert math.isclose(actual,float(c['mapped_area_ha']),abs_tol=1e-6)
    assert math.isclose(float(c['mapped_area_ha'])+float(c['unmapped_area_ha']),float(c['region_area_ha']),abs_tol=1e-6)
    for dimension in ('임상','수종','경급','영급'):
        total=sum(float(r['area_ha']) for r in areas if r['year']==year and r['scope']==c['scope'] and r['dimension']==dimension)
        assert math.isclose(total,actual,abs_tol=1e-6),(year,c['scope'],dimension,total,actual)

changes=read_csv(OUT/'forest_change_2013_2025.csv')
common=shapely.union_all(geometries['2013']).intersection(shapely.union_all(geometries['2025']))
for scope,region in [('금정산·캠퍼스 전체',CLIP),('부산대 캠퍼스',CAMPUS_GEOM)]:
    recorded=sum(float(r['area_ha']) for r in changes if r['scope']==scope)
    assert math.isclose(recorded,common.intersection(region).area/10000,abs_tol=1e-5)

comparisons=read_csv(OUT/'seo2011_vs_forestmap.csv')
assert len(comparisons)==50 and len({r['plot_id'] for r in comparisons})==50
for year in ('2013','2025'):
    p=f'fm{year}_'; tree=shapely.STRtree(geometries[year])
    forest=[g for g,a in zip(geometries[year],attrs_by_year[year]) if stocked(a)]
    forest_tree=shapely.STRtree(forest)
    attrs={int(a['FID_U']):a for a in attrs_by_year[year]}
    by_id={int(a['FID_U']):g for g,a in zip(geometries[year],attrs_by_year[year])}
    for row in comparisons:
        point=Point(*TO_XY.transform(float(row['lon']),float(row['lat'])))
        fid=row[p+'polygon_id']
        if fid:
            a=attrs[int(fid)]
            assert by_id[int(fid)].intersects(point)
            assert row[p+'species_code']==a['KOFTR_GROU']
            assert row[p+'in_stocked_forest']==str(stocked(a))
        else:
            assert not len(tree.query(point,predicate='intersects'))
        nearest=forest[int(forest_tree.nearest(point))]
        assert math.isclose(float(row[p+'nearest_forest_distance_m']),point.distance(nearest),abs_tol=1e-6)

visits=read_csv(OUT/'field_visit_combined.csv'); candidates=read_csv(OUT/'field_sites.csv')
assert len(visits)==75 and len({r['visit_id'] for r in visits})==75 and len(candidates)==25
original=read_csv(OUT/'field_sites_gedi_original.csv')
assert [r['site_id'] for r in original]==[r['site_id'] for r in candidates]
for a,b in zip(original,candidates):
    assert all(a[key]==b[key] for key in a),a['site_id']
assert sum(r['record_type']=='야장 조사구' for r in visits)==50
assert sum(r['record_type']=='GEDI 후보' for r in visits)==25
root=ET.parse(OUT/'field_visit_combined.kml'); ns={'k':'http://www.opengis.net/kml/2.2'}
marks=root.findall('.//k:Placemark',ns)
assert len(marks)==75
assert {m.findtext('k:name',namespaces=ns) for m in marks}=={r['visit_id'] for r in visits}
geod=Geod(ellps='WGS84')
plots={r['visit_id']:r for r in visits if r['record_type']=='야장 조사구'}
gedi={r['site_id']:r for r in candidates}; pairs=read_csv(OUT/'field_visit_pairs_100m.csv')
expected_pairs=set()
for name,p in plots.items():
    for gid,g in gedi.items():
        d=geod.inv(float(p['lon']),float(p['lat']),float(g['lon']),float(g['lat']))[2]
        if d<=100:
            expected_pairs.add((name,gid))
assert {(r['plot_id'],r['gedi_site_id']) for r in pairs}==expected_pairs
for r in pairs:
    p,g=plots[r['plot_id']],gedi[r['gedi_site_id']]
    d=geod.inv(float(p['lon']),float(p['lat']),float(g['lon']),float(g['lat']))[2]
    assert math.isclose(d,float(r['distance_m']),abs_tol=1e-8)
assert not {'pandas','geopandas','rasterio','pyogrio'}&set(sys.modules)
result={'validation':'PASS','polygon_parts':meta['layer_counts'],'trails':len(geometries['trails']),
        'plot_comparisons':50,'combined_visits':75,'kml_placemarks':75,
        'candidate_plot_pairs_100m':len(pairs),'roundtrip_area_errors':errors}
save_json(OUT/'forest_analysis_validation.json',result)
print(json.dumps(result,ensure_ascii=False,indent=2))
