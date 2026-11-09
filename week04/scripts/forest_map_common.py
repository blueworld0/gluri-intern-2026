"""Forest-map I/O and documented code tables; no GDAL or dataframe libraries."""
from pathlib import Path
from collections import defaultdict
import csv
import json
import math
import warnings

import numpy as np
import shapefile
import shapely
from shapely.geometry import shape, box, Polygon, LineString
from shapely.geometry.polygon import orient
from shapely.ops import transform
from pyproj import CRS, Transformer

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/forest_map'
MERGED = RAW / 'merged'
OUT = ROOT / 'outputs'
BBOX = (129.00, 35.20, 129.12, 35.34)
CAMPUS = (129.0738, 35.2299, 129.0847, 35.2386)
TO_XY = Transformer.from_crs(4326, 5179, always_xy=True)
TO_LL = Transformer.from_crs(5179, 4326, always_xy=True)
CLIP = box(*TO_XY.transform_bounds(*BBOX, densify_pts=41))

def geographic_box(bounds):
    w, s, e, n = bounds
    k = np.linspace(0, 1, 41)
    ring = [(w+(e-w)*v,s) for v in k] + [(e,s+(n-s)*v) for v in k]
    ring += [(e-(e-w)*v,n) for v in k] + [(w,n-(n-s)*v) for v in k]
    return transform(TO_XY.transform, Polygon(ring))

CAMPUS_GEOM = geographic_box(CAMPUS)
# DATA003 definition, PDF pages 4-8. Names in *_NM are retained when present.
FRTP = {'0':'무립목지/비산림','1':'침엽수림','2':'활엽수림','3':'혼효림','4':'죽림'}
SPECIES = {'11':'소나무','12':'잣나무','13':'낙엽송','14':'리기다소나무','15':'곰솔',
 '16':'잔나무','17':'편백나무','18':'삼나무','19':'가문비나무','20':'비자나무','21':'은행나무',
 '10':'기타침엽수','31':'상수리나무','32':'신갈나무','33':'굴참나무','34':'기타참나무류',
 '35':'오리나무','36':'고로쇠나무','37':'자작나무','38':'박달나무','39':'밤나무','40':'물푸레나무',
 '41':'서어나무','42':'때죽나무','43':'호두나무','44':'백합나무','45':'포플러','46':'벚나무',
 '47':'느티나무','48':'층층나무','49':'아까시나무','30':'기타활엽수','61':'가시나무',
 '62':'구실잣밤나무','63':'녹나무','64':'굴거리나무','65':'황칠나무','66':'사스레피나무',
 '67':'후박나무','68':'새덕이','60':'기타상록활엽수','77':'침활혼효림','78':'죽림',
 '81':'미립목지','82':'제지','83':'관목덤불','91':'주거지','92':'초지','93':'경작지',
 '94':'수체','95':'과수원','99':'기타'}
DIAMETER = {'0':'치수','1':'소경목','2':'중경목','3':'대경목'}
AGE = {str(i): f'{i}영급' for i in range(1,10)}
DENSITY = {'A':'소','B':'중','C':'밀'}
MAPS = [('FRTP_CD','FRTP_NM',FRTP),('KOFTR_GROU','KOFTR_NM',SPECIES),
        ('DMCLS_CD','DMCLS_NM',DIAMETER),('AGCLS_CD','AGCLS_NM',AGE),('DNST_CD','DNST_NM',DENSITY)]

def clean(value):
    return '' if value is None else str(value).strip()

def stocked(attrs):
    return clean(attrs.get('FRTP_CD')) in ('1','2','3','4')

def parts(geom, kind='Polygon'):
    if geom.is_empty:
        return []
    if geom.geom_type == kind:
        return [geom]
    if hasattr(geom, 'geoms'):
        return [p for g in geom.geoms for p in parts(g, kind)]
    return []

def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def save_csv(path, rows, columns=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(dict.fromkeys(key for row in rows for key in row))
    pending = path.with_suffix(path.suffix+'.tmp')
    with pending.open('w',encoding='utf-8-sig',newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    pending.replace(path)

def save_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')

def normalize(attrs):
    for code, label, mapping in MAPS:
        value = clean(attrs.get(code))
        if value and value not in mapping:
            raise ValueError(f'Unknown code {code}={value}')
        attrs[code] = value
        attrs[label] = clean(attrs.get(label)) or mapping.get(value,'')
    attrs.setdefault('HEIGHT','')
    attrs.setdefault('HEIGHT_NM','')
    attrs.setdefault('갱신년도','')
    return attrs

def load_layer(period, kind='Polygon'):
    features, schema, log = [], {}, []
    for district in ('금정구','북구','양산시'):
        files = sorted((RAW / period / district).glob('*.shp'))
        if not files:
            raise FileNotFoundError(f'{period}/{district}')
        for path in files:
            crs = CRS.from_wkt(path.with_suffix('.prj').read_text(encoding='utf-8-sig'))
            if crs.to_epsg() != 5179:
                raise ValueError(f'Expected EPSG:5179: {path}')
            counts = dict(path=str(path.relative_to(ROOT)), input_records=0, bbox_candidates=0,
                          clipped_parts=0, repaired_geometries=0, dropped_area_m2=0.)
            with shapefile.Reader(str(path),encoding='cp949') as reader:
                counts['input_records'] = len(reader)
                for field in reader.fields[1:]:
                    schema[field[0]] = tuple(field)
                for sr in reader.iterShapeRecords(bbox=CLIP.bounds):
                    counts['bbox_candidates'] += 1
                    geom = shape(sr.shape.__geo_interface__)
                    if not geom.is_valid:
                        counts['repaired_geometries'] += 1
                        geom = shapely.make_valid(geom)
                    geom = geom.intersection(CLIP)
                    attrs = sr.record.as_dict()
                    if kind == 'Polygon':
                        attrs = normalize(attrs)
                    for part_id, part in enumerate(parts(geom,kind)):
                        if part.is_empty or (part.area <= 1e-6 if kind=='Polygon' else part.length <= 1e-6):
                            continue
                        record = dict(attrs, SRC_DIST=district, SRC_ID=int(sr.record.oid), PART_ID=part_id)
                        features.append({'geometry':part, 'attrs':record})
                        counts['clipped_parts'] += 1
            log.append(counts)
    if kind == 'Polygon':
        # Remove positive-area overlaps in a stable district/record order.
        # Touching boundaries are retained. Only cleaned geometry is analyzed.
        original = [f['geometry'] for f in features]
        tree = shapely.STRtree(original)
        cleaned, overlap_area, overlap_details = [], 0., []
        accepted = {}
        for i, feature in enumerate(features):
            geom = feature['geometry']
            neighbor_ids = [int(j) for j in tree.query(geom, predicate='intersects') if j < i and j in accepted]
            neighbors = [accepted[j] for j in neighbor_ids]
            if neighbors:
                cut = shapely.union_all(neighbors)
                overlap = geom.intersection(cut).area
                if overlap > 1e-6:
                    for j in neighbor_ids:
                        shared = geom.intersection(accepted[j]).area
                        if shared <= 1e-6:
                            continue
                        a,b = features[j]['attrs'],feature['attrs']
                        overlap_details.append(dict(layer=period,kept_district=a['SRC_DIST'],kept_src_id=a['SRC_ID'],
                            removed_district=b['SRC_DIST'],removed_src_id=b['SRC_ID'],area_m2=shared,
                            kept_type=a['FRTP_CD'],removed_type=b['FRTP_CD'],
                            kept_species=a['KOFTR_GROU'],removed_species=b['KOFTR_GROU'],
                            classification_conflict=(a['FRTP_CD'],a['KOFTR_GROU'],a['DMCLS_CD'],a['AGCLS_CD']) !=
                                                    (b['FRTP_CD'],b['KOFTR_GROU'],b['DMCLS_CD'],b['AGCLS_CD'])))
                    geom = geom.difference(cut)
                    overlap_area += overlap
            accepted[i] = geom
            for subpart, part in enumerate(parts(geom)):
                if part.area > 1e-6:
                    cleaned.append({'geometry':part,'attrs':dict(feature['attrs'], PART_ID=subpart)})
        features = cleaned
        log.append({'positive_overlap_removed_m2':overlap_area,
                    'overlap_priority':'금정구→북구→양산시, 각 원본 레코드 순; 앞선 속성 유지',
                    'conflicting_overlap_m2':sum(r['area_m2'] for r in overlap_details if r['classification_conflict'])})
        save_csv(OUT/f'forest_overlap_{period[:4]}.csv',overlap_details,
                 ['layer','kept_district','kept_src_id','removed_district','removed_src_id','area_m2',
                  'kept_type','removed_type','kept_species','removed_species','classification_conflict'])
    for i, feature in enumerate(features):
        feature['attrs']['FID_U'] = i+1
        feature['attrs']['Shape_Leng'] = feature['geometry'].length
        if kind == 'Polygon':
            feature['attrs']['Shape_Area'] = feature['geometry'].area
    print(f'{period}: {len(features):,} clipped {kind} features',flush=True)
    return features, schema, log

def write_layer(features, schema, path, kind='Polygon'):
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    schema = dict(schema)
    if kind == 'Polygon':
        for _, name, _ in MAPS:
            schema.setdefault(name,(name,'C',120,0))
        schema.setdefault('HEIGHT',('HEIGHT','C',2,0))
        schema.setdefault('HEIGHT_NM',('HEIGHT_NM','C',120,0))
        schema.setdefault('갱신년도',('갱신년도','C',8,0))
    schema.update(SRC_DIST=('SRC_DIST','C',20,0),SRC_ID=('SRC_ID','N',10,0),
                  PART_ID=('PART_ID','N',8,0),FID_U=('FID_U','N',10,0))
    with shapefile.Writer(str(path),shapeType=shapefile.POLYGON if kind=='Polygon' else shapefile.POLYLINE,
                          encoding='cp949') as writer:
        for field in schema.values():
            writer.field(*field)
        for feature in features:
            geom, attrs = feature['geometry'], feature['attrs']
            if kind == 'Polygon':
                poly = orient(geom,sign=-1)
                writer.poly([list(poly.exterior.coords)]+[list(r.coords) for r in poly.interiors])
            else:
                writer.line([list(geom.coords)])
            writer.record(*[attrs.get(key,'') for key in schema])
    path.with_suffix('.prj').write_text(CRS.from_epsg(5179).to_wkt(version='WKT1_ESRI'),encoding='utf-8')
    path.with_suffix('.cpg').write_text('CP949',encoding='ascii')

class ForestIndex:
    def __init__(self,features):
        self.features = features
        self.geoms = np.array([f['geometry'] for f in features],dtype=object)
        self.tree = shapely.STRtree(self.geoms)
        self.forest_ids = np.array([i for i,f in enumerate(features) if stocked(f['attrs'])],dtype=int)
        self.forest_tree = shapely.STRtree(self.geoms[self.forest_ids])

    def locate(self, point):
        hits = sorted(int(i) for i in self.tree.query(point,predicate='intersects'))
        # Boundaries are deterministic, with stocked forest preferred at a tie.
        hits.sort(key=lambda i:(not stocked(self.features[i]['attrs']),i))
        inside = hits[0] if hits else None
        near = int(self.forest_ids[int(self.forest_tree.nearest(point))])
        distance = float(point.distance(self.geoms[near]))
        return inside, near, distance, len(hits)

def diameter_code(value):
    return '0' if value < 6 else '1' if value < 18 else '2' if value < 30 else '3'

def height_bounds(code):
    if not clean(code):
        return None
    h = int(code)
    if h == 0:
        return 0.,1.
    return (float(h-1),math.inf if h==40 else float(h+1))

def attach_info(index,point,prefix):
    inside, nearest, distance, hits = index.locate(point)
    attrs = index.features[inside]['attrs'] if inside is not None else {}
    near = index.features[nearest]['attrs']
    result = {prefix+'polygon_id':attrs.get('FID_U',''),prefix+'in_stocked_forest':stocked(attrs),
              prefix+'covering_polygons':hits,prefix+'nearest_forest_distance_m':distance,
              prefix+'nearest_forest_id':near['FID_U'],prefix+'nearest_forest_species':near['KOFTR_NM']}
    for source,target in [('FRTP_CD','type_code'),('FRTP_NM','type'),('KOFTR_GROU','species_code'),
                          ('KOFTR_NM','species'),('DMCLS_CD','diameter_code'),('DMCLS_NM','diameter_class'),
                          ('AGCLS_CD','age_code'),('AGCLS_NM','age_class'),('DNST_NM','density'),
                          ('HEIGHT','height_code'),('HEIGHT_NM','height_class'),('갱신년도','update_year')]:
        result[prefix+target] = attrs.get(source,'')
    bounds = height_bounds(attrs.get('HEIGHT',''))
    result[prefix+'height_midpoint_m'] = int(attrs['HEIGHT']) if bounds else ''
    result[prefix+'height_lower_m'] = bounds[0] if bounds else ''
    result[prefix+'height_upper_m'] = bounds[1] if bounds and math.isfinite(bounds[1]) else ''
    return result, inside

def write_code_mapping():
    rows=[]
    page_map={'FRTP_CD':'4','KOFTR_GROU':'4-6','DMCLS_CD':'6-7','AGCLS_CD':'7','DNST_CD':'7-8'}
    for field, label, mapping in MAPS:
        for code, name in mapping.items():
            lower=upper=''
            if field=='DMCLS_CD':
                lower={'0':0,'1':6,'2':18,'3':30}[code]
                upper={'0':6,'1':18,'2':30,'3':''}[code]
            rows.append(dict(field=field,code=code,label=name,lower_inclusive=lower,upper_exclusive=upper,
                             definition_pdf_pages=page_map[field]))
    for h in range(0,41,2):
        lo,hi=height_bounds(str(h))
        rows.append(dict(field='HEIGHT',code=f'{h:02}',label=f'{lo:g}m 이상 '+(f'{hi:g}m 미만' if math.isfinite(hi) else ''),
                         lower_inclusive=lo,upper_exclusive=hi if math.isfinite(hi) else '',definition_pdf_pages='8-9'))
    save_csv(OUT/'forest_code_mapping.csv',rows)
