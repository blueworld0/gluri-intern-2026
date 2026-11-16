"""Read raw SHP files. Code tables are a CSV snapshot of week 4, not imports."""
from pathlib import Path
import re
import pandas as pd
import geopandas as gpd
import shapely
from .aoi import polygonal

def code(value):
    if pd.isna(value): return ''
    s=str(value).strip()
    return s[:-2] if s.endswith('.0') else s

def read_forest_maps(config, boundary=None):
    ref=Path(config['reference_dir'])
    codes=pd.read_csv(ref/'forest_map_codes.csv',dtype=str).fillna('')
    sido=pd.read_csv(ref/'sido_codes.csv',dtype=str).set_index('input_code')
    frames=[]; audit={'files':[], 'overlap_removed_ha':0., 'fully_overlapped_count':0}
    for path in config['inputs']['forest_maps']:
        p=Path(path)
        if not re.fullmatch(r'\d{5}',p.stem): raise ValueError(f'SHP filename must be a 5-digit district code: {p}')
        prefix=p.stem[:2]
        if prefix not in sido.index: raise ValueError(f'Unknown sido code {prefix}')
        g=gpd.read_file(p,encoding=config['inputs']['encoding'],engine='pyogrio')
        if g.crs is None: raise ValueError(f'CRS missing: {p}')
        g=g.to_crs(config['spatial']['crs'])
        g['polygon_id']=[f'{p.stem}:{i}' for i in range(len(g))]
        g['source_row']=range(len(g)); g['district_code']=p.stem
        g['source_file']=str(p)
        g['sido_code']=sido.loc[prefix,'canonical_code']; g['sido_name']=sido.loc[prefix,'sido_name']
        g.geometry=shapely.make_valid(g.geometry.values)
        g.geometry=g.geometry.map(polygonal)
        if boundary is not None: g=g[g.intersects(boundary)].copy()
        for field, label in [('FRTP_CD','forest_type'),('KOFTR_GROU','species_group'),('DMCLS_CD','diameter_class'),('AGCLS_CD','age_label'),('DNST_CD','density_class')]:
            lookup=codes[codes.field==field].set_index('code').label.to_dict()
            if field not in g: g[field]=''
            g[field]=g[field].map(code)
            # Week 4 prefers nonempty source *_NM values. Unknown codes remain unclassified.
            known=g[field].map(lookup)
            source_label=codes[codes.field==field].label_field.iloc[0]
            if source_label in g:
                provided=g[source_label].fillna('').astype(str).str.strip()
                known=provided.where(provided.ne('') & known.notna(),known)
            g[label]=known.fillna('미분류')
        g['age_class']=pd.to_numeric(g['AGCLS_CD'],errors='coerce')
        g['eligible']=g.FRTP_CD.isin(config['classification']['eligible_codes'])
        g['classification_flag']=g[['forest_type','species_group','diameter_class','age_label','density_class']].eq('미분류').any(axis=1).map({True:'미분류',False:'해석됨'})
        audit['files'].append({'path':str(p),'retained_candidates':len(g)})
        frames.append(g)
    g=gpd.GeoDataFrame(pd.concat(frames,ignore_index=True),crs=config['spatial']['crs'])
    # Stable configured file order and source record order: first geometry wins.
    original=g.geometry.to_numpy(); tree=shapely.STRtree(original); accepted={}
    cleaned=[]
    for i,geom in enumerate(original):
        prior=[accepted[int(j)] for j in tree.query(geom,predicate='intersects') if int(j)<i and int(j) in accepted]
        result=polygonal(geom.difference(shapely.union_all(prior))) if prior else geom
        audit['overlap_removed_ha']+=(geom.area-result.area)/config['units']['m2_per_ha']
        audit['fully_overlapped_count']+=int(result.is_empty)
        accepted[i]=result; cleaned.append(result)
    g.geometry=cleaned  # Empty records retained, including their source IDs.
    return g,audit
