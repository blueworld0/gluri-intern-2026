"""AOI selection and area in the configured metric working projection."""
import io
import zipfile
from pathlib import Path
import geopandas as gpd
import shapely

def polygonal(g):
    if g is None or g.is_empty: return shapely.Polygon()
    if g.geom_type in ('Polygon','MultiPolygon'): return g
    return shapely.union_all([polygonal(p) for p in g.geoms]) if hasattr(g,'geoms') else shapely.Polygon()

def read_boundary(config):
    b=config['boundary']
    if b['path'] in (None,'none'): return None
    path=str(b['path'])
    archive=path.removeprefix('zip://').split('!')[0] if path.startswith('zip://') else path
    if archive.lower().endswith('.zip'):
        from .volume.tables import zip_members
        members=list(zip_members(archive))
        selected=b.get('member') or (path.split('!',1)[1] if path.startswith('zip://') and '!' in path else None)
        shapes=[name for name,_ in members if name.lower().endswith('.shp')]
        if selected: shapes=[name for name in shapes if name==selected or Path(name).name==selected]
        if len(shapes)!=1: raise ValueError(f'Boundary ZIP requires one unambiguous SHP; use boundary.member: {shapes}')
        stem=str(Path(shapes[0]).with_suffix(''))
        # Flatten only this SHP's sidecars into an in-memory UTF-8 ZIP. GDAL cannot
        # reliably traverse CP949 archive directory names. The source stays read-only.
        memory=io.BytesIO()
        with zipfile.ZipFile(memory,'w') as z:
            for name,blob in members:
                if str(Path(name).with_suffix(''))==stem: z.writestr(Path(name).name,blob)
        g=gpd.read_file(memory.getvalue(),engine='pyogrio',encoding=b.get('encoding',config['inputs']['encoding']))
    else:
        options={'encoding':b['encoding']} if b.get('encoding') else {}
        g=gpd.read_file(path,engine='pyogrio',**options)
    if g.crs is None: raise ValueError('AOI CRS missing')
    if b.get('select'):
        for field,value in b['select'].items(): g=g[g[field]==value]
    if g.empty: raise ValueError('AOI selection returned no polygons')
    result=polygonal(shapely.union_all(shapely.make_valid(g.to_crs(config['spatial']['crs']).geometry.values)))
    if result.is_empty: raise ValueError('AOI has no polygonal area')
    return result

def clip(g, boundary, config):
    g=g.copy()
    if boundary is not None:
        g=g[g.intersects(boundary)].copy()
        g.geometry=g.geometry.intersection(boundary).map(polygonal)
    g['area_ha']=g.geometry.area/config['units']['m2_per_ha']
    if g.empty: raise ValueError('No forest-map records intersect the AOI')
    return g.reset_index(drop=True)
