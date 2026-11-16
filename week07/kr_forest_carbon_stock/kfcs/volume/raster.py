"""Exact pixel/polygon intersections; physical areas measured in working CRS."""
from pathlib import Path
import math
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from pyproj import Transformer
from shapely.geometry import Polygon
from shapely.ops import transform
from .base import VolumeProvider, file_record

class RasterProvider(VolumeProvider):
    def __init__(self, config):
        self.config=config; self.v=config['volume']; self.path=Path(self.v['raster_path'])
        self.audit={}

    def provide(self, polygons):
        out=[]; scale=self.config['units']['m2_per_ha']
        working=polygons.to_crs(self.config['spatial']['crs'])
        with rasterio.open(self.path) as src:
            if src.crs is None: raise ValueError('Raster CRS missing')
            to_raster=Transformer.from_crs(working.crs,src.crs,always_xy=True).transform
            to_working=Transformer.from_crs(src.crs,working.crs,always_xy=True).transform
            inv=~src.transform
            for (_,row),geom in zip(polygons.iterrows(),working.geometry):
                result=np.nan; valid_area=0.; weighted=0.
                if not row.eligible:
                    out.append((0.,f'raster:{self.v["edition"]}','직접',0.)); continue
                if not geom.is_empty:
                    rg=transform(to_raster,geom)
                    x0,y0,x1,y1=rg.bounds
                    ij=[inv*(x,y) for x in (x0,x1) for y in (y0,y1)]
                    c0=max(0,math.floor(min(p[0] for p in ij))); c1=min(src.width,math.ceil(max(p[0] for p in ij)))
                    r0=max(0,math.floor(min(p[1] for p in ij))); r1=min(src.height,math.ceil(max(p[1] for p in ij)))
                    if c1>c0 and r1>r0:
                        data=src.read(self.v.get('band',1),window=Window(c0,r0,c1-c0,r1-r0),masked=True)
                        for rr,cc in np.ndindex(data.shape):
                            value=data[rr,cc]
                            if np.ma.is_masked(value) or not np.isfinite(value): continue
                            if value<0: raise ValueError('Negative valid volume pixel; set proper nodata')
                            r,c=r0+rr,c0+cc
                            pixel=Polygon([src.transform*(c,r),src.transform*(c+1,r),src.transform*(c+1,r+1),src.transform*(c,r+1)])
                            # Densify edges before a non-linear CRS transform.
                            step=max(pixel.length/self.v['edge_segments'],np.finfo(float).eps)
                            pixel=transform(to_working,pixel.segmentize(step))
                            overlap=geom.intersection(pixel).area
                            if overlap<=0: continue
                            density=float(value) if self.v['raster_unit']=='m3_per_ha' else float(value)*scale/pixel.area
                            weighted+=density*overlap; valid_area+=overlap
                if valid_area>0: result=weighted/valid_area
                out.append((result,f'raster:{self.v["edition"]}','직접' if np.isfinite(result) else '결측',valid_area/scale))
        result=pd.DataFrame(out,index=polygons.index,columns=['m3_per_ha','volume_source','volume_flag','raster_valid_area_ha'])
        self.audit={'missing_count':int(result.volume_flag.eq('결측').sum()),'missing_area_ha':float(polygons.loc[result.volume_flag.eq('결측'),'area_ha'].sum()),
                    'valid_pixel_overlap_ha':float(result.raster_valid_area_ha.sum()),
                    'partial_coverage_area_ha':float((polygons.area_ha-result.raster_valid_area_ha).clip(lower=0)[polygons.eligible].sum())}
        return result

    def provenance(self):
        return {'provider':'raster','input':file_record(self.path,edition=self.v['edition']), 'rules':self.v,'audit':self.audit,
                'area_crs':self.config['spatial']['crs'],'method':'exact pixel intersection / valid intersecting area; no imputation'}
