"""Independently check center matches with raw-ring even/odd ray casting."""
from pathlib import Path
import json
import struct
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent


def inside_rings(x,y,rings):
    inside=np.zeros(len(x),dtype=bool)
    for ring in rings:
        ax,ay=ring[:-1].T
        bx,by=ring[1:].T
        nonhorizontal=by!=ay
        ax,ay,bx,by=[v[nonhorizontal] for v in [ax,ay,bx,by]]
        crossing=(ay[None]>y[:,None])!=(by[None]>y[:,None])
        intersection=ax[None]+(y[:,None]-ay[None])*(bx-ax)[None]/(by-ay)[None]
        inside ^= np.logical_xor.reduce(crossing & (x[:,None]<intersection),axis=1)
    return inside


def main():
    # Exterior, hole, and disjoint exterior; independent of ring orientation.
    outer=np.array([[0,0],[10,0],[10,10],[0,10],[0,0]],float)
    hole=np.array([[3,3],[7,3],[7,7],[3,7],[3,3]],float)
    island=outer+20
    assert inside_rings(np.array([1,5,21,15]),np.array([1,5,21,15]),[outer,hole,island]).tolist()==[True,False,True,False]
    out=HERE/'outputs'
    data=pd.read_csv(out/'landcover_points.csv')
    matches=pd.read_csv(out/'point_polygon_matches.csv')
    x,y=data.landcover_x_m.to_numpy(),data.landcover_y_m.to_numpy()
    expected=set()
    n_records=0
    for shp in sorted((HERE/'source').rglob('*.shp')):
        tile=shp.stem.split('_')[1]
        raw=shp.read_bytes()
        at=100
        while at<len(raw):
            record_id,size=struct.unpack_from('>II',raw,at)
            r=raw[at+8:at+8+size*2]
            at+=8+size*2
            n_records+=1
            lo,bot,hi,top=struct.unpack_from('<dddd',r,4)
            candidates=np.flatnonzero((x>=lo)&(x<=hi)&(y>=bot)&(y<=top))
            if not len(candidates):
                continue
            parts,n=struct.unpack_from('<II',r,36)
            starts=struct.unpack_from('<'+'I'*parts,r,44)
            xy=np.frombuffer(r,dtype='<f8',count=n*2,offset=44+parts*4).reshape(-1,2)
            rings=[xy[s:e] for s,e in zip(starts,list(starts[1:])+[n])]
            keep=inside_rings(x[candidates],y[candidates],rings)
            expected.update((data.point_id.iloc[i],tile,record_id) for i in candidates[keep])
    observed=set(zip(matches.point_id,matches.tile.astype(str),matches.record_id))
    assert expected==observed
    assert len(expected)==441 and matches.boundary_distance_m.gt(.001).all()
    forest=data.L2_CODE.astype(str).isin(['310','320','330'])
    assert data.landcover_status.eq('forest').equals(forest)
    assert data.point_id.is_unique and len(data)==441
    result=dict(status='passed',method='Raw SHP rings, even/odd ray casting; no Shapely geometry used',
                source_records_checked=n_records,matches=441,all_point_polygon_matches_identical=True,
                synthetic_hole_and_multipart_check=True,code_status_mapping_verified=True)
    (out/'independent_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
