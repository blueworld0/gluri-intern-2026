"""Extract and structurally verify six supplied EGIS SHPs, then optionally delete ZIPs.

Uses stdlib for SHP/SHX/DBF structure checks and rasterio for CRS. It does
not claim topological geometry validity or classify landcover from tile bounds.
"""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import zipfile
import zlib

import numpy as np
import pandas as pd
from rasterio.crs import CRS
from rasterio.warp import transform

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TILES = ['358161','358162','358163','358164','359131','359133']
EXTENSIONS = {'.shp','.shx','.dbf','.prj'}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dbf_info(data):
    n = struct.unpack_from('<I',data,4)[0]
    header, length = struct.unpack_from('<HH',data,8)
    assert data[header-1] == 13 and len(data) >= header+n*length
    fields, offset = [], 1
    for at in range(32,header-1,32):
        f = data[at:at+32]
        name = f[:11].split(b'\0')[0].decode('ascii')
        fields.append(dict(name=name,type=chr(f[11]),length=f[16],decimals=f[17],offset=offset))
        offset += f[16]
    assert offset == length
    # Codes and schema are ASCII. No guessed encoding is written as a .cpg.
    values = {f['name']:set() for f in fields}
    deleted = 0
    for i in range(n):
        row = data[header+i*length:header+(i+1)*length]
        assert row[0] in [32,42]
        deleted += int(row[0]==42)
        for f in fields:
            v=row[f['offset']:f['offset']+f['length']].strip()
            try:
                v=v.decode('ascii')
            except UnicodeDecodeError:
                v='<non-ASCII text; original bytes preserved>'
            if len(values[f['name']]) < 100:
                values[f['name']].add(v)
    return dict(record_n=n,deleted_n=deleted,fields=fields,
                distinct_value_sample={k:sorted(v) for k,v in values.items()},
                language_driver_id=data[29],text_encoding='not independently confirmed; original DBF preserved')


def shape_info(shp,shx):
    for data in [shp,shx]:
        assert struct.unpack_from('>I',data,0)[0]==9994
        assert struct.unpack_from('>I',data,24)[0]*2==len(data)
        assert struct.unpack_from('<I',data,28)[0]==1000
    shape_type=struct.unpack_from('<I',shp,32)[0]
    assert shape_type in [5,15,25] and shape_type==struct.unpack_from('<I',shx,32)[0]
    assert (len(shx)-100)%8==0
    n=(len(shx)-100)//8
    null_n,part_n,point_n=0,0,0
    expected=100
    for i in range(n):
        offset, words=struct.unpack_from('>II',shx,100+8*i)
        offset*=2
        assert offset==expected
        number,content_words=struct.unpack_from('>II',shp,offset)
        assert number==i+1 and words==content_words
        record=shp[offset+8:offset+8+words*2]
        assert len(record)==words*2
        expected=offset+8+words*2
        kind=struct.unpack_from('<I',record,0)[0]
        if kind==0:
            null_n+=1
            continue
        assert kind==shape_type
        parts,points=struct.unpack_from('<II',record,36)
        assert parts>0 and points>=4 and len(record)>=44+parts*4+points*16
        starts=struct.unpack_from('<'+'I'*parts,record,44)
        assert starts[0]==0 and all(a<b for a,b in zip(starts,starts[1:])) and starts[-1]<points
        xy=np.frombuffer(record,dtype='<f8',count=points*2,offset=44+parts*4).reshape(-1,2)
        assert np.isfinite(xy).all()
        for start,end in zip(starts,list(starts[1:])+[points]):
            assert end-start>=4 and np.array_equal(xy[start],xy[end-1])
        part_n+=parts
        point_n+=points
    assert expected==len(shp)
    return dict(record_n=n,null_geometry_n=null_n,part_n=part_n,vertex_n=point_n,
                shape_type=shape_type,bounds=list(struct.unpack_from('<dddd',shp,36)),
                checks='header lengths, every SHX offset, record count, finite XY and ring closure; not topology validation')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--delete-archives',action='store_true')
    args=parser.parse_args()
    workspace=ROOT.resolve()
    source=(HERE/'source').resolve()
    assert source.is_relative_to(workspace) and source.parent==HERE.resolve()
    if source.exists():
        raise RuntimeError('Source already exists; refusing overwrite or archive deletion.')
    archives=[workspace/'week07'/f'SG04_{tile}_20251113.zip' for tile in TILES]
    # Check all archives before any extraction or deletion.
    for tile,path in zip(TILES,archives):
        assert path.resolve().is_relative_to(workspace/'week07') and path.is_file()
        with zipfile.ZipFile(path) as z:
            assert z.testzip() is None
            expected={f'lv2_{tile}{ext}' for ext in EXTENSIONS}
            assert len(z.infolist())==4 and set(z.namelist())==expected
    manifest=[]
    source.mkdir()
    points=pd.read_csv(workspace/'week06/step5/outputs/vegetation_candidates.csv',usecols=['point_id','longitude','latitude'])
    covered=np.zeros(len(points),bool)
    for tile,path in zip(TILES,archives):
        tile_dir=source/tile
        tile_dir.mkdir()
        members=[]
        with zipfile.ZipFile(path) as z:
            for item in z.infolist():
                target=(tile_dir/item.filename).resolve()
                assert target.parent==tile_dir.resolve() and target.suffix in EXTENSIONS
                data=z.read(item)
                assert zlib.crc32(data)==item.CRC
                with target.open('xb') as f:
                    f.write(data)
                assert target.read_bytes()==data
                members.append(dict(path=str(target.relative_to(workspace)),bytes=len(data),
                                    zip_crc32=f'{item.CRC:08x}',sha256=sha(data)))
        stem=tile_dir/f'lv2_{tile}'
        shape=shape_info(stem.with_suffix('.shp').read_bytes(),stem.with_suffix('.shx').read_bytes())
        dbf=dbf_info(stem.with_suffix('.dbf').read_bytes())
        assert shape['record_n']==dbf['record_n'] and dbf['deleted_n']==0
        crs=CRS.from_wkt(stem.with_suffix('.prj').read_text(encoding='utf-8'))
        assert crs.is_projected
        x,y=map(np.asarray,transform('EPSG:4326',crs,points.longitude.tolist(),points.latitude.tolist()))
        left,bottom,right,top=shape['bounds']
        inside=(x>=left)&(x<=right)&(y>=bottom)&(y<=top)
        covered |= inside
        manifest.append(dict(tile=tile,archive=str(path.relative_to(workspace)),archive_bytes=path.stat().st_size,
            archive_sha256=sha(path.read_bytes()),files=members,shape=shape,dbf=dbf,
            crs_wkt=crs.to_wkt(),epsg=crs.to_epsg(),points_within_file_bbox_n=int(inside.sum())))
        print(json.dumps(dict(tile=tile,features=shape['record_n'],epsg=crs.to_epsg(),
                             fields=[f['name'] for f in dbf['fields']],code_samples=dbf['distinct_value_sample']),ensure_ascii=False),flush=True)
    assert covered.all()
    provenance=dict(status='passed',dataset='User-supplied EGIS 2025 medium landcover (lv2)',
        year_basis='2025 selected in user download UI; ZIP suffix is not assumed to be imagery acquisition date',
        source_file_n=24,tiles=manifest,all441_within_combined_file_bboxes=True,
        actual_point_polygon_classification='not yet performed',topological_geometry_validation='not performed',
        archives_deleted=[],kept='All original SHP/SHX/DBF/PRJ bytes; no member discarded')
    destination=HERE/'provenance.json'
    def write():
        destination.write_text(json.dumps(provenance,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    write()
    if args.delete_archives:
        # Exact six files only; no recursive or computed-directory deletion.
        # Recheck retained bytes AND each archive hash immediately before unlink.
        for entry,path in zip(manifest,archives):
            for member in entry['files']:
                assert sha((workspace/member['path']).read_bytes())==member['sha256']
            assert path.resolve().parent==(workspace/'week07').resolve()
            assert sha(path.read_bytes())==entry['archive_sha256']
        for entry,path in zip(manifest,archives):
            path.unlink()
            assert not path.exists()
            provenance['archives_deleted'].append(entry['archive'])
            write()
    (HERE/'README.md').write_text('''# EGIS 토지피복지도 준비

사용자가 받은 2025년 중분류 6개 도엽을 `source/<도엽번호>/`에 보존했다. 각 도엽의 `.shp`, `.shx`, `.dbf`, `.prj`는 모두 필요하므로 24개 파일을 남겼다. ZIP 외에 버린 내부 파일은 없다.

ZIP CRC, 추출 전후 바이트 일치, SHP/SHX 레코드·오프셋·좌표 유한성·링 폐합, DBF 레코드 수와 필드, PRJ 좌표계를 검사했다. 좌표계와 속성 목록·코드 샘플은 `provenance.json`에 기록했다. 도형 위상 유효성을 검사한 것은 아니다.

441점 전체가 6개 파일의 결합 bbox 안에 포함됨을 확인했다. 이는 실제 polygon 포함이나 산림 판정을 의미하지 않는다. 피복 분류와 경계 검토는 다음 분석에서 수행한다. 촬영·갱신 시기는 ZIP 이름만으로 확정하지 않았다. DBF 문자 인코딩은 확인 없이 새 CPG를 만들지 않고 원본을 보존했다.

사용자 요청에 따라 검증 완료 후 이번 ZIP 6개를 삭제했다. 원본 ZIP의 파일명·크기·SHA-256과 보존 파일의 SHA-256·CRC를 `provenance.json`에 남겼다. 기존 주차 데이터·분석 결과는 변경하지 않았다.
''',encoding='utf-8')
    print(f'Kept 24 required files; deleted {len(provenance["archives_deleted"])} archives; bbox coverage 441/441.')


if __name__=='__main__':
    main()
