"""Join supplied EGIS polygons to frozen GEDI points; no training or filtering."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import struct

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import rasterio
from rasterio.crs import CRS
from rasterio.warp import transform
import shapely
from shapely.geometry import Polygon, MultiPolygon, Point, box
from shapely.strtree import STRtree

FOREST={'310','320','330'}
LABELS={'110':'주거지역','120':'공업지역','130':'상업지역','140':'문화·체육·휴양지역',
 '150':'교통지역','160':'공공시설지역','210':'논','220':'밭','230':'시설재배지',
 '240':'과수원','250':'기타재배지','310':'활엽수림','320':'침엽수림','330':'혼효림',
 '410':'자연초지','420':'인공초지','510':'내륙습지','520':'연안습지','610':'자연나지',
 '620':'인공나지','710':'내륙수','720':'해양수'}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def md(df):
    def fmt(v):
        return f'{v:.3f}' if isinstance(v,(float,np.floating)) else str(v)
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join(
        '| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))


def read_dbf(path):
    data=path.read_bytes()
    n=struct.unpack_from('<I',data,4)[0]
    header,length=struct.unpack_from('<HH',data,8)
    fields=[]
    offset=1
    for at in range(32,header-1,32):
        f=data[at:at+32]
        fields.append((f[:11].split(b'\0')[0].decode('ascii'),offset,f[16]))
        offset+=f[16]
    rows=[]
    for i in range(n):
        row=data[header+i*length:header+(i+1)*length]
        assert row[0]==32
        rows.append({name:row[start:start+size].rstrip(b' \0').decode('cp949') for name,start,size in fields})
    return rows


def polygon_from_record(record):
    kind=struct.unpack_from('<I',record,0)[0]
    assert kind==5
    parts,n=struct.unpack_from('<II',record,36)
    starts=struct.unpack_from('<'+'I'*parts,record,44)
    xy=np.frombuffer(record,dtype='<f8',count=2*n,offset=44+4*parts).reshape(-1,2)
    rings=[xy[start:end] for start,end in zip(starts,list(starts[1:])+[n])]
    shells,holes=[],[]
    for ring in rings:
        signed=np.sum(ring[:-1,0]*ring[1:,1]-ring[1:,0]*ring[:-1,1])
        (shells if signed<0 else holes).append(ring)
    if not shells:
        return None,'no_clockwise_shell'
    # SHP convention: clockwise exterior rings, counterclockwise holes.
    shell_polys=[Polygon(ring) for ring in shells]
    assigned=[[] for _ in shells]
    for ring in holes:
        probe=Polygon(ring).representative_point()
        containers=[i for i,p in enumerate(shell_polys) if p.covers(probe)]
        if not containers:
            return None,'unassigned_hole'
        chosen=min(containers,key=lambda i:shell_polys[i].area)
        assigned[chosen].append(ring)
    polys=[Polygon(shell,h) for shell,h in zip(shells,assigned)]
    geom=polys[0] if len(polys)==1 else MultiPolygon(polys)
    return geom,str(shapely.is_valid_reason(geom))


def metric(g):
    e=g.agbd-g.prediction
    den=((g.agbd-g.agbd.mean())**2).sum()
    return dict(n=len(g),rmse=float(np.sqrt((e**2).mean())),mae=float(e.abs().mean()),
                r2=float(1-(e**2).sum()/den) if len(g)>1 and den>0 else None,
                mean_residual=float(e.mean()),mean_agbd=float(g.agbd.mean()))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=HERE/'outputs')
    out=parser.parse_args().output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use a new empty --output-dir.')
    provenance=json.loads((HERE/'provenance.json').read_text(encoding='utf-8'))
    paths=[HERE/'provenance.json',ROOT/'week07/step3/outputs/evaluation_input.csv',
        ROOT/'week07/step3/outputs/predictions.csv',ROOT/'week07/step1/outputs/top20_absolute_errors.csv']
    paths += [ROOT/f['path'] for t in provenance['tiles'] for f in t['files']]
    hashes={str(p.relative_to(ROOT)):sha(p) for p in paths}
    for t in provenance['tiles']:
        assert all(sha(ROOT/f['path'])==f['sha256'] for f in t['files'])
    data=pd.read_csv(paths[1],dtype={'shot_number':str})
    top=pd.read_csv(paths[3])[['point_id','error_rank']]
    data=data.rename(columns={'forest_status':'previous_forest_status'})
    data=data.merge(top,on='point_id',how='left',validate='one_to_one')
    assert len(data)==441 and data.point_id.is_unique
    geoms,attrs,bad_boxes,audit=[],[],[],[]
    crs=None
    for t in provenance['tiles']:
        stem=HERE/'source'/t['tile']/f'lv2_{t["tile"]}'
        tile_crs=CRS.from_wkt(stem.with_suffix('.prj').read_text(encoding='utf-8'))
        assert tile_crs.to_epsg()==5186 and (crs is None or tile_crs==crs)
        crs=tile_crs
        records=read_dbf(stem.with_suffix('.dbf'))
        raw=stem.with_suffix('.shp').read_bytes()
        at=100
        for i,row in enumerate(records):
            number,length=struct.unpack_from('>II',raw,at)
            assert number==i+1
            record=raw[at+8:at+8+length*2]
            at+=8+length*2
            code=row['L2_CODE']
            assert code in LABELS and row['L2_NAME']==LABELS[code]
            geom,reason=polygon_from_record(record)
            valid=geom is not None and geom.is_valid and not geom.is_empty
            audit.append(dict(tile=t['tile'],record_id=number,L2_CODE=code,L2_NAME=row['L2_NAME'],
                              valid=valid,reason=reason,IMG_DATE=row['IMG_DATE'],IMG_NAME=row['IMG_NAME']))
            if valid:
                geoms.append(geom)
                attrs.append(dict(tile=t['tile'],record_id=number,**row))
            else:
                bad_boxes.append(box(*struct.unpack_from('<dddd',record,4)))
        assert at==len(raw)
        print(f'Validated polygon topology and attributes: {t["tile"]}',flush=True)
    tree=STRtree(geoms)
    bad_tree=STRtree(bad_boxes) if bad_boxes else None
    xx,yy=map(np.asarray,transform('EPSG:4326',crs,data.longitude.tolist(),data.latitude.tolist()))
    lon,lat=transform(crs,'EPSG:4326',xx.tolist(),yy.tolist())
    assert np.allclose(lon,data.longitude,atol=1e-8,rtol=0) and np.allclose(lat,data.latitude,atol=1e-8,rtol=0)
    data['landcover_x_m'],data['landcover_y_m']=xx,yy
    point_rows,matches=[],[]
    for row,x,y in zip(data.itertuples(index=False),xx,yy):
        pt=Point(x,y)
        hits=tree.query(pt,predicate='intersects').tolist()
        nearby_bad=bool(bad_tree is not None and len(bad_tree.query(pt,predicate='intersects')))
        boundary=False
        for idx in hits:
            a=attrs[idx]
            dist=pt.distance(geoms[idx].boundary)
            boundary |= dist<=.001
            matches.append(dict(point_id=row.point_id,polygon_index=idx,tile=a['tile'],record_id=a['record_id'],
                L2_CODE=a['L2_CODE'],L2_NAME=a['L2_NAME'],IMG_DATE=a['IMG_DATE'],boundary_distance_m=dist))
        status='unclassified'
        reason='no_polygon_match'
        selected=None
        if nearby_bad:
            reason='invalid_geometry_bbox_at_point'
        elif len(hits)>1:
            reason='multiple_polygon_matches'
        elif boundary:
            reason='within_1mm_polygon_boundary'
        elif len(hits)==1:
            selected=attrs[hits[0]]
            status='forest' if selected['L2_CODE'] in FOREST else 'nonforest'
            reason='single_valid_polygon_interior'
        # A 12.5m-radius circle is a diagnostic approximation; it is not an
        # observed, energy-weighted GEDI footprint or a claim about geolocation.
        disk=pt.buffer(12.5,quad_segs=32)
        local=tree.query(disk,predicate='intersects').tolist()
        forest_parts=[geoms[i].intersection(disk) for i in local if attrs[i]['L2_CODE'] in FOREST]
        other_parts=[geoms[i].intersection(disk) for i in local if attrs[i]['L2_CODE'] not in FOREST]
        f=shapely.union_all(forest_parts)
        n=shapely.union_all(other_parts)
        union=shapely.union_all([f,n])
        overlap=f.intersection(n).area/disk.area
        unknown=max(0.,1-union.area/disk.area)
        disk_bad=bool(bad_tree is not None and len(bad_tree.query(disk,predicate='intersects')))
        point_rows.append(dict(point_id=row.point_id,landcover_status=status,classification_reason=reason,
            center_match_n=len(hits),L2_CODE=selected['L2_CODE'] if selected else '',
            L2_NAME=selected['L2_NAME'] if selected else '',landcover_tile=selected['tile'] if selected else '',
            IMG_DATE=selected['IMG_DATE'] if selected else '',IMG_NAME=selected['IMG_NAME'] if selected else '',
            polygon_boundary_distance_m=float(pt.distance(geoms[hits[0]].boundary)) if len(hits)==1 else np.nan,
            circle_r12p5_forest_fraction=float(f.area/disk.area),circle_r12p5_nonforest_fraction=float(n.area/disk.area),
            circle_r12p5_uncovered_fraction=unknown,circle_r12p5_conflicting_overlap_fraction=float(overlap),
            circle_r12p5_invalid_geometry_bbox_nearby=disk_bad,
            circle_r12p5_mixed_forest_nonforest=bool(f.area/disk.area>1e-6 and n.area/disk.area>1e-6)))
    result=data.merge(pd.DataFrame(point_rows),on='point_id',validate='one_to_one')
    assert len(result)==441
    assert np.allclose(result.circle_r12p5_forest_fraction+result.circle_r12p5_nonforest_fraction-
        result.circle_r12p5_conflicting_overlap_fraction+result.circle_r12p5_uncovered_fraction,1,atol=1e-8)
    pred=pd.read_csv(paths[2])
    joined=pred.merge(result[['point_id','landcover_status','L2_CODE','error_rank']],on='point_id',validate='many_to_one')
    summaries=[]
    for scope,g in [('all441',result),('agbd_ge100',result[result.agbd.ge(100)]),('previous_top20',result[result.error_rank.notna()])]:
        for status in ['forest','nonforest','unclassified']:
            h=g[g.landcover_status.eq(status)]
            summaries.append(dict(scope=scope,status=status,n=len(h),mean_agbd=float(h.agbd.mean()) if len(h) else None,
                circle_mixed_n=int(h.circle_r12p5_mixed_forest_nonforest.sum()),
                mean_forest_circle_fraction=float(h.circle_r12p5_forest_fraction.mean()) if len(h) else None))
    summary=pd.DataFrame(summaries)
    codes=result.groupby(['landcover_status','L2_CODE','L2_NAME'],dropna=False).agg(n=('point_id','size'),
        high_agbd_n=('agbd',lambda x:int(x.ge(100).sum())),top20_n=('error_rank',lambda x:int(x.notna().sum())),
        mean_agbd=('agbd','mean')).reset_index()
    submetrics=pd.DataFrame([dict(feature_set=k,status=s,**metric(g)) for (k,s),g in joined.groupby(['feature_set','landcover_status'])])
    out.mkdir(parents=True,exist_ok=True)
    for name,frame in [('landcover_points',result),('point_polygon_matches',pd.DataFrame(matches)),
        ('classification_summary',summary),('landcover_code_counts',codes),('existing_oof_by_landcover',submetrics),
        ('geometry_audit',pd.DataFrame(audit))]:
        frame.to_csv(out/f'{name}.csv',index=False,encoding='utf-8-sig')
    counts=result.landcover_status.value_counts().to_dict()
    counts={k:int(counts.get(k,0)) for k in ['forest','nonforest','unclassified']}
    invalid=sum(not a['valid'] for a in audit)
    high_summary=summary[summary.scope.eq('agbd_ge100')]
    top_summary=summary[summary.scope.eq('previous_top20')]
    report=f'''# EGIS 중분류 토지피복과 GEDI 441점 결합

## 결과

6개 도엽의 원본 좌표계 EPSG:5186으로 441점 좌표를 변환하고 실제 polygon과 결합했다. 중심점 기준 **산림 {counts['forest']}점, 비산림 {counts['nonforest']}점, 판정 불가 {counts['unclassified']}점**이다. 여기서 산림은 토지피복 코드 310·320·330에 포함되는 점이며 법적 산림 범위 또는 독립 현장 검증을 의미하지 않는다.

{md(summary)}

{md(codes)}

식생 후보 중 257/441점(58.3%)이 토지피복 비산림으로 분류됐다. 비산림에는 인공초지 90점, 밭 48점, 논 44점 등이 포함된다. SCL 식생 분류와 산림 분류를 구분할 필요가 실제 자료에서 확인됐다. 고 AGBD 56점 중 산림은 51점, 비산림은 5점이며, 기존 오차 상위 20점 중 산림은 16점이다. 따라서 비산림 혼입을 정리하는 것만으로 높은 AGBD 과소 예측이 해결된다고 볼 수 없다. 산림 중심점 184점 중 22점은 반경 12.5m 근사 원에 비산림 면적도 포함돼 주변 혼합 플래그가 남는다.

## 고 AGBD 및 큰 오차점

고 AGBD(100 Mg/ha 이상) 56점과 기존 개선 RF 절대오차 상위 20점은 고정한 동일 ID로 살폈다. 피복 분류에 평가 AGBD나 예측 오차를 사용하지 않았다.

{md(result[result.error_rank.notna()].sort_values('error_rank')[['error_rank','agbd','landcover_status','L2_NAME','circle_r12p5_forest_fraction','circle_r12p5_nonforest_fraction']])}

## 기존 예측의 피복별 진단

{md(submetrics)}

위 값은 **전체 식생 후보를 학습한 기존 모델의 저장 OOF 예측을 피복별로 나눈 기술적 지표**다. 산림 표본만으로 재학습한 모델의 성능이 아니다. 피복별 AGBD 분포와 표본 수가 달라 지표를 단순 비교해 피복 효과라고 해석하지 않는다. 기존 모델 선택에 사용한 CV를 다시 진단한 것이며 최종 외부 성능으로 보고하지 않는다.

## 결합 규칙과 주변 피복

DBF `L2_CODE`를 사용하고 CP949로 읽은 `L2_NAME`이 중분류 코드명과 일치하는지 전체 레코드에서 확인했다. 원본 파일·CPG는 수정하지 않았다. SHP 링 방향에 따라 외곽·홀을 구성하고 Shapely {shapely.__version__}로 유효성을 확인했다. 총 {len(audit)}개 도형 중 유효하지 않거나 구성할 수 없는 도형은 {invalid}개다. 도형을 임의 수정하지 않고, 해당 bbox에 걸리는 점은 판정 불가로 처리한다.

중심점이 단일 유효 도형 내부에 있을 때만 코드로 판정한다. 다중 도형 중첩, 경계에서 1mm 이내, 무자료는 별도 사유로 남긴다. 도형 경계까지 거리는 동일 피복 사이 경계나 도엽 경계도 포함하므로 산림 경계 거리로 해석하지 않는다.

반경 12.5m 원형 근사 영역을 실제 polygon과 교차하여 산림·비산림 면적 비율, 무자료·서로 다른 피복의 중첩 비율을 계산했다. 원은 128개 선분으로 근사했다. 중첩은 union으로 중복 면적을 제거하며 산림/비산림 상충 면적은 따로 보존한다. 가까운 무효 도형 bbox도 표시했다. 이는 실제 GEDI footprint 모양·에너지 가중치·위치 오차를 복원한 값이 아니다. 중심점 산림 판정이 관측 영역 전체 산림임을 보장하지 않는다. 혼합 플래그는 양쪽 면적 비율이 각각 1e-6보다 클 때이며 필터링에 사용하지 않았다.

## 자료 시점과 한계

사용자가 신청 화면에서 선택한 2025년 중분류 자료다. 원본 `IMG_DATE`의 값은 geometry_audit.csv와 개별점 CSV에 그대로 남겼다. `20241231` 같은 연말 일자는 메타데이터 필드값이며 실제 각 도형의 정확한 촬영일이라고 단정하지 않는다. 2025-04-25 GEDI 관측 당시 피복과의 일치를 독립 검증하지 않았다. 중분류 일반화와 경계·위치 오차를 고려해야 한다.

## 지도

![피복별 표본 위치와 구성](landcover_overview.png)

## 다음 분석의 기준

이번 단계에서는 441점을 유지하고 피복 판정 열만 추가했다. 산림 표본만의 평가를 진행하려면 기존 분할을 유지할지와 학습 표본 수를 먼저 고정하고, 산림/비산림·주변 혼합 조건별 표본 수를 보고해야 한다. 현재 저장 예측의 부분집합 점수를 산림 전용 학습 결과로 바꿔 부르지 않는다.

## 재현

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\landcover\\classify_points.py --output-dir .\\week07\\landcover\\reproduction_outputs
```

기존 Python 환경에 Shapely {shapely.__version__}를 추가했다. 원본 토지피복 파일의 보존 해시와 이번 분석 입력 해시를 확인했으며, 판정과 저장 결과를 validation.json에 기록했다. 모델 재학습·추가 다운로드·원본 수정은 하지 않았다.
'''
    (out/'report.md').write_text(report,encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(12,5.5),layout='constrained')
    colors={'forest':'#228b22','nonforest':'#d99432','unclassified':'#777777'}
    for status,col in colors.items():
        h=result[result.landcover_status.eq(status)]
        axes[0].scatter(h.landcover_x_m/1000,h.landcover_y_m/1000,s=12,c=col,label=f'{status} ({len(h)})')
    h=result[result.error_rank.notna()]
    axes[0].scatter(h.landcover_x_m/1000,h.landcover_y_m/1000,s=65,facecolors='none',edgecolors='black',linewidth=.8,label='Previous top20 errors')
    axes[0].set(xlabel='EPSG:5186 easting (km)',ylabel='Northing (km)',title='Center-point landcover classification')
    axes[0].set_aspect('equal')
    axes[0].legend(fontsize=8)
    bottom=np.zeros(3)
    scopes=['all441','agbd_ge100','previous_top20']
    for status,col in colors.items():
        vals=[int(summary[(summary.scope==s)&(summary.status==status)].n.iloc[0]) for s in scopes]
        axes[1].bar(scopes,vals,bottom=bottom,color=col,label=status)
        for i,v in enumerate(vals):
            if v:
                axes[1].text(i,bottom[i]+v/2,str(v),ha='center',va='center',fontsize=9)
        bottom+=np.array(vals)
    axes[1].set(ylabel='Point count',title='Composition of fixed diagnostic groups')
    axes[1].legend()
    fig.savefig(out/'landcover_overview.png',dpi=170)
    plt.close(fig)
    saved=pd.read_csv(out/'landcover_points.csv')
    assert len(saved)==441 and saved.point_id.is_unique and set(saved.point_id)==set(data.point_id)
    assert all(sha(ROOT/p)==h for p,h in hashes.items())
    validation=dict(status='passed',n=441,counts=counts,source_epsg=5186,source_geometry_n=len(audit),
        invalid_geometry_n=invalid,unique_point_ids=True,coordinate_roundtrip_passed=True,
        circle_area_accounting_passed=True,source_hashes_unchanged=True,input_sha256=hashes,
        model_refit=False,previous_predictions_unchanged=True,shapely_version=shapely.__version__,
        raw_encoding='CP949 strict decode; code/name pairs checked',script_sha256=sha(Path(__file__)))
    (out/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(counts))
    print(summary.to_string(index=False))


if __name__=='__main__':
    main()
