"""Independent checks from saved CSVs and source shape records (no pipeline import)."""
from pathlib import Path
from collections import defaultdict
import csv, math, json
import numpy as np
import shapefile, shapely
from shapely.geometry import shape
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def run():
    checks=[]
    def test(name,ok,detail):
        checks.append(dict(test=name,passed=bool(ok),detail=str(detail)))
        if not ok:raise AssertionError((name,detail))
    def close(a,b):return math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-9)
    trees=read(OUT/'carbon_A_trees.csv'); plots=read(OUT/'carbon_A_plots.csv')
    test('A 제외수·계산수',len(trees)==1750 and sum(r['included']=='True' for r in trees)==1694
         and sum(float(r['흉고직경_cm'])<6 for r in trees)==47 and sum(r['수종']=='등나무' for r in trees)==9,'1750−47−9=1694')
    included=[r for r in trees if r['included']=='True']
    eq={(r['species'],r['form'],r['part']):r for r in read(ROOT/'data/reference/nifos2014_allometry.csv')}
    pine=[r for r in included if r['수종']=='소나무']
    errors=[]
    for r in pine:
        d=float(r['흉고직경_cm']); target='강원지방소나무' if d>40 else '중부지방소나무'
        # Independently reconstruct a directly matched species, both forms.
        for form,col in [('a*D^b','carbon_kg'),('a*D^b*H^c','carbon_DH_kg')]:
            total=0
            for part in ('줄기','가지','잎','뿌리'):
                c=eq[target,form,part]; mass=float(c['a'])*d**float(c['b'])
                if form.endswith('H^c'):mass*=float(r['수고_m'])**float(c['c'])
                total+=mass
            errors.append(abs(total*.5-float(r[col])))
    test('소나무 모든나무 두식 독립재계산',max(errors)<1e-8,(len(pine),max(errors)))
    test('40cm 초과 소나무 전환',sum(float(r['흉고직경_cm'])>40 for r in pine)==7,'7그루')
    test('A 50구 합계 역산',all(close(sum(float(t['carbon_kg']) for t in included if t['plot_id']==p['plot_id'])/40,p['carbon_tC_ha']) for p in plots),'kg/(1000×0.04)=kg/40')
    bsum=read(OUT/'carbon_B_summary.csv')
    for year in (2013,2025):
        bs=read(OUT/f'carbon_B_polygons_{year}.csv')
        test(f'B{year} 모든 폴리곤 재적·AGB·뿌리·탄소',all(
            close(float(r['volume_m3_ha'])*float(r['WD'])*float(r['BEF']),r['AGB_t_ha']) and
            close(float(r['AGB_t_ha'])*float(r['R']),r['root_t_ha']) and
            close((float(r['AGB_t_ha'])+float(r['root_t_ha']))/2,r['carbon_tC_ha']) and
            close(float(r['carbon_total_tC'])*44/12,r['CO2_total_tCO2']) for r in bs),len(bs))
        shp=shapefile.Reader(str(ROOT/f'data/forest_map/merged/forest{year}_geumjeong.shp'),encoding='cp949')
        native={str(sr.record.as_dict()['FID_U']):shape(sr.shape.__geo_interface__).area/10000 for sr in shp.iterShapeRecords()}
        test(f'B{year} SHP 원도형면적 역대조',all(close(native[r['polygon_id']],r['area_ha']) for r in bs),len(native))
        for scope,tag in [('금정산','geumjeong'),('부산대 부산캠퍼스','campus'),('클리핑 전체','bbox')]:
            total=next(r for r in bsum if r['year']==str(year) and r['scope']==scope and r['group_dimension']=='total')
            test(f'B{year} {scope} 폴리곤→요약',close(sum(float(r[tag+'_carbon_tC']) for r in bs),total['carbon_total_tC']),'tC합')
    records=[]
    reader=shapefile.Reader(str(ROOT/'data/forest_map/merged/forest2025_geumjeong.shp'),encoding='cp949')
    for sr in reader.iterShapeRecords():
        a=sr.record.as_dict()
        if str(a['FRTP_CD']) in ('1','2','3'):records.append((str(a['FID_U']),shape(sr.shape.__geo_interface__)))
    foresttree=shapely.STRtree(np.array([g for fid,g in records],dtype=object))
    raw=read(OUT/'gedi_v3_quality_points.csv'); q=read(OUT/'carbon_C_points_in_forest.csv')
    to_xy=Transformer.from_crs(4326,5179,always_xy=True)
    x,y=to_xy.transform(np.array([float(r['lon']) for r in raw]),np.array([float(r['lat']) for r in raw]))
    points=shapely.points(x,y); hits=foresttree.query(points,predicate='intersects')
    valid_ids=set(hits[0]); expected={raw[i]['shot_number'] for i in valid_ids}
    test('C 원자료+SHP 독립 숲필터',expected=={r['shot_number'] for r in q},len(expected))
    membership=read(OUT/'carbon_boundary_membership.csv'); selected={r['polygon_id'] for r in membership if r['selected_geumjeong_20m']=='True'}
    test('타 산 제외 경계조건',not any(r['selected_geumjeong_20m']=='True' and r['excluded_other_mountain']=='True' for r in membership),'타 산 사전제외 성분0개 포함')
    csum=read(OUT/'carbon_C_summary.csv')
    for scope in ('금정산','부산대 부산캠퍼스','클리핑 전체'):
        subset=[r for r in q if (r['polygon_id'] in selected if scope=='금정산' else
               (129.0738<=float(r['lon'])<=129.0847 and 35.2299<=float(r['lat'])<=35.2386 if scope=='부산대 부산캠퍼스' else True))]
        r=next(r for r in csum if r['scope']==scope and r['forest_type']=='전체')
        mean=sum(float(r['agbd']) for r in subset)/len(subset)
        test('C '+scope+' 점평균/개수',len(subset)==int(r['n_points']) and close(mean,r['AGB_t_ha']),(len(subset),mean))
    old=read(OUT/'seo2011_vs_forestmap.csv'); new={r['plot_id']:r for r in read(OUT/'carbon_compare_plots.csv')}
    for year in (2013,2025):
        pairs=[(r['fm'+str(year)+'_polygon_id'] if r['fm'+str(year)+'_in_stocked_forest']=='True' else '',new[r['plot_id']].get(f'B{year}_polygon_id','')) for r in old]
        test(f'기존 공간중첩과{year} ID 일치',all(a==b for a,b in pairs),len(pairs))
    summary=read(OUT/'carbon_compare_methods.csv')
    test('6개 주요 결과 단위',all(close(float(r['carbon_tC_ha'])*float(r['forest_area_ha']),r['carbon_total_tC']) and
         close(float(r['carbon_total_tC'])*44/12,r['CO2_total_tCO2']) for r in summary),'tC/ha×ha=tC;44/12')
    test('입력 야장캠퍼스 내부0',not any(129.0738<=float(p['경도'])<=129.0847 and 35.2299<=float(p['위도'])<=35.2386 for p in plots),'캠퍼스A는실측표본아님')
    out=OUT/'carbon_independent_validation.csv'
    with out.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['test','passed','detail']);writer.writeheader();writer.writerows(checks)
    print(json.dumps({'independent_checks':len(checks),'all_passed':all(r['passed'] for r in checks)},ensure_ascii=False))
    return checks
if __name__=='__main__':run()
