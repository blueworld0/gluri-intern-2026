"""Local forest-map clipping, field comparison, transitions and visit planning."""
from collections import Counter,defaultdict
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import shapely
from shapely.geometry import Point,LineString
from pyproj import Geod

from forest_map_common import *

def area_summary(layers):
    rows,coverage=[],[]
    for year,features in layers.items():
        union=shapely.union_all([f['geometry'] for f in features])
        for scope,region in [('금정산·캠퍼스 전체',CLIP),('부산대 캠퍼스',CAMPUS_GEOM)]:
            covered=union.intersection(region).area
            wooded=shapely.union_all([f['geometry'] for f in features if stocked(f['attrs'])]).intersection(region).area
            coverage.append(dict(year=year,scope=scope,region_area_ha=region.area/10000,
                                 mapped_area_ha=covered/10000,stocked_forest_area_ha=wooded/10000,
                                 unmapped_area_ha=(region.area-covered)/10000))
            for code_key,label_key,dimension in [('FRTP_CD','FRTP_NM','임상'),('KOFTR_GROU','KOFTR_NM','수종'),
                                                 ('DMCLS_CD','DMCLS_NM','경급'),('AGCLS_CD','AGCLS_NM','영급')]:
                totals=defaultdict(float); counts=Counter()
                for feature in features:
                    area=feature['geometry'].intersection(region).area
                    if area<=1e-6:
                        continue
                    attrs=feature['attrs']
                    key=(attrs[code_key],attrs[label_key] or '미분류/해당없음')
                    totals[key]+=area; counts[key]+=1
                for (code,label),area in sorted(totals.items()):
                    rows.append(dict(year=year,scope=scope,dimension=dimension,code=code,label=label,
                                     area_ha=area/10000,feature_parts=counts[(code,label)]))
    save_csv(OUT/'forest_area_summary.csv',rows)
    save_csv(OUT/'forest_coverage_summary.csv',coverage)
    return rows,coverage

CONIFERS={'소나무','곰솔','리기다소나무','잣나무','편백','편백나무','삼나무','낙엽송','노간주나무'}
OAKS={'상수리나무','신갈나무','굴참나무','떡갈나무','졸참나무','갈참나무'}
MEMBERS={'11':{'소나무'},'12':{'잣나무'},'13':{'낙엽송'},'14':{'리기다소나무'},'15':{'곰솔'},
         '17':{'편백','편백나무'},'18':{'삼나무'},'31':{'상수리나무'},'32':{'신갈나무'},
         '33':{'굴참나무'},'34':{'갈참나무','떡갈나무','졸참나무'},
         '35':{'오리나무','물오리나무','사방오리나무'},'36':{'고로쇠나무'},'39':{'밤나무'},
         '40':{'물푸레나무','들메나무'},'41':{'서어나무','개서어나무'},
         '42':{'때죽나무'},'46':{'벚나무','산벚나무'},'47':{'느티나무'},'48':{'층층나무'},
         '49':{'아까시나무'},'66':{'사스레피나무'}}

def species_comparison(name,code,in_forest):
    if not in_forest:
        return '비교불가(숲 폴리곤 밖)','비교불가'
    names=set(name.split(' / '))
    if code=='77':
        return '비교불가(혼효림·단일수종 아님)','혼효림(직접비교불가)'
    if code in ('30','60','10'):
        compatible=any((n in CONIFERS)==(code=='10') for n in names)
        return '비교불가(기타수종 묶음)', '상위군 호환' if compatible else '상위군 불일치'
    if code in MEMBERS:
        included=bool(names & MEMBERS[code])
        if code=='34':
            return '비교불가(기타참나무류 묶음)', '정의서 그룹 포함' if included else ('참나무 상위군 호환' if names & OAKS else '상위군 불일치')
        return ('일치' if included else '불일치'), ('정의서 수종군 포함' if included else '불일치')
    return '비교불가(대응표 미정)', '비교불가'

def compare_plots(indices,plots,summaries):
    lookup={r['plot_id']:r for r in summaries}
    rows=[]
    for plot in plots:
        stats=lookup[plot['plot_id']]
        point=Point(*TO_XY.transform(float(plot['경도']),float(plot['위도'])))
        dbh=float(stats['평균_흉고직경_cm']); height=float(stats['평균_수고_m'])
        mean_code=diameter_code(dbh)
        row=dict(plot_id=plot['plot_id'],zone=plot['zone'],lat=float(plot['위도']),lon=float(plot['경도']),
                 field_survey_period='2009-05/2010-09',field_dominant_species=stats['우점_수종'],
                 field_mean_dbh_cm=dbh,field_mean_height_m=height,field_mean_dbh_class=DIAMETER[mean_code],
                 field_stem_count=int(plot['나무 수']))
        for year,index in indices.items():
            prefix=f'fm{year}_'
            attrs,hit=attach_info(index,point,prefix)
            row.update(attrs)
            inside=attrs[prefix+'in_stocked_forest']
            code=attrs[prefix+'species_code']
            exact,group=species_comparison(stats['우점_수종'],code,inside)
            row[prefix+'species_exact_comparison']=exact
            row[prefix+'species_group_comparison']=group
            mapped=attrs[prefix+'diameter_code']
            row[prefix+'mean_dbh_comparison']=('일치' if mean_code==mapped else '불일치') if inside and mapped else '비교불가'
            if year==2025:
                bounds=height_bounds(attrs[prefix+'height_code']) if inside else None
                row[prefix+'mean_height_comparison']=('일치' if bounds[0]<=height<bounds[1] else '불일치') if bounds else '비교불가'
                row[prefix+'height_midpoint_minus_field_m']=attrs[prefix+'height_midpoint_m']-height if bounds else ''
            else:
                row[prefix+'mean_height_comparison']='비교불가(HEIGHT 미제공)'
        rows.append(row)
    save_csv(OUT/'seo2011_vs_forestmap.csv',rows)
    counts=[]
    for year in indices:
        for metric,key in [('정밀 수종','species_exact_comparison'),('수종 그룹','species_group_comparison'),
                           ('평균 흉고직경 경급','mean_dbh_comparison'),('평균 수고 등급','mean_height_comparison')]:
            values=Counter(r[f'fm{year}_{key}'] for r in rows)
            for status,count in sorted(values.items()):
                counts.append(dict(year=year,metric=metric,status=status,plots=count))
    save_csv(OUT/'forest_field_agreement_summary.csv',counts)
    return rows,counts

def transitions(old,new):
    totals=defaultdict(float)
    geoms=np.array([f['geometry'] for f in new],dtype=object)
    tree=shapely.STRtree(geoms)
    overlay_area=0.
    for i,feature in enumerate(old):
        hits=tree.query(feature['geometry'],predicate='intersects')
        if not len(hits):
            continue
        intersections=shapely.intersection(geoms[hits],feature['geometry'])
        for j,geom in zip(hits,intersections):
            if geom.area<=1e-6:
                continue
            a,b=feature['attrs'],new[int(j)]['attrs']
            key=(a['FRTP_CD'],a['KOFTR_GROU'],b['FRTP_CD'],b['KOFTR_GROU'])
            totals[('금정산·캠퍼스 전체',)+key]+=geom.area
            overlay_area+=geom.area
            campus_area=geom.intersection(CAMPUS_GEOM).area
            if campus_area>1e-6:
                totals[('부산대 캠퍼스',)+key]+=campus_area
        if (i+1)%500==0:
            print(f'Change overlay: {i+1}/{len(old)} old polygons',flush=True)
    rows=[]
    for (scope,t1,s1,t2,s2),area in sorted(totals.items()):
        both=t1 in ('1','2','3','4') and t2 in ('1','2','3','4')
        type_changed=t1!=t2
        # Nonforest subclass reclassification is kept separate from a forest-species change.
        species_changed=both and s1!=s2
        rows.append(dict(scope=scope,type2013_code=t1,type2013=FRTP[t1],species2013_code=s1,species2013=SPECIES[s1],
                         type2025_code=t2,type2025=FRTP[t2],species2025_code=s2,species2025=SPECIES[s2],
                         area_ha=area/10000,type_changed=type_changed,species_code_changed=s1!=s2,
                         both_stocked_forest=both,forest_species_changed=species_changed,
                         type_or_species_code_changed=type_changed or s1!=s2,
                         type_or_forest_species_changed=type_changed or species_changed,
                         nonforest_code_only_changed=(not both and not type_changed and s1!=s2)))
    save_csv(OUT/'forest_change_2013_2025.csv',rows)
    summary=[]
    for scope in ('금정산·캠퍼스 전체','부산대 캠퍼스'):
        subset=[r for r in rows if r['scope']==scope]
        summary.append(dict(scope=scope,common_coverage_ha=sum(r['area_ha'] for r in subset),
                            type_changed_ha=sum(r['area_ha'] for r in subset if r['type_changed']),
                            forest_species_changed_ha=sum(r['area_ha'] for r in subset if r['forest_species_changed']),
                            species_code_changed_ha=sum(r['area_ha'] for r in subset if r['species_code_changed']),
                            type_or_species_code_changed_ha=sum(r['area_ha'] for r in subset if r['type_or_species_code_changed']),
                            type_or_forest_species_changed_ha=sum(r['area_ha'] for r in subset if r['type_or_forest_species_changed']),
                            nonforest_code_only_changed_ha=sum(r['area_ha'] for r in subset if r['nonforest_code_only_changed'])))
    save_csv(OUT/'forest_change_summary.csv',summary)
    return rows,summary

def components(geoms,tolerance):
    parent=list(range(len(geoms)))
    def find(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]]; a=parent[a]
        return a
    tree=shapely.STRtree(geoms)
    edges=tree.query(geoms,predicate='dwithin',distance=tolerance) if tolerance else tree.query(geoms,predicate='intersects')
    for a,b in zip(*edges):
        a,b=find(int(a)),find(int(b))
        if a!=b:
            parent[max(a,b)]=min(a,b)
    return np.array([find(i) for i in range(len(geoms))],dtype=int)

def ridge_model(index,plots):
    central=sorted([p for p in plots if p['zone']=='중앙'],key=lambda p:float(p['위도']))
    # Southern control is an explicitly approximate South-Gate-area reference,
    # not a surveyed gate coordinate. Central plots provide the ridge transect.
    controls=[(129.061,35.235)]+[(float(p['경도']),float(p['위도'])) for p in central]+[(129.0506,35.2801)]
    line=LineString([TO_XY.transform(lon,lat) for lon,lat in controls])
    corridor=line.buffer(100)
    forest_geoms=index.geoms[index.forest_ids]
    seed_mask=shapely.intersects(forest_geoms,corridor)
    strict=components(forest_geoms,0)
    relaxed=components(forest_geoms,10)
    seed0=set(strict[seed_mask]); seed10=set(relaxed[seed_mask])
    by_id={int(fid):i for i,fid in enumerate(index.forest_ids)}
    model=dict(line=line,corridor=corridor,strict=strict,relaxed=relaxed,seed0=seed0,seed10=seed10,
               local_ids=by_id)
    metadata={'reference_controls_lon_lat':controls,'corridor_buffer_m':100,'strict_gap_tolerance_m':0,
              'sensitivity_gap_tolerance_m':10,'strict_component_count':len(set(strict)),
              'sensitivity_component_count':len(set(relaxed)),'seed_polygon_parts':int(seed_mask.sum()),
              'strict_seed_components':sorted(map(int,seed0)),'sensitivity_seed_components':sorted(map(int,seed10)),
              'south_reference_note':'35.235,129.061은 남문 일대 근사 참고점이며 성문 실측 좌표가 아님',
              'reference_method':'남문 일대 참고점-중앙 야장 조사구열-사용자 제공 고당봉 위치를 연결한 참고선. 지형 능선의 실측 선형이 아님',
              'limitations':'폴리곤 연결은 지형 능선 동일성을 보증하지 않음. 10m 결과는 단절 민감도이며 자동 연결 판정에 사용하지 않음.'}
    save_json(OUT/'ridge_connectivity_method.json',metadata)
    save_json(MERGED/'ridge_reference.geojson',{'type':'FeatureCollection','features':[
        {'type':'Feature','properties':{'note':metadata['reference_method']},
         'geometry':{'type':'LineString','coordinates':[list(p) for p in line.coords]}}],
         'crs':{'type':'name','properties':{'name':'EPSG:5179'}}})
    return model,metadata

def ridge_status(index,point,model,hit=None):
    if hit is None:
        hit=index.locate(point)[0]
    forest=hit is not None and stocked(index.features[hit]['attrs'])
    if not forest:
        return dict(ridge_connected_0m=False,ridge_connected_10m=False,ridge_component_id='',
                    ridge_status='2025 숲 폴리곤 밖')
    local=model['local_ids'][hit]
    a,b=int(model['strict'][local]),int(model['relaxed'][local])
    s0,s10=a in model['seed0'],b in model['seed10']
    return dict(ridge_connected_0m=s0,ridge_connected_10m=s10,ridge_component_id=a,
                ridge_status='주능선 참고선 연결' if s0 else ('10m 이내 틈 민감·현장확인' if s10 else '주능선 참고선 비연결'))

def enrich_candidates(index,trails,model):
    source=OUT/'field_sites.csv'
    backup=OUT/'field_sites_gedi_original.csv'
    if not backup.exists():
        shutil.copy2(source,backup)
    # Always regenerate from the original 25, avoiding iterative enrichment drift.
    candidates=read_csv(backup)
    geoms=np.array([f['geometry'] for f in trails],dtype=object)
    tree=shapely.STRtree(geoms)
    rows=[]
    for c in candidates:
        point=Point(*TO_XY.transform(float(c['lon']),float(c['lat'])))
        info,hit=attach_info(index,point,'fm2025_')
        ridge=ridge_status(index,point,model,hit)
        nearest=int(tree.nearest(point)); attrs=trails[nearest]['attrs']
        mountain=clean(attrs.get('MNTN_NM'))
        other=mountain not in ('금정산','인재개발원화명수목원','')
        if c['zone']=='부산대 부산캠퍼스':
            flag='캠퍼스(타 산 판정 제외)'
            recommend=True
            reason='캠퍼스 비교 유지; 임상도 숲 여부는 별도 확인'
        elif other:
            flag=f'다른 산 후보({mountain})'+('·숲 연결은 있음' if ridge['ridge_connected_0m'] else '')
            recommend=False; reason='산림청 등산로 산명 기준으로 금정산 외 산 후보'
        elif not info['fm2025_in_stocked_forest']:
            flag='임상도 숲 밖·GPS/경계 확인'; recommend=False; reason='2025 입목지/죽림 폴리곤 밖'
        elif ridge['ridge_connected_0m']:
            flag='금정산 주능선 연결'; recommend=True; reason='2025 숲 및 주능선 참고선 연결'
        else:
            flag='단절 민감·현장확인' if ridge['ridge_connected_10m'] else '주능선 비연결·타 산/섬숲 확인'
            recommend=False; reason=ridge['ridge_status']
        rows.append(dict(c,**info,**ridge,nearest_trail_mountain=mountain,
                         nearest_trail_distance_m=float(point.distance(geoms[nearest])),
                         other_mountain_flag=flag,recommended_for_visit=recommend,recommendation_reason=reason))
    save_csv(source,rows)
    save_csv(OUT/'field_sites_reselected.csv',[r for r in rows if r['recommended_for_visit']],list(rows[0]))
    return rows

def combine_visits(plots,summaries,comparisons,candidates,index,model):
    lookup={s['plot_id']:s for s in summaries}; comparison={r['plot_id']:r for r in comparisons}
    plot_xy=shapely.points(*[np.array(v) for v in zip(*[TO_XY.transform(float(p['경도']),float(p['위도'])) for p in plots])])
    candidate_xy=shapely.points(*[np.array(v) for v in zip(*[TO_XY.transform(float(p['lon']),float(p['lat'])) for p in candidates])])
    quality=read_csv(OUT/'gedi_v3_quality_points.csv')
    quality_xy=shapely.points(*TO_XY.transform(np.array([float(r['lon']) for r in quality]),
                                              np.array([float(r['lat']) for r in quality])))
    qtree=shapely.STRtree(quality_xy)
    ctree=shapely.STRtree(candidate_xy); ptree=shapely.STRtree(plot_xy)
    geod=Geod(ellps='WGS84')
    def actual_distance(lon,lat,otherlon,otherlat):
        return geod.inv(lon,lat,otherlon,otherlat)[2]
    rows=[]
    pairs=[]
    for i,p in enumerate(plots):
        lon,lat=float(p['경도']),float(p['위도'])
        candidate_ids=[int(j) for j in ctree.query(plot_xy[i],predicate='dwithin',distance=101)
                       if actual_distance(lon,lat,float(candidates[j]['lon']),float(candidates[j]['lat']))<=100]
        quality_ids=[int(j) for j in qtree.query(plot_xy[i],predicate='dwithin',distance=101)
                     if actual_distance(lon,lat,float(quality[j]['lon']),float(quality[j]['lat']))<=100]
        nearest_q=int(qtree.nearest(plot_xy[i])); nearest_c=int(ctree.nearest(plot_xy[i]))
        for j in candidate_ids:
            pairs.append(dict(plot_id=p['plot_id'],gedi_site_id=candidates[j]['site_id'],
                              distance_m=actual_distance(lon,lat,float(candidates[j]['lon']),float(candidates[j]['lat']))))
        s=lookup[p['plot_id']]; comp=comparison[p['plot_id']]
        info={k:v for k,v in comp.items() if k.startswith('fm2025_')}
        ridge=ridge_status(index,plot_xy[i],model)
        rows.append(dict(visit_id=p['plot_id'],record_type='야장 조사구',zone=p['zone'],lat=lat,lon=lon,
                         source_period='2009-2010',dominant_species=s['우점_수종'],mean_dbh_cm=float(s['평균_흉고직경_cm']),
                         mean_height_m=float(s['평균_수고_m']),agbd='',recommended_for_visit=True,
                         comparison_label='GEDI+야장 동시 비교 가능' if quality_ids else '100m 이내 품질 GEDI 없음',
                         comparison_gedi_source='선정 후보25곳' if candidate_ids else ('후보 외 품질점' if quality_ids else ''),
                         gedi_candidate_ids_100m=';'.join(candidates[j]['site_id'] for j in sorted(candidate_ids)),
                         gedi_candidate_count_100m=len(candidate_ids),gedi_quality_count_100m=len(quality_ids),
                         gedi_quality_shots_100m=';'.join(quality[j]['shot_number'] for j in quality_ids),
                         nearest_gedi_candidate_m=actual_distance(lon,lat,float(candidates[nearest_c]['lon']),float(candidates[nearest_c]['lat'])),
                         nearest_gedi_quality_m=actual_distance(lon,lat,float(quality[nearest_q]['lon']),float(quality[nearest_q]['lat'])),
                         nearest_gedi_quality_shot=quality[nearest_q]['shot_number'],paired_plot_ids=p['plot_id'],
                         other_mountain_flag='',google_maps=f'https://maps.google.com/?q={lat:.6f},{lon:.6f}',**info,**ridge))
    for i,c in enumerate(candidates):
        lon,lat=float(c['lon']),float(c['lat'])
        near=[int(j) for j in ptree.query(candidate_xy[i],predicate='dwithin',distance=101)
              if actual_distance(lon,lat,float(plots[j]['경도']),float(plots[j]['위도']))<=100]
        info={k:v for k,v in c.items() if k.startswith('fm2025_') or k.startswith('ridge_')}
        rows.append(dict(visit_id=c['site_id'],record_type='GEDI 후보',zone=c['zone'],lat=lat,lon=lon,
                         source_period=c['date'],dominant_species='',mean_dbh_cm='',mean_height_m='',agbd=float(c['agbd']),
                         recommended_for_visit=c['recommended_for_visit'],comparison_label='GEDI+야장 동시 비교 가능' if near else '100m 이내 야장 없음',
                         comparison_gedi_source='선정 후보25곳' if near else '',gedi_candidate_ids_100m=c['site_id'],
                         gedi_candidate_count_100m=1,gedi_quality_count_100m='',gedi_quality_shots_100m=c['shot_number'],
                         nearest_gedi_candidate_m=0,nearest_gedi_quality_m=0,nearest_gedi_quality_shot=c['shot_number'],
                         paired_plot_ids=';'.join(plots[j]['plot_id'] for j in near),other_mountain_flag=c['other_mountain_flag'],
                         google_maps=c['google_maps'],**info))
    save_csv(OUT/'field_visit_combined.csv',rows)
    save_csv(OUT/'field_visit_pairs_100m.csv',pairs,['plot_id','gedi_site_id','distance_m'])
    return rows,pairs

def combined_kml(rows):
    ns='http://www.opengis.net/kml/2.2'; ET.register_namespace('',ns)
    tag=lambda x:f'{{{ns}}}{x}'
    root=ET.Element(tag('kml')); doc=ET.SubElement(root,tag('Document'))
    ET.SubElement(doc,tag('name')).text='금정산 통합 현장 방문 목록: 야장50 + GEDI25'
    ET.SubElement(doc,tag('description')).text=('야장2009-2010, GEDI2019-2025, 임상도2013/2025의 시점 차이에 유의. '
       '야장 좌표는 원문 측지계 미기재. 숲 연결은 근사 주능선 참고선 기준이며 지형학적 산 경계가 아님.')
    for style_id,color in [('plot','ffff8800'),('gedi','ff0000ff'),('review','ff00aaff')]:
        style=ET.SubElement(doc,tag('Style'),id=style_id); icon=ET.SubElement(style,tag('IconStyle'))
        ET.SubElement(icon,tag('color')).text=color; ET.SubElement(icon,tag('scale')).text='0.8'
    for typ in ('야장 조사구','GEDI 후보'):
        folder=ET.SubElement(doc,tag('Folder')); ET.SubElement(folder,tag('name')).text=typ
        for r in rows:
            if r['record_type']!=typ:
                continue
            mark=ET.SubElement(folder,tag('Placemark')); ET.SubElement(mark,tag('name')).text=r['visit_id']
            style='plot' if typ=='야장 조사구' else ('gedi' if r['recommended_for_visit'] else 'review')
            ET.SubElement(mark,tag('styleUrl')).text='#'+style
            text=[f'{typ}: {r["visit_id"]} / {r["zone"]}',r['comparison_label'],
                  f'연계 야장: {r["paired_plot_ids"]}',f'연계 GEDI 후보: {r["gedi_candidate_ids_100m"]}',
                  f'2025 임상도: {r["fm2025_type"]}, {r["fm2025_species"]}, '
                  f'{r["fm2025_diameter_class"]}, {r["fm2025_age_class"]}, 임분고 코드 {r["fm2025_height_code"]}',
                  f'주능선 연결: {r["ridge_status"]}',r['other_mountain_flag']]
            if typ=='야장 조사구':
                text += [f'야장 우점 {r["dominant_species"]}, 평균 DBH {r["mean_dbh_cm"]:.2f}cm, 평균 수고 {r["mean_height_m"]:.2f}m',
                         f'100m 이내 품질 GEDI {r["gedi_quality_count_100m"]}점 / {r["comparison_gedi_source"]}']
            else:
                text += [f'GEDI AGBD {r["agbd"]:.1f} Mg/ha ({r["source_period"]})']
            ET.SubElement(mark,tag('description')).text='\n'.join(text)
            point=ET.SubElement(mark,tag('Point'))
            ET.SubElement(point,tag('coordinates')).text=f'{r["lon"]:.9f},{r["lat"]:.9f},0'
    ET.ElementTree(root).write(OUT/'field_visit_combined.kml',encoding='utf-8',xml_declaration=True)

def map_figure(features,trails,plots,candidates,model):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection,LineCollection
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    from shapely.plotting import patch_from_polygon
    plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
    totals=defaultdict(float)
    for f in features:
        if stocked(f['attrs']):
            totals[f['attrs']['KOFTR_NM']]+=f['geometry'].area
    names=sorted(totals,key=totals.get,reverse=True)
    palette=list(plt.get_cmap('tab20').colors)+list(plt.get_cmap('tab20b').colors)
    colors={name:palette[i%len(palette)] for i,name in enumerate(names)}
    fig,ax=plt.subplots(figsize=(13,12))
    patches=[]; fills=[]
    for feature in features:
        patches.append(patch_from_polygon(feature['geometry']))
        fills.append(colors[feature['attrs']['KOFTR_NM']] if stocked(feature['attrs']) else '#e4e4e4')
    ax.add_collection(PatchCollection(patches,facecolors=fills,edgecolors='#ffffff',linewidths=.08,zorder=1))
    ax.add_collection(LineCollection([np.array(f['geometry'].coords) for f in trails],colors='#5b4130',linewidths=.55,alpha=.7,zorder=2))
    for zone in ('중앙','동측','서측','남측','북측'):
        subset=[p for p in plots if p['zone']==zone]
        x,y=TO_XY.transform(np.array([float(p['경도']) for p in subset]),np.array([float(p['위도']) for p in subset]))
        ax.scatter(x,y,s=30,marker='^',c='black',edgecolors='white',linewidths=.55,zorder=5)
    x,y=TO_XY.transform(np.array([float(c['lon']) for c in candidates]),np.array([float(c['lat']) for c in candidates]))
    ax.scatter(x,y,s=36,facecolors='none',edgecolors='red',linewidths=1,zorder=6)
    cx,cy=CAMPUS_GEOM.exterior.xy; ax.plot(cx,cy,color='black',ls='--',lw=1.5,zorder=4)
    rx,ry=model['line'].xy; ax.plot(rx,ry,color='white',ls='--',lw=1,zorder=3)
    ax.set(xlim=(CLIP.bounds[0],CLIP.bounds[2]),ylim=(CLIP.bounds[1],CLIP.bounds[3]),
           title='2025 산림청 임상도 · 금정산·부산대 일대\n수종, 등산로, 서정범 조사구 50개와 GEDI 후보 25곳',
           xlabel='동서 좌표 (EPSG:5179, km)',ylabel='남북 좌표 (EPSG:5179, km)')
    ax.set_aspect('equal')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v/1000:.0f}'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v/1000:.0f}'))
    handles=[Patch(facecolor=colors[n],label=f'{n} ({totals[n]/10000:,.0f} ha)') for n in names]
    handles += [Patch(facecolor='#e4e4e4',label='무립목지/비산림'),Line2D([],[],color='#5b4130',label='산림청 등산로'),
                Line2D([],[],marker='^',color='black',linestyle='None',label='야장 조사구 50개'),
                Line2D([],[],marker='o',mfc='none',color='red',linestyle='None',label='GEDI 후보 25곳'),
                Line2D([],[],ls='--',color='black',label='부산대 캠퍼스 범위')]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.015,1),fontsize=8,title='수종·표시',frameon=False)
    fig.text(.5,.025,'폴리곤의 갱신년도는 서로 다름 · 흰색 파선은 주능선 연결 판단용 근사 참고선',ha='center',fontsize=9)
    fig.savefig(OUT/'forest_map_2025.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

def report(metadata,coverage,agreement,change,candidates,visits,pairs):
    def count(year,metric,status):
        return sum(int(r['plots']) for r in agreement if int(r['year'])==year and r['metric']==metric and r['status']==status)
    exact13=count(2013,'정밀 수종','일치'); exact25=count(2025,'정밀 수종','일치')
    denominator13=exact13+count(2013,'정밀 수종','불일치')
    denominator25=exact25+count(2025,'정밀 수종','불일치')
    d13=count(2013,'평균 흉고직경 경급','일치'); d25=count(2025,'평균 흉고직경 경급','일치')
    dn13=d13+count(2013,'평균 흉고직경 경급','불일치'); dn25=d25+count(2025,'평균 흉고직경 경급','불일치')
    h25=count(2025,'평균 수고 등급','일치'); hn25=h25+count(2025,'평균 수고 등급','불일치')
    lines=['# 금정산 임상도·야장·GEDI 통합 분석','',
           '탄소·바이오매스 계산은 하지 않았다. 원본 SHP는 보존하고 잘라낸 결과만 merged 폴더에 저장했다.','',
           '## 한국어 요약','',
           f'- 수종명 직접 비교가 가능한 곳은 2013년 {exact13}/{denominator13}곳, 2025년 {exact25}/{denominator25}곳이 일치했다. 기타수종·혼효림 묶음은 정밀 일치율에서 제외했다.',
           f'- 야장 평균 흉고직경을 정의서 경급 구간과 비교하면 2013년 {d13}/{dn13}곳, 2025년 {d25}/{dn25}곳이 일치했다.',
           f'- 야장 평균 수고가 2025 임분고 등급 범위 안에 든 곳은 {h25}/{hn25}곳이었다. 2013 자료에는 HEIGHT가 없다.',
           '- 야장(2009–2010)과 지도 갱신시점 및 수관점유 기준·매목 평균의 차이가 있어, 불일치를 지도 오류로 단정할 수 없다.','',
           '## 임상도와 현장 측정 비교','',
           '| 연도 | 비교 항목 | 결과 | 조사구 수 |','|---|---|---|---:|']
    for r in agreement:
        lines.append(f'| {r["year"]} | {r["metric"]} | {r["status"]} | {r["plots"]} |')
    lines += ['', '수종의 정밀 일치와 기타활엽수·기타참나무류 등의 상위군 호환은 구분했다. '
              '혼효림은 단일 우점 수종과 직접 비교하지 않았다.',
              '야장 수종은 흉고단면적 기준 우점 수종이다. 원문 오탈자 후보는 교정하지 않았다.',
              '경급은 원래 수관점유면적 51% 이상 입목의 직경구간이므로 전체 매목 평균 DBH와의 일치는 참고 지표다.',
              '경급 경계: 치수 <6cm, 소경목 6≤DBH<18, 중경목 18≤DBH<30, 대경목 ≥30.',
              '임분고 02~38은 코드±1m 범위, 00은 <1m, 40은 ≥39m이다. 2013 HEIGHT는 결측으로 유지했다.',
              '야장은 2009–2010년, 지도는 2013/2025 배포본이다. 2025 파일의 개별 갱신년도는 2015–2025년 등이 섞여 있다.',
              '따라서 이 결과는 동일 시점·동일 공간 지지범위의 정확도 검증이 아니며, 불일치를 지도 오류로 단정할 수 없다.', '',
              '## 면적·피복','', '| 연도 | 범위 | 범위 면적(ha) | 지도 피복(ha) | 입목지·죽림(ha) | 미피복(ha) |',
              '|---|---|---:|---:|---:|---:|']
    for r in coverage:
        lines.append(f'| {r["year"]} | {r["scope"]} | {r["region_area_ha"]:.2f} | {r["mapped_area_ha"]:.2f} | '
                     f'{r["stocked_forest_area_ha"]:.2f} | {r["unmapped_area_ha"]:.2f} |')
    lines += ['', '면적은 EPSG:5179 도형으로 다시 계산했다. 네 분류 항목의 면적 합계를 서로 더하면 중복 집계가 된다.',
              '전체 자르기 범위는 WGS84 BBOX를 변환한 EPSG:5179 축 정렬 외접 사각형이며, 캠퍼스 면적은 경위도 사각형 경계를 투영한 도형으로 계산했다.',
              '세 구·시 자료에 없는 지역을 비산림으로 간주하지 않고 미피복으로 남겼다.', '', '## 2013→2025 변화', '',
              '| 범위 | 공통 피복(ha) | 임상 변화(ha) | 양쪽 숲의 수종 코드 변화(ha) | 임상·전체 수종코드 중 하나 변화(ha) | 비입목 코드만 변경 제외(ha) |',
              '|---|---:|---:|---:|---:|---:|']
    for r in change:
        lines.append(f'| {r["scope"]} | {r["common_coverage_ha"]:.2f} | {r["type_changed_ha"]:.2f} | '
                     f'{r["forest_species_changed_ha"]:.2f} | {r["type_or_species_code_changed_ha"]:.2f} | '
                     f'{r["type_or_forest_species_changed_ha"]:.2f} |')
    lines += ['', '같은 위치의 폴리곤 교차면적으로 계산했다. 2013 주거지·경작지 등과 2025 제지 등의 비입목 코드 재분류는 별도 열에 분리했다.',
              '공통 피복 밖은 변화량 분모에서 제외했다. 지도 경계·분류체계·갱신시점 차이가 포함되므로 실제 벌채·천이 면적과 같지 않다.', '',
              '## 후보 재검토와 연결 판단', '',
              '기존 25곳은 field_sites_gedi_original.csv에 보존하고 field_sites.csv에 2025 속성과 판단 열을 붙였다.',
              '재선정은 기존 후보 중에서만 수행했으며 새로운 GEDI 지점을 보충하지 않았다.',
              f'- 재검토 우선 후보: {sum(r["recommended_for_visit"] for r in candidates)}곳 '
              f'(캠퍼스는 비교 목적상 유지, 임상도 숲 밖 여부 별도 표시)',
              f'- 기존 GEDI 후보와 야장 사이 100m 이내 쌍: {len(pairs)}쌍',
              f'- 100m 이내 품질 통과 GEDI가 있는 야장 조사구: '
              f'{sum(r["record_type"]=="야장 조사구" and r["comparison_label"]=="GEDI+야장 동시 비교 가능" for r in visits)}곳',
              '- 100m 표시는 후보25곳과 전체 품질 통과 GEDI를 구분한 열을 제공한다. 후보 외 점은 통합 75행에 추가하지 않고 shot number로 기록했다.',
              '- 거리 100m는 WGS84 측지거리로 최종 판정했다. STRtree 101m 검색 뒤 정확한 거리를 확인했다.',
              '- 야장 20×20m와 GEDI 약25m footprint는 중심 100m 이내라도 정확히 같은 면적이 아니다.',
              '- 원문 야장 GPS는 초 단위이며 측지계 미기재이므로 위치 대응은 근사다.', '',
              '연결 판단은 2025 입목지·죽림 폴리곤의 접촉 연결(0m)을 기본으로 했다. 10m 이내 틈 허용 결과는 민감도 열이며 기본 판정을 대체하지 않았다.',
              '주능선 참고선은 고당봉(사용자 제공), 중앙 야장 조사구열, 남문 일대 근사 참고점(35.235,129.061)을 연결하고 100m 완충한 것이다. '
              '남문 참고점은 성문 실측 좌표가 아니며 이 선은 지형학적 능선을 정밀 추출한 결과가 아니다.',
              '산림청 등산로의 가장 가까운 산명(윤산·백양산 등)을 별도로 대조했다. 숲이 연결돼 있어도 산명이 다르면 다른 산 후보로 표시했다.',
              '산림청 자료에 나타난 산명은 지점의 산 소속을 확정하는 행정경계가 아니며, 주능선 비연결만으로 다른 산임을 단정하지 않는다.', '',
              '| 후보 | 구역 | 숲 안 | 산림청 인접 등산로 산명 | 연결 | 재검토 |', '|---|---|---|---|---|---|']
    for r in candidates:
        lines.append(f'| {r["site_id"]} | {r["zone"]} | {r["fm2025_in_stocked_forest"]} | {r["nearest_trail_mountain"]} | '
                     f'{r["ridge_status"]} | {r["other_mountain_flag"]} |')
    lines += ['', '지리 참고: [금정구 금정산성 안내](https://www.geumjeong.go.kr/index.geumj?menuCd=DOM_000000141002003000)',
              '코드 출처: data/forest_map/2013_임상도/금정구/DATA003_임상도(5000).pdf, 4–9쪽 및 2025 동명 정의서.', '',
              '상세 파일: forest_area_summary.csv, forest_coverage_summary.csv, seo2011_vs_forestmap.csv, '
              'forest_field_agreement_summary.csv, forest_change_2013_2025.csv, forest_change_summary.csv, '
              'field_sites.csv, field_sites_reselected.csv, field_visit_combined.csv/.kml, field_visit_pairs_100m.csv.', '',
              '공간 처리·연결 기준·검증 정보는 forest_analysis_metadata.json과 ridge_connectivity_method.json에 보존했다.']
    lines += ['', '## 중복 피복 처리', '',
              '면적 중복은 금정구→북구→양산시 및 원본 레코드 순으로 앞선 속성을 유지해 제거했다. '
              '중복 면적을 더하지 않기 위한 처리이며 속성의 우열을 보증하지 않는다. '
              '겹침별 속성 충돌은 forest_overlap_2013.csv, forest_overlap_2025.csv에 기록했다.']
    for year in ('2013','2025'):
        overlap=metadata['processing_logs'][year][-1]
        lines.append(f'- {year}: 중복 제거 {overlap["positive_overlap_removed_m2"]/10000:.4f}ha, '
                     f'분류 속성 충돌 영역 {overlap["conflicting_overlap_m2"]/10000:.4f}ha.')
    (OUT/'forest_analysis_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
    sys.stdout.reconfigure(encoding='utf-8')
    start=time.monotonic()
    MERGED.mkdir(parents=True,exist_ok=True)
    layers={}; logs={}
    for year in (2013,2025):
        features,schema,log=load_layer(f'{year}_임상도')
        write_layer(features,schema,MERGED/f'forest{year}_geumjeong.shp')
        layers[year]=features; logs[str(year)]=log
    trails,schema,log=load_layer('등산로','LineString')
    write_layer(trails,schema,MERGED/'trails_geumjeong.shp','LineString'); logs['trails']=log
    write_code_mapping()
    _,coverage=area_summary(layers)
    indices={year:ForestIndex(features) for year,features in layers.items()}
    plots=read_csv(ROOT/'data/field/seo2011_plots.csv')
    summaries=read_csv(ROOT/'data/field/seo2011_plot_summary.csv')
    comparison,agreement=compare_plots(indices,plots,summaries)
    print('50 field plots compared',flush=True)
    _,change=transitions(layers[2013],layers[2025])
    print('2013/2025 transitions complete',flush=True)
    model,ridge_meta=ridge_model(indices[2025],plots)
    candidates=enrich_candidates(indices[2025],trails,model)
    visits,pairs=combine_visits(plots,summaries,comparison,candidates,indices[2025],model)
    combined_kml(visits)
    map_figure(layers[2025],trails,plots,candidates,model)
    metadata=dict(bbox_lon_lat=BBOX,clip_epsg5179_bounds=CLIP.bounds,campus_lon_lat=CAMPUS,
                  layer_counts={str(y):len(fs) for y,fs in layers.items()},trail_parts=len(trails),processing_logs=logs,
                  field_plots=len(plots),gedi_candidates=len(candidates),combined_visits=len(visits),pair_count_100m=len(pairs),
                  field_plots_near_quality_gedi=sum(r['record_type']=='야장 조사구' and r['gedi_quality_count_100m']>0 for r in visits),
                  recommended_gedi_candidates=sum(r['recommended_for_visit'] for r in candidates),
                  ridge_method=ridge_meta,elapsed_seconds=time.monotonic()-start,
                  definitions={str(y):str(next((RAW/f'{y}_임상도'/'금정구').glob('DATA003*.pdf')).relative_to(ROOT)) for y in layers},
                  input_shp_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                    for period in ('2013_임상도','2025_임상도','등산로') for p in (RAW/period).glob('*/*.shp')},
                  carbon_or_biomass_calculated=False)
    save_json(OUT/'forest_analysis_metadata.json',metadata)
    report(metadata,coverage,agreement,change,candidates,visits,pairs)
    assert len(plots)==50 and len(candidates)==25 and len(visits)==75
    assert not {'pandas','geopandas','rasterio','pyogrio'} & set(sys.modules)
    print(json.dumps({'layer_counts':metadata['layer_counts'],'trails':len(trails),'comparison':agreement,
                      'change_summary':change,'recommended':metadata['recommended_gedi_candidates'],
                      'combined':len(visits),'candidate_plot_pairs_100m':len(pairs),
                      'plots_near_quality_gedi':metadata['field_plots_near_quality_gedi'],
                      'elapsed_seconds':metadata['elapsed_seconds']},ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    main()
