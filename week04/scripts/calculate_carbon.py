"""Local forest-carbon comparison. No dataframe/GDAL libraries or uploads.

Allometry: kg dry matter, D cm/H m. Maps: EPSG:5179, m2 -> ha.
The CSVs supplied in data/reference and their README are authoritative inputs.
"""
from pathlib import Path
from collections import defaultdict, Counter
from functools import lru_cache
import csv, json, math, sys, hashlib, platform
import numpy as np
import scipy
import shapely
from shapely.geometry import Point, LineString, mapping
from shapely.ops import transform
from pyproj import Geod
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch, Patch
from matplotlib.path import Path as MplPath
from forest_map_common import ROOT, MERGED, OUT, CLIP, CAMPUS_GEOM, TO_XY, TO_LL, read_csv, save_csv, save_json, parts, stocked
from analyze_forest_maps import components
from carbon_audit import read_merged

REF=ROOT/'data/reference'
PLOT_AREA_HA=0.04
G='금정산'; U='부산대 부산캠퍼스'; W='클리핑 전체'
TYPES={'1':'침엽수림','2':'활엽수림','3':'혼효림'}
CO2=44/12
assumptions=[]; tests=[]
def assume(category,key,value,reason,source='data/reference/README.md'):
    assumptions.append(dict(category=category,key=key,value=value,reason=reason,source=source))
def check(name,condition,detail):
    tests.append(dict(test=name,passed=bool(condition),detail=str(detail)))
    if not condition:
        save_csv(OUT/'carbon_validation.csv',tests)
        raise AssertionError(name+': '+str(detail))
def number(x):
    return float(x) if x not in ('',None) else None
def stats(values):
    a=np.asarray(values,dtype=float)
    return dict(mean=float(a.mean()),median=float(np.median(a)),min=float(a.min()),max=float(a.max()))

FACT={r['species']:r for r in read_csv(REF/'nifos2014_carbon_factors.csv')}
CF={r['species']:float(r['CF']) for r in read_csv(REF/'nifos2010_carbon_fraction.csv')}
EQ={(r['species'],r['form'],r['part']):r for r in read_csv(REF/'nifos2014_allometry.csv')}
VOLUME={(r['region'],r['forest_type']):r for r in read_csv(REF/'kosis_nfi2020_volume_per_ha.csv')}
CON=[s for s,r in FACT.items() if r['forest_type']=='침엽수']
DEC=[s for s,r in FACT.items() if r['forest_type']=='낙엽활엽수']
OAK=['상수리나무','굴참나무','신갈나무','졸참나무']
ALIASES={'낙엽송':'일본잎갈나무','편백나무':'편백','개서어나무':'서어나무','소사나무':'서어나무',
         '백합나무':'튤립나무','포플러':'현사시나무'}
# 백합나무 is the official Korean synonym for 튤립나무; generic 포플러 is
# deliberately NOT mapped to 현사시, whose species identity is more specific.
ALIASES.pop('포플러')

def targets(species,dbh=None,forest_code=None):
    if species=='소나무':
        return ['강원지방소나무' if dbh is not None and dbh>40 else '중부지방소나무'],'소나무 직경 규칙'
    if species in ('떡갈나무','기타참나무류'):
        return OAK,'참나무4종 평균'
    canonical=ALIASES.get(species,species)
    if canonical in FACT:
        return [canonical],'직접/명시 별칭'
    if forest_code=='3' or species in ('침활혼효림','혼효림'):
        # Equal group weights, not equal weights across all 18 species.
        return None,'침엽·낙엽활엽 그룹 평균의 평균'
    if forest_code=='1' or species in ('노간주나무','기타침엽수','가문비나무','비자나무','은행나무','잔나무'):
        return CON,'기타 침엽수 평균'
    return DEC,'기타 활엽수→낙엽활엽수 평균'

def weights(species,dbh=None,forest_code=None):
    ts,rule=targets(species,dbh,forest_code)
    if ts is None:
        return [(s,0.5/len(CON)) for s in CON]+[(s,0.5/len(DEC)) for s in DEC],rule
    return [(s,1/len(ts)) for s in ts],rule

@lru_cache(None)
def factors(species,forest_code):
    ws,rule=weights(species,None,forest_code)
    result={k:sum(w*float(FACT[s][k]) for s,w in ws) for k in ('WD','BEF','R')}
    result['CF_measured']=sum(w*CF.get(s,0.5) for s,w in ws)
    result['CF_fallback_weight']=sum(w for s,w in ws if s not in CF)
    result['mapping']=rule
    result['targets']=';'.join(f'{s}:{w:.8g}' for s,w in ws)
    return result

def allometry(species,d,h,form):
    ws,rule=weights(species,d)
    masses={}; flagged=[]
    for part in ('줄기','가지','잎','뿌리'):
        mass=0
        for s,w in ws:
            row=EQ[s,form,part]
            if not float(row['DBH_min_cm'])<=d<=float(row['DBH_max_cm']):
                flagged.append(s)
            v=float(row['a'])*d**float(row['b'])
            if form=='a*D^b*H^c': v*=h**float(row['c'])
            mass+=w*v
        masses[part]=mass
    ag=sum(masses[k] for k in ('줄기','가지','잎')); root=masses['뿌리']
    cf=sum(w*CF.get(s,0.5) for s,w in ws)
    return dict(ag=ag,root=root,carbon=(ag+root)*0.5,carbon_cf=(ag+root)*cf,
                cf=cf,cf_fallback=sum(w for s,w in ws if s not in CF),
                mapping=rule,targets=';'.join(s for s,w in ws),flag=';'.join(sorted(set(flagged))),**masses)

def calculate_A(plots):
    raw=read_csv(ROOT/'data/field/seo2011_trees.csv'); rows=[]
    for i,t in enumerate(raw,1):
        s=t['수종']; d=float(t['흉고직경_cm']); h=float(t['수고_m'])
        excluded='등나무(덩굴)' if s=='등나무' else ('DBH<6cm' if d<6 else '')
        r=dict(tree_id=i,**t,included=not bool(excluded),exclusion_reason=excluded)
        if not excluded:
            a=allometry(s,d,h,'a*D^b'); b=allometry(s,d,h,'a*D^b*H^c')
            r.update(mapping=a['mapping'],target_species=a['targets'],out_of_range_targets=a['flag'],
                     CF_measured_or_fallback=a['cf'],CF_fallback_weight=a['cf_fallback'])
            for key,part in [('stem_kg','줄기'),('branch_kg','가지'),('leaf_kg','잎'),('root_kg','뿌리')]:
                r[key]=a[part]
            r.update(AGB_kg=a['ag'],carbon_kg=a['carbon'],carbon_measured_CF_kg=a['carbon_cf'],
                     AGB_DH_kg=b['ag'],root_DH_kg=b['root'],carbon_DH_kg=b['carbon'],carbon_DH_measured_CF_kg=b['carbon_cf'])
        rows.append(r)
    output=[]
    for p in plots:
        members=[t for t in rows if t['plot_id']==p['plot_id']]; used=[t for t in members if t['included']]
        r=dict(p,plot_area_ha=PLOT_AREA_HA,raw_tree_count=len(members),included_tree_count=len(used),
               small_dbh_count=sum(float(t['흉고직경_cm'])<6 for t in members),
               liana_count=sum(t['수종']=='등나무' for t in members),
               extrapolated_tree_count=sum(bool(t['out_of_range_targets']) for t in used))
        for column,outname in [('AGB_kg','AGB_t_ha'),('root_kg','root_t_ha'),('carbon_kg','carbon_tC_ha'),
              ('carbon_measured_CF_kg','carbon_measured_CF_tC_ha'),('AGB_DH_kg','AGB_DH_t_ha'),
              ('root_DH_kg','root_DH_t_ha'),('carbon_DH_kg','carbon_DH_tC_ha'),
              ('carbon_DH_measured_CF_kg','carbon_DH_measured_CF_tC_ha')]:
            r[outname]=sum(t[column] for t in used)/1000/PLOT_AREA_HA
        r['carbon_total_tC']=r['carbon_tC_ha']*PLOT_AREA_HA
        r['CO2_tCO2_ha']=r['carbon_tC_ha']*CO2
        r['CO2_total_tCO2']=r['carbon_total_tC']*CO2
        r['AGB_total_t']=r['AGB_t_ha']*PLOT_AREA_HA
        output.append(r)
    save_csv(OUT/'carbon_A_trees.csv',rows); save_csv(OUT/'carbon_A_plots.csv',output)
    save_csv(OUT/'carbon_A_out_of_range.csv',[r for r in rows if r['included'] and r['out_of_range_targets']],list(rows[0])+[k for k in rows[next(i for i,r in enumerate(rows) if r['included'])] if k not in rows[0]])
    group=[]
    selections={'전체50':output,**{z:[r for r in output if r['zone']==z] for z in ('중앙','동측','서측','남측','북측')},
                '장전동 인접 남측01~05':[r for r in output if r['plot_id'] in [f'남측-{i:02}' for i in range(1,6)]]}
    for label,rs in selections.items():
        for col in ('AGB_t_ha','carbon_tC_ha','carbon_DH_tC_ha','carbon_measured_CF_tC_ha'):
            group.append(dict(group=label,metric=col,n_plots=len(rs),**stats([r[col] for r in rs]),
                              measured_area_ha=len(rs)*PLOT_AREA_HA,
                              measured_carbon_total_tC=sum(r['carbon_total_tC'] for r in rs)))
    save_csv(OUT/'carbon_A_summary.csv',group)
    mapping_rows=[]
    for s in sorted({r['수종'] for r in rows}):
        used=[r for r in rows if r['수종']==s and r['included']]
        mapping_rows.append(dict(species=s,raw_count=sum(r['수종']==s for r in rows),included_count=len(used),
            carbon_tC=sum(r['carbon_kg'] for r in used)/1000,
            carbon_share_pct=100*sum(r['carbon_kg'] for r in used)/sum(r['carbon_kg'] for r in rows if r['included']),
            mapping=';'.join(sorted({r['mapping'] for r in used})),
            target_species=';'.join(sorted({r['target_species'] for r in used}))))
        assume('A 수종 대응',s,mapping_rows[-1]['target_species'] or '제외',mapping_rows[-1]['mapping'] or '등나무/소경목 제외')
    save_csv(OUT/'carbon_A_species_contributions.csv',mapping_rows)
    check('A 조사구50·원자료1750',len(output)==50 and len(raw)==1750,(len(output),len(raw)))
    check('A 나무kg→조사구tC 일치',math.isclose(sum(r['carbon_kg'] for r in rows if r['included'])/1000,
          sum(r['carbon_total_tC'] for r in output),rel_tol=1e-12),'나무합/1000=조사구합')
    return rows,output,group,mapping_rows

def boundary(features,trails,plots):
    ids=[i for i,f in enumerate(features) if stocked(f['attrs'])]
    gs=np.array([features[i]['geometry'] for i in ids],dtype=object)
    strict=components(gs,0); tree=shapely.STRtree(gs)
    votes=defaultdict(lambda:defaultdict(float))
    for t in trails:
        name=t['attrs'].get('MNTN_NM','').strip()
        for i in tree.query(t['geometry'],predicate='intersects'):
            votes[int(strict[i])][name]+=float(gs[i].intersection(t['geometry']).length)
    meta=json.loads((OUT/'ridge_connectivity_method.json').read_text(encoding='utf-8'))
    line=LineString([TO_XY.transform(*xy) for xy in meta['reference_controls_lon_lat']])
    seeds=shapely.intersects(gs,line.buffer(100)); seed_components=set(strict[seeds])
    dominant={c:max(v,key=v.get) if v else '' for c,v in votes.items()}
    other={c for c,name in dominant.items() if name not in ('','금정산','인재개발원화명수목원') and c not in seed_components}
    eligible=np.array([i for i,c in enumerate(strict) if c not in other],dtype=int)
    relaxed=components(gs[eligible],20); selected_components=set(relaxed[seeds[eligible]])
    selected_local=eligible[np.isin(relaxed,list(selected_components))]
    selected_ids={ids[i] for i in selected_local}
    geom=shapely.union_all(gs[selected_local])
    bare=components(gs,20); bare_ids=np.isin(bare,list(set(bare[seeds])))
    strict_geom=shapely.union_all(gs[np.isin(strict,list(seed_components))])
    memberships=[]
    for j,fid in enumerate(ids):
        memberships.append(dict(polygon_id=features[fid]['attrs']['FID_U'],strict_component=int(strict[j]),
            component_dominant_trail_mountain=dominant.get(int(strict[j]),''),
            excluded_other_mountain=int(strict[j]) in other,selected_geumjeong_20m=fid in selected_ids,
            area_ha=gs[j].area/10000))
    save_csv(OUT/'carbon_boundary_membership.csv',memberships)
    excluded_rows=[]
    for c in sorted(other):
        excluded_rows.append(dict(component=int(c),dominant_mountain=dominant[c],area_ha=float(shapely.area(gs[strict==c]).sum())/10000,
                                  intersecting_trail_lengths_m=json.dumps(votes[c],ensure_ascii=False)))
    save_csv(OUT/'carbon_boundary_excluded_components.csv',excluded_rows)
    metadata=dict(method='2025 숲0m 연결 성분의 등산로 교차길이 최다 산명으로 타 산 제외 후 20m 연결',
        seed_controls_lon_lat=meta['reference_controls_lon_lat'],seed_line_buffer_m=100,gap_m=20,
        strict_seed_area_ha=strict_geom.area/10000,selected_area_including_bamboo_ha=geom.area/10000,
        bare_20m_area_before_other_mountain_exclusion_ha=float(shapely.area(gs[bare_ids]).sum())/10000,
        excluded_component_count=len(other),selected_polygon_count=len(selected_ids),
        limitations='공식 공원/지형 경계 아님; 산명 없는 작은 성분은20m 연결로 포함될 수 있음. 2013도2025 고정경계를 사용.')
    save_json(OUT/'carbon_boundary_method.json',metadata)
    geo={'type':'FeatureCollection','features':[{'type':'Feature','properties':{'name':'금정산 근사 숲 경계','gap_m':20},
              'geometry':mapping(transform(TO_LL.transform,geom))}]}
    save_json(OUT/'carbon_geumjeong_boundary.geojson',geo)
    assume('공간','금정산 경계',json.dumps(metadata,ensure_ascii=False),'타 산 산명 성분을 사전 제외. 빈20m 틈 자체 면적은 더하지 않음','outputs/carbon_boundary_method.json')
    return geom,strict_geom,metadata,selected_ids

@lru_cache(None)
def density(species,forest_code,age_code,district,volume_override=None):
    if forest_code not in TYPES:
        return dict(AGB_t_ha=0.,root_t_ha=0.,carbon_tC_ha=0.,carbon_measured_CF_tC_ha=0.,volume_m3_ha=0.,
                    R=0.,WD=0.,BEF=0.,CF_measured=0.5,CF_fallback_weight=0.,age_substitution=False,
                    region='',mapping='죽림/무립목지/비산림 0',targets='')
    if not age_code: raise ValueError('숲인데 영급 없음: '+str((species,district)))
    age=int(age_code); age_used=min(6,max(2,age))
    if district not in ('금정구','북구','양산시'): raise ValueError('Unknown district '+district)
    region=volume_override or ('경상남도' if district=='양산시' else '부산광역시')
    volrow=VOLUME[region,TYPES[forest_code]]
    v=float(volrow['age6plus_m3_ha' if age_used==6 else f'age{age_used}_m3_ha'])
    f=factors(species,forest_code); ag=v*f['WD']*f['BEF']; root=ag*f['R']
    return dict(f,AGB_t_ha=ag,root_t_ha=root,carbon_tC_ha=(ag+root)*0.5,
                carbon_measured_CF_tC_ha=(ag+root)*f['CF_measured'],volume_m3_ha=v,
                age_substitution=age==1,age_used='6+' if age_used==6 else str(age_used),region=region)

def feature_density(f,override=None):
    a=f['attrs']
    return density(a['KOFTR_NM'],a['FRTP_CD'],a['AGCLS_CD'],a['SRC_DIST'],override)

def calculate_B(layers,scopes,selected):
    records={}; summaries=[]; totals={}; area_arrays={}; mapping_shares=[]
    ridge_parts=np.array([layers[2025][i]['geometry'] for i in sorted(selected)],dtype=object)
    ridge_tree=shapely.STRtree(ridge_parts)
    for year,fs in layers.items():
        print('B polygons',year,flush=True)
        areas={s:np.array([f['geometry'].intersection(g).area/10000 if f['geometry'].intersects(g) else 0 for f in fs]) for s,g in scopes.items() if s!=G}
        if year==2025:
            areas[G]=np.array([f['geometry'].area/10000 if i in selected else 0 for i,f in enumerate(fs)])
        else:
            areas[G]=np.array([sum(f['geometry'].intersection(ridge_parts[int(i)]).area for i in
                             ridge_tree.query(f['geometry'],predicate='intersects'))/10000 for f in fs])
        area_arrays[year]=areas
        ds=[feature_density(f) for f in fs]
        rows=[]
        for i,(f,d) in enumerate(zip(fs,ds)):
            a=f['attrs']; r=dict(year=year,polygon_id=a['FID_U'],district=a['SRC_DIST'],forest_code=a['FRTP_CD'],
                forest_type=a['FRTP_NM'],species=a['KOFTR_NM'],age_class=a['AGCLS_NM'],diameter_class=a['DMCLS_NM'],
                map_update_year=a['갱신년도'],height_code=a['HEIGHT'],area_ha=f['geometry'].area/10000,**d)
            r.update(AGB_total_t=d['AGB_t_ha']*r['area_ha'],root_total_t=d['root_t_ha']*r['area_ha'],
                     volume_total_m3=d['volume_m3_ha']*r['area_ha'],carbon_total_tC=d['carbon_tC_ha']*r['area_ha'],
                     CO2_total_tCO2=d['carbon_tC_ha']*r['area_ha']*CO2,CO2_tCO2_ha=d['carbon_tC_ha']*CO2)
            for s,tag in ((G,'geumjeong'),(U,'campus'),(W,'bbox')):
                r[tag+'_area_ha']=areas[s][i]
                r[tag+'_carbon_tC']=areas[s][i]*d['carbon_tC_ha']
            rows.append(r)
        records[year]=rows
        save_csv(OUT/f'carbon_B_polygons_{year}.csv',rows)
        for scope in scopes:
            for dim,key in [('total',None),('forest_type','forest_type'),('species','species')]:
                cats=['전체'] if key is None else sorted({r[key] for i,r in enumerate(rows) if areas[scope][i]>0})
                for cat in cats:
                    ii=[i for i,r in enumerate(rows) if areas[scope][i]>0 and (key is None or r[key]==cat)]
                    if not ii: continue
                    area=sum(areas[scope][i] for i in ii)
                    forest_area=sum(areas[scope][i] for i in ii if rows[i]['forest_code'] in TYPES)
                    bamboo=sum(areas[scope][i] for i in ii if rows[i]['forest_code']=='4')
                    r=dict(year=year,scope=scope,group_dimension=dim,group=cat,mapped_area_ha=area,
                           forest_area_ha=forest_area,bamboo_zero_area_ha=bamboo,polygon_count=len(ii))
                    for col,out in [('volume_m3_ha','volume_total_m3'),('AGB_t_ha','AGB_total_t'),('root_t_ha','root_total_t'),
                         ('carbon_tC_ha','carbon_total_tC'),('carbon_measured_CF_tC_ha','carbon_measured_CF_total_tC')]:
                        r[out]=sum(areas[scope][i]*rows[i][col] for i in ii)
                    r['AGB_t_ha']=r['AGB_total_t']/forest_area if forest_area else ''
                    r['carbon_tC_ha']=r['carbon_total_tC']/forest_area if forest_area else ''
                    r['CO2_total_tCO2']=r['carbon_total_tC']*CO2
                    r['CO2_tCO2_ha']=r['carbon_tC_ha']*CO2 if forest_area else ''
                    r['R_area_weighted']=sum(areas[scope][i]*rows[i]['R'] for i in ii if rows[i]['forest_code'] in TYPES)/forest_area if forest_area else ''
                    r['CF_measured_area_weighted']=sum(areas[scope][i]*rows[i]['CF_measured'] for i in ii if rows[i]['forest_code'] in TYPES)/forest_area if forest_area else ''
                    r['age1_substitution_area_ha']=sum(areas[scope][i] for i in ii if rows[i]['age_substitution'])
                    r['CF_fallback_equivalent_area_ha']=sum(areas[scope][i]*rows[i]['CF_fallback_weight'] for i in ii if rows[i]['forest_code'] in TYPES)
                    summaries.append(r)
                    if dim=='total': totals[year,scope]=r
        for spec,fc in sorted({(f['attrs']['KOFTR_NM'],f['attrs']['FRTP_CD']) for f in fs}):
            d=next(d for f,d in zip(fs,ds) if f['attrs']['KOFTR_NM']==spec and f['attrs']['FRTP_CD']==fc)
            assume('B 수종 대응',f'{year}:{fc}:{spec}',d['targets'],d['mapping'],'nifos2014_carbon_factors.csv')
    save_csv(OUT/'carbon_B_summary.csv',summaries)
    for year,rows in records.items():
        for scope in scopes:
            total=totals[year,scope]
            for rule in sorted({r['mapping'] for r in rows if r['forest_code'] in TYPES}):
                ids=[i for i,r in enumerate(rows) if r['mapping']==rule and r['forest_code'] in TYPES]
                ar=sum(area_arrays[year][scope][i] for i in ids)
                carbon=sum(area_arrays[year][scope][i]*rows[i]['carbon_tC_ha'] for i in ids)
                mapping_shares.append(dict(year=year,scope=scope,mapping=rule,forest_area_ha=ar,
                    area_share_pct=100*ar/total['forest_area_ha'],carbon_total_tC=carbon,
                    carbon_share_pct=100*carbon/total['carbon_total_tC']))
    save_csv(OUT/'carbon_B_mapping_contributions.csv',mapping_shares)
    changes=[]
    for scope in scopes:
        old=totals[2013,scope]; new=totals[2025,scope]
        changes.append(dict(scope=scope,forest2013_ha=old['forest_area_ha'],forest2025_ha=new['forest_area_ha'],
            carbon2013_tC=old['carbon_total_tC'],carbon2025_tC=new['carbon_total_tC'],
            difference_tC=new['carbon_total_tC']-old['carbon_total_tC'],difference_pct=100*(new['carbon_total_tC']/old['carbon_total_tC']-1),
            CO2_2013_tCO2=old['CO2_total_tCO2'],CO2_2025_tCO2=new['CO2_total_tCO2'],
            caveat='지도 분류·갱신·피복 차이와 공통2020 축적 적용; 실측 성장량/흡수량 아님; 금정산은2025 고정 숲 경계'))
    save_csv(OUT/'carbon_B_2013_2025.csv',changes)
    return records,summaries,totals,area_arrays

def calculate_C(fs,scopes,totals,areas):
    gs=np.array([f['geometry'] for f in fs],dtype=object); tree=shapely.STRtree(gs)
    raw=read_csv(OUT/'gedi_v3_quality_points.csv'); inside=[]
    for q in raw:
        ag=float(q['agbd']); point=Point(*TO_XY.transform(float(q['lon']),float(q['lat'])))
        hits=[int(i) for i in tree.query(point,predicate='intersects') if fs[i]['attrs']['FRTP_CD'] in TYPES]
        if not hits: continue
        check('GEDI 유효AGBD '+q['shot_number'],math.isfinite(ag) and ag>=0,ag) if not(math.isfinite(ag) and ag>=0) else None
        i=min(hits); f=fs[i]
        inside.append(dict(q,polygon_id=f['attrs']['FID_U'],forest_type=TYPES[f['attrs']['FRTP_CD']],
                           geometry=point,scope_geumjeong=scopes[G].intersects(point),scope_campus=scopes[U].intersects(point)))
    summaries=[]; lookup={}
    for scope in scopes:
        used=[r for r in inside if scopes[scope].intersects(r['geometry'])]
        for typ in ['전체']+list(TYPES.values()):
            qs=[r for r in used if typ=='전체' or r['forest_type']==typ]
            if typ=='전체': b=totals[2025,scope]
            else:
                ii=[i for i,f in enumerate(fs) if f['attrs']['FRTP_CD'] in TYPES and TYPES[f['attrs']['FRTP_CD']]==typ]
                ar=areas[2025][scope][ii]
                ds=[feature_density(fs[i]) for i in ii]; area=float(ar.sum())
                b=dict(forest_area_ha=area,R_area_weighted=float(sum(a*d['R'] for a,d in zip(ar,ds))/area) if area else 0,
                       CF_measured_area_weighted=float(sum(a*d['CF_measured'] for a,d in zip(ar,ds))/area) if area else .5)
            area=b['forest_area_ha']; root=b['R_area_weighted']; cf=b['CF_measured_area_weighted']
            r=dict(scope=scope,forest_type=typ,n_points=len(qs),n_unique_shots=len({q['shot_number'] for q in qs}),
                   n_granules=len({q['granule'] for q in qs}),years=';'.join(sorted({q['date'][:4] for q in qs})),
                   forest_area_ha=area,R_area_weighted=root,CF_measured_area_weighted=cf,
                   estimator='점 단순평균×2025 숲면적; 궤도·연도 중복가중 보정 없음')
            if qs:
                ag=float(np.mean([float(q['agbd']) for q in qs])); c=ag*(1+root)*.5
                r.update(AGB_t_ha=ag,AGB_total_t=ag*area,root_t_ha=ag*root,root_total_t=ag*root*area,
                    carbon_tC_ha=c,carbon_total_tC=c*area,CO2_total_tCO2=c*area*CO2,CO2_tCO2_ha=c*CO2,
                    carbon_measured_CF_total_tC=ag*(1+root)*cf*area,
                    AGB_SD_t_ha=float(np.std([float(q['agbd']) for q in qs],ddof=1)) if len(qs)>1 else '',
                    V21_reference_carbon_tC=c*area/.74)
            summaries.append(r)
            if typ=='전체': lookup[scope]=r
    save_csv(OUT/'carbon_C_summary.csv',summaries)
    save_csv(OUT/'carbon_C_points_in_forest.csv',[{k:v for k,v in q.items() if k!='geometry'} for q in inside])
    check('C shot_number 중복 없음',len({r['shot_number'] for r in inside})==len(inside),len(inside))
    return inside,summaries,lookup

def compare_locations(a_plots,layers,b_records,qs,totals):
    trees={year:shapely.STRtree(np.array([f['geometry'] for f in fs],dtype=object)) for year,fs in layers.items()}
    qtree=shapely.STRtree(np.array([q['geometry'] for q in qs],dtype=object)); geod=Geod(ellps='WGS84')
    result=[]
    for a in a_plots:
        lon=float(a['경도']); lat=float(a['위도']); point=Point(*TO_XY.transform(lon,lat))
        row=dict(plot_id=a['plot_id'],zone=a['zone'],lon=lon,lat=lat,A_AGB_t_ha=a['AGB_t_ha'],A_carbon_tC_ha=a['carbon_tC_ha'])
        local_R=None
        for year,fs in layers.items():
            hits=[int(i) for i in trees[year].query(point,predicate='intersects') if fs[i]['attrs']['FRTP_CD'] in TYPES]
            if hits:
                i=min(hits); d=b_records[year][i]
                row.update({f'B{year}_polygon_id':d['polygon_id'],f'B{year}_carbon_tC_ha':d['carbon_tC_ha'],
                            f'B{year}_AGB_t_ha':d['AGB_t_ha'],f'B{year}_species':d['species']})
                if year==2025: local_R=d['R']
            else:
                forestids=[i for i,f in enumerate(fs) if f['attrs']['FRTP_CD'] in TYPES]
                foresttree=shapely.STRtree(np.array([fs[i]['geometry'] for i in forestids],dtype=object))
                nearest=forestids[int(foresttree.nearest(point))]
                row[f'B{year}_nearest_forest_m']=point.distance(fs[nearest]['geometry'])
        near=[int(i) for i in qtree.query(point,predicate='dwithin',distance=101)
              if geod.inv(lon,lat,float(qs[i]['lon']),float(qs[i]['lat']))[2]<=100]
        row['C_n_points_100m']=len(near)
        row['C_shots_100m']=';'.join(qs[i]['shot_number'] for i in near)
        row['C_years_100m']=';'.join(sorted({qs[i]['date'][:4] for i in near}))
        if near:
            ag=float(np.mean([float(qs[i]['agbd']) for i in near]))
            rr=local_R if local_R is not None else totals[2025,G]['R_area_weighted']
            row.update(C_AGB_t_ha=ag,C_R=rr,C_R_basis='2025 해당폴리곤' if local_R is not None else '금정산 면적가중R 대체',
                       C_carbon_tC_ha=ag*(1+rr)*.5)
        for label in ('B2013','B2025','C'):
            value=row.get(label+'_carbon_tC_ha')
            if value is not None:
                row[label+'_minus_A_tC_ha']=value-a['carbon_tC_ha']
                row[label+'_relative_difference_pct']=100*(value/a['carbon_tC_ha']-1)
        result.append(row)
    save_csv(OUT/'carbon_compare_plots.csv',result)
    metrics=[]
    for label in ('B2013','B2025','C'):
        for subset in ('해당방법 가능한 전체','B2013와 C 공통 조사구'):
            used=[r for r in result if label+'_carbon_tC_ha' in r and (subset=='해당방법 가능한 전체' or all(k in r for k in ('B2013_carbon_tC_ha','C_carbon_tC_ha')))]
            a=np.array([r['A_carbon_tC_ha'] for r in used]); b=np.array([r[label+'_carbon_tC_ha'] for r in used]); delta=b-a
            metrics.append(dict(comparison='A vs '+label,subset=subset,n=len(a),Pearson_r=float(np.corrcoef(a,b)[0,1]) if len(a)>1 else '',
                 mean_A_tC_ha=float(a.mean()) if len(a) else '',mean_other_tC_ha=float(b.mean()) if len(a) else '',
                 mean_difference_tC_ha=float(delta.mean()) if len(a) else '',
                 ratio_of_means_difference_pct=100*float(b.mean()/a.mean()-1) if len(a) else '',
                 mean_relative_difference_pct=float(np.mean(delta/a)*100) if len(a) else '',
                 MAE_tC_ha=float(np.mean(abs(delta))) if len(a) else '',RMSE_tC_ha=float(np.sqrt(np.mean(delta**2))) if len(a) else ''))
    save_csv(OUT/'carbon_comparison_statistics.csv',metrics)
    return result,metrics

def compare_and_sensitivity(a,btot,ctot,layers,areas):
    comparisons=[]; sensitivity=[]
    for scope in (G,U,W):
        sample=a if scope!=U else [p for p in a if p['plot_id'] in [f'남측-{i:02}' for i in range(1,6)]]
        b=btot[2025,scope]; c=ctot[scope]; area=b['forest_area_ha']
        for method in ('A','B','C'):
            if method=='A':
                ag=float(np.mean([p['AGB_t_ha'] for p in sample])); carbon=float(np.mean([p['carbon_tC_ha'] for p in sample])); total=carbon*area
                row=dict(scope=scope,method=method,AGB_t_ha=ag,carbon_tC_ha=carbon,carbon_total_tC=total,
                         sample_n=len(sample),sample_area_ha=len(sample)*.04,sample_carbon_total_tC=sum(p['carbon_total_tC'] for p in sample),
                         total_basis='표본평균×B2025숲면적 참고외삽(대표성 미확인)',
                         sample_basis='장전동 인접 남측01~05; 캠퍼스 내부 표본 아님' if scope==U else '금정산 일대 야장50 전체(등산로 주변)')
                scenarios={'기준_D_CF0.5':carbon,'수종별실측CF_없는수종0.5':float(np.mean([p['carbon_measured_CF_tC_ha'] for p in sample])),
                           'D와H식_CF0.5':float(np.mean([p['carbon_DH_tC_ha'] for p in sample]))}
            elif method=='B':
                row=dict(scope=scope,method=method,AGB_t_ha=b['AGB_t_ha'],carbon_tC_ha=b['carbon_tC_ha'],carbon_total_tC=b['carbon_total_tC'],
                         total_basis='2025 폴리곤 면적합×지역·임상·영급2020 축적×국가계수')
                scenarios={'기준_지자체축적_CF0.5':b['carbon_tC_ha'],'수종별실측CF_없는수종0.5':b['carbon_measured_CF_total_tC']/area}
                for region in ('부산광역시','경상남도','전국'):
                    total=sum(ar*feature_density(f,region)['carbon_tC_ha'] for ar,f in zip(areas[2025][scope],layers[2025]))
                    scenarios['전구역_'+region+'_축적']=total/area
            else:
                row=dict(scope=scope,method=method,AGB_t_ha=c.get('AGB_t_ha',''),carbon_tC_ha=c.get('carbon_tC_ha',''),carbon_total_tC=c.get('carbon_total_tC',''),
                         sample_n=c['n_points'],total_basis='2025숲 안 품질GEDI V3 점평균×숲면적; R 면적가중')
                scenarios={'기준_V3_CF0.5':c['carbon_tC_ha'],'수종별실측CF_없는수종0.5':c['carbon_measured_CF_total_tC']/area,
                           'V2.1참고환산_나누기0.74':c['carbon_tC_ha']/.74} if c['n_points'] else {}
            row['forest_area_ha']=area
            row['CO2_total_tCO2']=row['carbon_total_tC']*CO2 if row['carbon_total_tC']!='' else ''
            row['CO2_tCO2_ha']=row['carbon_tC_ha']*CO2 if row['carbon_tC_ha']!='' else ''
            row['AGB_total_t']=row['AGB_t_ha']*area if row['AGB_t_ha']!='' else ''
            if scenarios:
                base=row['carbon_tC_ha']
                for scenario,value in scenarios.items():
                    sensitivity.append(dict(scope=scope,method=method,scenario=scenario,forest_area_ha=area,
                        carbon_tC_ha=value,carbon_total_tC=value*area,CO2_total_tCO2=value*area*CO2,
                        difference_tC=(value-base)*area,difference_pct=100*(value/base-1),
                        total_basis=row['total_basis'],range_interpretation='1요인씩 바꾼 시나리오 범위; 신뢰구간 아님'))
                row.update(sensitivity_min_tC_ha=min(scenarios.values()),sensitivity_max_tC_ha=max(scenarios.values()),
                           sensitivity_min_total_tC=min(scenarios.values())*area,sensitivity_max_total_tC=max(scenarios.values())*area,
                           sensitivity_min_total_tCO2=min(scenarios.values())*area*CO2,sensitivity_max_total_tCO2=max(scenarios.values())*area*CO2)
            comparisons.append(row)
    save_csv(OUT/'carbon_compare_methods.csv',comparisons); save_csv(OUT/'carbon_sensitivity.csv',sensitivity)
    return comparisons,sensitivity

def patch_polygon(p):
    vertices=[]; codes=[]
    # Orient holes opposite to shells so they remain unfilled.
    p=shapely.geometry.polygon.orient(p,sign=1)
    for ring in [p.exterior]+list(p.interiors):
        a=np.asarray(ring.coords); vertices.extend(a); codes.extend([MplPath.MOVETO]+[MplPath.LINETO]*(len(a)-2)+[MplPath.CLOSEPOLY])
    return PathPatch(MplPath(vertices,codes))

def figures(fs,records,boundary_geom,plots,location,comparison):
    plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False,'font.size':10})
    fig,ax=plt.subplots(figsize=(10,12))
    patches=[patch_polygon(f['geometry']) for f in fs]; values=np.array([r['carbon_tC_ha'] for r in records])
    collection=PatchCollection(patches,cmap='YlGn',edgecolor='none',rasterized=True)
    collection.set_array(values); collection.set_clim(0,float(np.percentile(values[values>0],99)))
    ax.add_collection(collection)
    for poly in parts(boundary_geom):
        x,y=poly.exterior.xy; ax.plot(x,y,color='#b34241',lw=.65,alpha=.85)
    for poly in parts(CAMPUS_GEOM):
        x,y=poly.exterior.xy; ax.plot(x,y,color='#6742ad',lw=2)
    xy=np.array([TO_XY.transform(float(p['경도']),float(p['위도'])) for p in plots]); ax.scatter(xy[:,0],xy[:,1],s=14,c='black',edgecolors='white',linewidths=.35,zorder=5)
    for p,(x,y) in zip(plots,xy): ax.annotate(p['plot_id'],(x,y),xytext=(3,2),textcoords='offset points',fontsize=5.3)
    ax.plot([],[],color='#b34241',label='금정산 근사 숲 경계 (20m 연결·타 산 제외)'); ax.plot([],[],color='#6742ad',label='부산대 캠퍼스 범위'); ax.scatter([],[],s=15,c='black',label='2009~2010 야장 조사구50')
    ax.legend(loc='upper right',fontsize=8); ax.set_aspect('equal'); ax.set_xlim(CLIP.bounds[0],CLIP.bounds[2]); ax.set_ylim(CLIP.bounds[1],CLIP.bounds[3]); ax.ticklabel_format(style='plain',useOffset=False)
    ax.set_xlabel('EPSG:5179 동향좌표 (m)'); ax.set_ylabel('북향좌표 (m)'); ax.set_title('방법 B · 2025 임상도 기반 생체 탄소밀도\n국가축적2020 + WD·BEF·R, CF=0.5')
    fig.colorbar(collection,ax=ax,label='지상부+뿌리 탄소 (tC/ha)',fraction=.035,pad=.025,extend='max')
    fig.text(.1,.025,'공식 공원경계 아님 · 지도 미피복은 미산정 · 죽림/무립목지 0 · 토양/고사목/낙엽 제외',fontsize=8)
    fig.savefig(OUT/'carbon_map_B_2025.png',dpi=190,bbox_inches='tight'); plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(15,5.4))
    for ax,label,title in zip(axes,('B2013','B2025','C'),('A vs B (2013 임상도)','A vs B (2025 임상도)','A vs C (GEDI 100m 이내)')):
        used=[r for r in location if label+'_carbon_tC_ha' in r]; a=np.array([r['A_carbon_tC_ha'] for r in used]); b=np.array([r[label+'_carbon_tC_ha'] for r in used])
        ax.scatter(a,b,s=30,alpha=.75,color='#2b8270' if label!='C' else '#ce8742'); limit=max(a.max(),b.max())*1.08
        ax.plot([0,limit],[0,limit],'--',color='.45',lw=1); ax.set_xlim(0,limit); ax.set_ylim(0,limit); ax.set_aspect('equal'); ax.grid(alpha=.2)
        ax.set_title(f'{title}\nn={len(a)}, r={np.corrcoef(a,b)[0,1]:.3f}, MAE={np.mean(abs(b-a)):.1f}'); ax.set_xlabel('A 야장 탄소밀도 (tC/ha)'); ax.set_ylabel('B/C 탄소밀도 (tC/ha)')
    fig.suptitle('같은 위치 비교 · A:2009~2010 / B:지도2013·2025+축적2020 / C:2019~2025',fontsize=12)
    fig.tight_layout(); fig.savefig(OUT/'carbon_compare_plots.png',dpi=180,bbox_inches='tight'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,6))
    palette=['#657ca5','#2c8b78','#d59a45']
    for j,method in enumerate(('A','B','C')):
        rs=[next(r for r in comparison if r['scope']==s and r['method']==method) for s in (G,U)]
        yy=np.array([r['carbon_tC_ha'] for r in rs]); lo=np.array([r['sensitivity_min_tC_ha'] for r in rs]); hi=np.array([r['sensitivity_max_tC_ha'] for r in rs]); xx=np.arange(2)+(j-1)*.24
        ax.bar(xx,yy,.23,color=palette[j],label={'A':'A 야장 표본평균 (캠퍼스는 인접5개 대용)','B':'B 임상도+국가통계','C':'C GEDI V3'}[method],yerr=np.vstack([yy-lo,hi-yy]),capsize=4)
        for x,y in zip(xx,yy): ax.text(x,y-7,f'{y:.1f}',ha='center',va='top',fontsize=10,color='white')
    ax.set_xticks([0,1],[G,U]); ax.set_ylabel('지상부+뿌리 탄소밀도 (tC/ha)'); ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True)
    ax.set_title('세 방법의 탄소밀도와 가정별 민감도'); ax.legend(fontsize=9,loc='upper right'); ax.set_ylim(0,max(r['sensitivity_max_tC_ha'] for r in comparison if r['scope'] in (G,U))*1.25)
    fig.text(.1,.015,'오차막대 = CF·식·축적지역·V2.1 참고환산을 각각 바꾼 범위 (통계적 신뢰구간 아님)\nA는 등산로 주변 표본; 캠퍼스 내부 조사구 없음. B/C와 측정 시점·공간 지지면적이 다름.',fontsize=8)
    fig.savefig(OUT/'carbon_compare_methods.png',dpi=190,bbox_inches='tight'); plt.close(fig)

def report(a_trees,a,agroups,species,btot,bsummary,csum,ctot,compare,locations,metrics,bmeta,validation):
    fmt=lambda v:f'{v:,.1f}'
    lines=['# 금정산·부산대 산림 탄소저장량 비교','',
      '생체 지상부+뿌리 탄소만 계산했다. 토양·고사목·낙엽층·지피식생은 포함하지 않는다. CF 기본값0.5, CO₂ 환산44/12. 결과는 흡수량이나 배출량이 아닌 저장량 추정이다.',
      '', '## 구역별 핵심 결과','',
      '| 구역·방법 | AGB(t/ha) | 탄소(tC/ha) | 탄소총량(tC) | CO₂(tCO₂) | 가정별 탄소총량 범위(tC) |',
      '|---|---:|---:|---:|---:|---:|']
    for s in (G,U):
        for m in ('A','B','C'):
            r=next(r for r in compare if r['scope']==s and r['method']==m)
            label=m+(' 참고외삽*' if m=='A' else '')
            lines.append(f'| {s} {label} | {fmt(r["AGB_t_ha"])} | {fmt(r["carbon_tC_ha"])} | {fmt(r["carbon_total_tC"])} | {fmt(r["CO2_total_tCO2"])} | {fmt(r["sensitivity_min_total_tC"])}–{fmt(r["sensitivity_max_total_tC"])} |')
    lines+=['','*A 총량은 표본 평균을 B2025의 숲 면적에 단순 곱한 참고 외삽이다. 산 전체/캠퍼스의 검증된 재고량으로 해석하지 않는다. '
            '금정산은 야장50개 평균, 캠퍼스는 캠퍼스 내부에 조사구가 없어 장전동 인접 남측01~05를 대용했다. '
            '가정별 범위는 한 번에 한 요인씩 바꾼 시나리오의 최저~최고이며 신뢰구간이 아니다.', '']
    lines+=['| 구역·방법 | 민감도 범위(tC/ha) | 민감도 범위(tCO₂ 총량) |','|---|---:|---:|']
    for s in (G,U):
        for m in ('A','B','C'):
            r=next(r for r in compare if r['scope']==s and r['method']==m)
            lines.append(f'| {s} {m} | {r["sensitivity_min_tC_ha"]:.1f}–{r["sensitivity_max_tC_ha"]:.1f} | {fmt(r["sensitivity_min_total_tCO2"])}–{fmt(r["sensitivity_max_total_tCO2"])} |')
    lines+=['']
    for s in (G,U,W):
        b=btot[2025,s]; c=ctot[s]
        lines.append(f'- {s}: 2025 산정 숲(침엽/활엽/혼효) {fmt(b["forest_area_ha"])}ha, 죽림0 처리 {fmt(b["bamboo_zero_area_ha"])}ha, GEDI {c["n_points"]:,}점/{c["n_granules"]}개 파일, 관측연도 {c["years"]}, 면적가중 R={b["R_area_weighted"]:.4f}.')
    lines += [f'- 캠퍼스 사각형 전체는 {CAMPUS_GEOM.area/10000:.2f}ha이나, 지도에 있는 산정 숲만 {btot[2025,U]["forest_area_ha"]:.2f}ha이다. 미피복 공간을 탄소0으로 간주하지 않았다.',
              f'- V2.1 참고환산(V3÷0.74): 금정산 {fmt(ctot[G]["carbon_total_tC"]/.74)}tC ({ctot[G]["carbon_tC_ha"]/.74:.1f}tC/ha), 캠퍼스 {fmt(ctot[U]["carbon_total_tC"]/.74)}tC ({ctot[U]["carbon_tC_ha"]/.74:.1f}tC/ha). 이 비율은 금정산에서 검증되지 않았다.',
              '', '## 세 방법 차이: 해석 요약','']
    bg=btot[2025,G]; cg=ctot[G]
    overall=next(r for r in agroups if r['group']=='전체50' and r['metric']=='carbon_tC_ha')
    common=[r for r in metrics if r['subset']=='B2013와 C 공통 조사구']
    mb=next(r for r in common if r['comparison']=='A vs B2013'); mc=next(r for r in common if r['comparison']=='A vs C')
    closest='B(2013)' if mb['MAE_tC_ha']<mc['MAE_tC_ha'] else 'C(GEDI)'
    lines += [f'금정산 기준 B는 {bg["carbon_tC_ha"]:.1f}, C는 {cg["carbon_tC_ha"]:.1f}tC/ha로 C가 B보다 {100*(cg["carbon_tC_ha"]/bg["carbon_tC_ha"]-1):+.1f}% 다르다. 야장50 평균은 {overall["mean"]:.1f}tC/ha이다.',
              'B는 임상·영급별 국가 평균 축적을 각 폴리곤에 적용하므로 개별 조사구의 직경 분포와 본수를 재현하지 못한다. A는 등산로 주변 표본이고 C는 위성 궤도 표본이어서 공간 대표성도 다르다.',
              f'B2013와 C를 모두 비교할 수 있는 동일 {mb["n"]}개 조사구에서 MAE는 B {mb["MAE_tC_ha"]:.1f}, C {mc["MAE_tC_ha"]:.1f}tC/ha로 {closest}가 A에 더 가깝다. 정확도 검증으로 일반화할 수 없다.',
              'A는2009~2010, B는 지도2013/2025에 축적2020을 결합했고 C는2019~2025이다. 최대10~15년가량의 시점 차이, 400㎡ 조사구와 약25m GEDI 발자국, 100m 매칭 반경의 차이를 감안해야 한다.',
              '민감도 중 GEDI ÷0.74 환산은 홍천의 경험적 비율을 옮긴 참고 시나리오이며 금정산에서 검증된 보정계수가 아니다. 수종 대용과 경계 정의의 불확실성은 이 범위보다 클 수 있다.',
              '', '## 금정산 주요 민감도','',
              '| 방법 | 변경 가정 | 탄소총량(tC) | 기준 대비(%) |','|---|---|---:|---:|']
    for r in read_csv(OUT/'carbon_sensitivity.csv'):
        if r['scope']==G and not r['scenario'].startswith('기준'):
            lines.append(f'| {r["method"]} | {r["scenario"]} | {fmt(float(r["carbon_total_tC"]))} | {float(r["difference_pct"]):+.2f} |')
    lines += ['', 'A 민감도 총량도 표본평균의 참고 외삽이다. B의 지역 변경은 원래 부산/양산을 구분하던 축적을 전 구역에 같은 지역 값으로 대체한 것이다.',
              '', '## A 현장 계산과 수종 대용','',
              f'- 원자료 {len(a_trees)}그루, DBH<6cm {sum(float(r["흉고직경_cm"])<6 for r in a_trees)}그루, 등나무 {sum(r["수종"]=="등나무" for r in a_trees)}그루 제외; 실제 계산 {sum(r["included"] for r in a_trees)}그루.',
              f'- 조사구50 ×0.04ha=2ha에서 계산한 생체 탄소 합계 {sum(p["carbon_total_tC"] for p in a):.3f}tC ({sum(p["CO2_total_tCO2"] for p in a):.3f}tCO₂).',
              f'- 탄소밀도 평균 {overall["mean"]:.2f}, 중앙값 {overall["median"]:.2f}, 범위 {overall["min"]:.2f}–{overall["max"]:.2f}tC/ha. 지상부 바이오매스 평균 {np.mean([p["AGB_t_ha"] for p in a]):.2f}t/ha.', '',
              '| 표본구역 | 조사구수 | 탄소평균(tC/ha) | D+H식 평균(tC/ha) |', '|---|---:|---:|---:|']
    for label in ('중앙','동측','서측','남측','북측','장전동 인접 남측01~05'):
        r=next(r for r in agroups if r['group']==label and r['metric']=='carbon_tC_ha'); h=next(r for r in agroups if r['group']==label and r['metric']=='carbon_DH_tC_ha')
        lines.append(f'| {label} | {r["n_plots"]} | {r["mean"]:.2f} | {h["mean"]:.2f} |')
    used=[r for r in a_trees if r['included']]; total=sum(r['carbon_kg'] for r in used)
    other=[r for r in used if r['mapping']=='기타 활엽수→낙엽활엽수 평균']
    flagged=[r for r in used if r['out_of_range_targets']]
    largepine=[r for r in used if r['수종']=='소나무' and float(r['흉고직경_cm'])>40]
    lines += ['',f'- 기타 활엽수 평균 대용 {len(other)}그루, A 탄소량의 {sum(r["carbon_kg"] for r in other)/total*100:.2f}%. 참나무4종 평균 대용은 별도로 {sum(r["mapping"]=="참나무4종 평균" for r in used)}그루이다.',
              f'- 소나무>40cm {len(largepine)}그루는 강원식으로 전환했다. 전환된 나무의 탄소 비중은 {sum(r["carbon_kg"] for r in largepine)/total*100:.2f}%. 적용범위 밖 계산 {len(flagged)}그루, 탄소 비중 {sum(r["carbon_kg"] for r in flagged)/total*100:.2f}%.',
              '- 적용범위 flag는 평균 대용의 구성 수종 중 하나라도 범위를 벗어나면 켰다. 나무별 파일에 대상 수종을 기록했다.',
              '- 원문 의심 표기 개어서나무·대팻집·고뢰쇠나무는 임의로 수종을 정정하지 않고 기타 활엽수 평균을 적용했다.',
              '- 수종별 실측 CF는 목재 실측값을 전체 생체 부위에 공통 적용한 민감도다. 실측 CF가 없는 계수 수종은0.5, 평균 대용은 구성수종 CF의 같은 가중평균이다.',
              '', '## B 수종 대용과 영급 대체의 비중','',
              '| 구역 | 기타활엽 평균 대용 면적비(%) | 그 대용의 탄소비(%) | 1영급→2영급 면적(ha) |','|---|---:|---:|---:|']
    shares=read_csv(OUT/'carbon_B_mapping_contributions.csv')
    for scope in (G,U,W):
        r=next(r for r in shares if r['year']=='2025' and r['scope']==scope and r['mapping']=='기타 활엽수→낙엽활엽수 평균')
        lines.append(f'| {scope} | {float(r["area_share_pct"]):.2f} | {float(r["carbon_share_pct"]):.2f} | {btot[2025,scope]["age1_substitution_area_ha"]:.3f} |')
    lines += ['','참나무4종 평균과 혼효림 그룹평균은 위 기타활엽 대용 비중에서 제외했다. 전체 수종 대응별 면적·탄소 비중은 carbon_B_mapping_contributions.csv에 따로 저장했다.',
              f'A에서 실측CF가 없는 계수수종의 대용 가중치로 배분한 탄소 비중은 {100*sum(r["carbon_kg"]*r["CF_fallback_weight"] for r in used)/total:.2f}%이며, 이 부분의 민감도CF는0.5를 유지한다.',
              '', '## 같은 위치 비교','',
              '| 비교 | 표본 | n | 상관계수 r | 평균차이(상대방법−A,tC/ha) | 평균의 비율 차이(%) | MAE(tC/ha) | RMSE(tC/ha) |',
              '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in metrics:
        lines.append(f'| {r["comparison"]} | {r["subset"]} | {r["n"]} | {r["Pearson_r"]:.3f} | {r["mean_difference_tC_ha"]:.2f} | {r["ratio_of_means_difference_pct"]:.1f} | {r["MAE_tC_ha"]:.2f} | {r["RMSE_tC_ha"]:.2f} |')
    lines += ['',f'- 100m 이내 숲 안 GEDI가 있는 조사구 {sum(r["C_n_points_100m"]>0 for r in locations)}곳. 매칭은 WGS84 타원체 거리로 재검사했다. 100m 근처 점평균은 조사구와 동일한 지표면 측정이 아니다.',
              '- 표의 차이율은 (상대방법 평균/A 평균−1)이다. 개별 조사구 차이율 평균도CSV에 저장했으나 A가 작은 조사구에서 크게 튈 수 있어 구분해서 해석해야 한다.',
              '- 조사구 중심이 임상도 숲 밖이면 B 값은 빈칸이며 가장 가까운 숲 거리만 기록했다. C 뿌리비는2025 중심폴리곤 R을 사용하고, 중심이 밖이면 금정산 면적가중R로 대체하고 표시했다.',
              '', '## 2013 vs 2025 (같은 통계·계수)','', '| 구역 | 2013(tC) | 2025(tC) | 차이(%) |','|---|---:|---:|---:|']
    for scope in (G,U,W):
        old=btot[2013,scope]['carbon_total_tC']; new=btot[2025,scope]['carbon_total_tC']
        lines.append(f'| {scope} | {fmt(old)} | {fmt(new)} | {100*(new/old-1):+.1f} |')
    lines += ['', '이 차이는 관측된 탄소 증가량이 아니다. 지도 영급·수종·피복·갱신 차이, 1영급 대체, NFI2020 평균 적용이 함께 작용한다. 금정산은2025 연결 숲 경계를2013에도 고정했으므로, 그 밖에서2013 이후 사라진 숲까지 포함하는 산림변화 분석은 아니다.',
              '', '## 경계·주요 가정과 한계','',
              f'- 0m 주능선 씨앗 숲 {bmeta["strict_seed_area_ha"]:.2f}ha → 타 산 제외+20m 연결 {bmeta["selected_area_including_bamboo_ha"]:.2f}ha(죽림 포함). 타 산 제외 없이20m만 적용하면 {bmeta["bare_20m_area_before_other_mountain_exclusion_ha"]:.2f}ha로 넓어진다. 이 차이는 지형 경계의 근사 불확실성이다.',
              '- 이전 참고선(고당봉~중앙 야장열~남문 일대 근사점)의100m 완충을 씨앗으로 사용했다. 남문 좌표는 실측 성문 위치가 아니다. 공식 공원경계나 DEM 능선 경계가 아니다.',
              '- 먼저0m 연결 성분과 산림청 등산로 교차길이를 계산하고, 최다 산명이 금정산/인재개발원화명수목원이 아닌 성분을 제외한 뒤20m 연결했다. 산명이 없는 작은 성분은 연결되면 포함될 수 있다. 제외 산명·면적은 carbon_boundary_excluded_components.csv에 모두 남겼다.',
              '- 산정 숲은 FRTP_CD1/2/3; 죽림4·무립목지·제지는 사용자 지시대로0이다. C도 같은 산정 숲 안 점만 사용해 죽림에 AGBD를 외삽하지 않았다. 면적 분모는 산정 숲이다.',
              '- B: 부산 금정구·북구는 부산광역시, 양산시는 경상남도 축적. 1영급은2영급, 6~9영급은6영급 이상 통계. WD·BEF·R은 각각 평균한 뒤 곱했다. 혼효림은 침엽/낙엽활엽 그룹 평균을1:1로 합쳤다.',
              '- B의 소나무는 개체 직경을 알 수 없으므로 중부지방 계수를 일괄 적용했다. 대경목 등급만으로 강원식 계수를 배정하지 않았다. 수종별 축적이 아니라 임상별 축적이라는 한계가 있다.',
              '- 기타 활엽수(상록활엽수의 미수록 수종도 포함)는 README 방침대로 낙엽활엽수 평균; 명시된 보고서 수종은 직접 대응. 임상도 백합나무는 보고서 튤립나무와 같은 수종명으로 대응했고, 포플러는 현사시로 단정하지 않았다.',
              '- 곰솔 BEF 불확도112.782%는 원문 계수표 값을 그대로 보존한다. 계수 불확도를 합성한 신뢰구간을 계산하지 않았고, 가정 민감도 범위는 이 불확도를 포함하지 않는다.',
              '- README 가지 검산 예시의0.069×24^1.995는39.6kg가 아니라39.117448kg이다. 산술 예시의 불일치를 기록하고, 제공CSV의 계수는 수정하지 않았다.',
              '- A 조사구의 소경목 제외·덩굴 제외와 수종 평균 대용, 수고식 선택으로 결과가 달라진다. 평균 수종 계수를 비선형식에 먼저 넣지 않고 각 수종 예측 건중량을 평균했다.',
              '- C는 관측점 단순평균으로, 반복궤도/관측연도/수종별 면적 층화 가중을 기본 추정에 적용하지 않았다. 점 수가 많아도 확률표본이 아니므로 작은 표준오차를 정확도의 증거로 제시하지 않았다.',
              '- GEDI 임상별 평균은 별도 저장했다. R과 민감도CF는 B의 각 구역·임상별 면적가중 계수이다. V2.1 ÷0.74는 참고만 가능하다.',
              '- 야장 GPS의 원 측지계 미기재·초 단위 반올림, 임상도 위치오차, 시기 차이와 GEDI 발자국의 부분 겹침이 있다. 재선정 후보 목록은 방문 참고자료이며 탄소 추정 표본을25곳으로 제한하지 않았다.',
              '- B에서 전체 사각형, 금정산, 캠퍼스는 서로 겹치는 집계이며 더하면 중복된다. 미피복을 탄소0으로 취급하지 않았다.',
              '', '## 단위·합계 검산','',f'계산 과정 {len(validation)}개 검증 모두 통과. 상세는 carbon_validation.csv.',
              f'저장 결과와 원래SHP·야장을 별도로 대조한 독립 검증 {len(read_csv(OUT/"carbon_independent_validation.csv"))}개도 모두 통과. 상세는 carbon_independent_validation.csv.',
              '', '## 자료 출처와 재현','',
              '- 제공된 data/reference/README.md 및4개 계수/통계CSV를 원본 그대로 사용했다. nifos2014_allometry.csv의 page는 국립산림과학원2014 보고서 PDF 쪽이다.',
              '- CF: 국립산림과학원2010 보고서 PDF50쪽 표10. 축적: KOSIS DT_136014_10210, 제7차 NFI2020, README의 확인 내역.',
              '- 야장: 서정범2011 동아대 박사논문 부록2,2009~2010 측정. 지도: 제공2013·2025 임상도와 산림청 등산로. GEDI: 로컬 품질통과V3 자료.',
              '- carbon_input_manifest.csv에 입력 SHA256, carbon_assumptions.csv에 대응·가정, carbon_validation.csv에 검산을 기록했다. 원자료를 수정하거나 외부로 업로드하지 않았다.',
              '', '## 생성 파일','']
    for p in sorted(set(OUT.glob('carbon_*'))|{OUT/'carbon_report.md'}):
        lines.append('- '+p.name)
    (OUT/'carbon_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
    check('금지 패키지 미사용',not any(k.split('.')[0] in ('pandas','geopandas','rasterio','pyogrio') for k in sys.modules),'sys.modules')
    for name,ok,detail in [('kg/t 변환',1000/1000==1,'1000kg=1t'),('면적변환',400/10000==.04,'400m2=0.04ha'),
        ('조사구 확장',math.isclose(1000/1000/.04,25),'1000kg/400m2=25t/ha'),
        ('재적·질량·뿌리·탄소',math.isclose(100*2*.5*1.5*(1+.2)*.5,90),'100m3/ha×2ha×WD.5×BEF1.5×(1+R.2)×CF.5=90tC'),
        ('CO2분자량換算',math.isclose(12*CO2,44),'12tC=44tCO2'),
        ('가지식 독립 수치검산',math.isclose(.069*24**1.995,39.117448363572755,rel_tol=1e-12),'.069×24^1.995=39.117448kg; README39.6과 불일치')]: check(name,ok,detail)
    check('계수22종·2식·4부위',len(FACT)==22 and len(EQ)==22*2*4,(len(FACT),len(EQ)))
    for s in FACT:
        for form in ('a*D^b','a*D^b*H^c'):
            check('식완비 '+s+' '+form,all((s,form,p) in EQ for p in ('줄기','가지','잎','뿌리')),'4부위')
    assume('단위','환산','kg/1000=t; m2/10000=ha; Mg=t; tCO2=tC×44/12','제공자료 단위 유지')
    assume('입력검산','README 가지 예시 불일치','39.6 표기 vs 식 계산39.117448kg','CSV계수 .069,1.995를 그대로 적용; README 산술예시를 맞추려고 계수 수정 안함')
    assume('CF','기준',.5,'국내 인벤토리 관행(사용자 지정)')
    assume('CF','실측없는수종',.5,'미측정CF를 추정해 만들지 않음; 그룹평균은 각 수종CF 동일가중','nifos2010_carbon_fraction.csv')
    assume('CF','실측CF 전체부위 적용','목재 CF를 가지·잎·뿌리에도 공통적용','부위별CF 자료 없음; 민감도에만 사용')
    assume('A','조사구',.04,'20m×20m=400m2; 관측2009~2010')
    assume('A','제외','DBH<6cm,등나무','대상 식 적용범위와 사용자 규칙')
    assume('A','직경범위초과','계산 유지+나무별flag','소나무>40은 강원식; 평균대용 구성수종 범위도 검사')
    assume('A','수종오타의심','개어서나무/대팻집/고뢰쇠나무→기타활엽 평균','원문을 임의 정정하지 않음')
    assume('A','평균대용','수종별 예측질량 평균','계수평균을 비선형식에 넣지 않음')
    assume('B','1영급','2영급 축적 대체+flag','원통계1영급 없음')
    assume('B','6영급이상','6영급이상 통계','영급6~9 한 구간')
    assume('B','혼효림','침엽평균과 낙엽활엽평균 각0.5','WD·BEF·R·CF 각각 같은 가중으로 평균; 이후 곱')
    assume('B','0 처리','죽림/무립목지/제지/비산림','사용자 지정; 미피복 공간과 구분')
    assume('B','소나무','중부지방소나무','개체DBH가 없는 지도에는40cm 분기 미적용')
    assume('B','지역','금정구·북구 부산광역시 / 양산시 경상남도','SRC_DIST 사용; 중첩 원본 우선순위는 이전 병합 결과 유지')
    assume('B','백합나무','튤립나무','동일수종 명칭; 포플러는현사시로 단정하지 않음')
    assume('C','AGBD','품질통과+2025 FRTP1/2/3 안; 단순점평균','추가수관피복필터 없음; 죽림0규칙 일관 적용')
    assume('C','R','B2025 구역별 산정숲 면적가중R','AGB×(1+R)×CF')
    assume('C','V2.1 참고','V3/0.74','홍천 경험비율; 금정산 검증된 보정 아님')
    assume('비교','A 총량','표본평균×B2025숲면적 참고외삽','금정산50, 캠퍼스 장전동남측5 대용; 대표성 미확인')
    assume('불확실성','계수 불확도','가정 민감도와 별개; 합성 신뢰구간 계산안함','곰솔BEF_unc112.782% 원값 유지')
    plots=read_csv(ROOT/'data/field/seo2011_plots.csv')
    print('A field allometry',flush=True)
    at,a,agroups,asp=calculate_A(plots)
    layers={year:read_merged(f'forest{year}_geumjeong.shp') for year in (2013,2025)}
    trails=read_merged('trails_geumjeong.shp')
    print('Boundary: named mountain components and 20m graph',flush=True)
    geom,strict,bmeta,selected=boundary(layers[2025],trails,plots)
    scopes={G:geom,U:CAMPUS_GEOM,W:CLIP}
    for scope in scopes.values(): shapely.prepare(scope)
    brec,bsum,btot,areas=calculate_B(layers,scopes,selected)
    print('C forest-screened GEDI',flush=True)
    qs,csum,ctot=calculate_C(layers[2025],scopes,btot,areas)
    print('Location comparisons and sensitivities',flush=True)
    locations,metrics=compare_locations(a,layers,brec,qs,btot)
    comparison,sensitivity=compare_and_sensitivity(a,btot,ctot,layers,areas)
    for year in layers:
        rs=brec[year]
        check(f'B{year} 전체면적',math.isclose(sum(r['area_ha'] for r in rs),btot[year,W]['mapped_area_ha'],rel_tol=1e-10),'전체폴리곤면적=구역합')
        check(f'B{year} 재적합',math.isclose(sum(r['volume_total_m3'] for r in rs),btot[year,W]['volume_total_m3'],rel_tol=1e-10),'m3/ha×ha')
        check(f'B{year} 탄소합',math.isclose(sum(r['carbon_total_tC'] for r in rs),btot[year,W]['carbon_total_tC'],rel_tol=1e-10),'tC/ha×ha')
        for scope in scopes:
            total=btot[year,scope]
            for dim in ('forest_type','species'):
                subtotal=sum(r['carbon_total_tC'] for r in bsum if r['year']==year and r['scope']==scope and r['group_dimension']==dim)
                check(f'B{year} {scope} {dim} 합',math.isclose(subtotal,total['carbon_total_tC'],rel_tol=1e-10),subtotal)
    for scope in scopes:
        c=ctot[scope]
        if c['n_points']:
            check('C '+scope+' 총량',math.isclose(c['carbon_total_tC'],c['AGB_t_ha']*c['forest_area_ha']*(1+c['R_area_weighted'])*.5,rel_tol=1e-12),'Mg/ha×ha×(1+R)×CF')
    check('금지패키지종료확인',not any(k.split('.')[0] in ('pandas','geopandas','rasterio','pyogrio') for k in sys.modules),'최종sys.modules')
    # Record actual inputs, not a mutable directory scan of newly generated files.
    files=[REF/'README.md',*sorted(REF.glob('*.csv')),ROOT/'data/field/seo2011_trees.csv',ROOT/'data/field/seo2011_plots.csv',
           OUT/'gedi_v3_quality_points.csv',OUT/'ridge_connectivity_method.json',OUT/'seo2011_vs_forestmap.csv',OUT/'field_sites_reselected.csv']
    for name in ('forest2013_geumjeong','forest2025_geumjeong','trails_geumjeong'):
        files.extend(MERGED/(name+ext) for ext in ('.shp','.shx','.dbf','.prj','.cpg'))
    save_csv(OUT/'carbon_input_manifest.csv',[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files if p.exists()])
    save_csv(OUT/'carbon_assumptions.csv',assumptions); save_csv(OUT/'carbon_validation.csv',tests)
    save_json(OUT/'carbon_environment.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,
        shapely=shapely.__version__,matplotlib=matplotlib.__version__,forbidden_imports_loaded=False))
    print('Figures',flush=True)
    figures(layers[2025],brec[2025],geom,a,locations,comparison)
    from carbon_independent_checks import run as independent_checks
    independent_checks()
    report(at,a,agroups,asp,btot,bsum,csum,ctot,comparison,locations,metrics,bmeta,tests)
    print(json.dumps(dict(status='complete',validation_tests=len(tests),all_passed=all(t['passed'] for t in tests),boundary=bmeta,
        main_results=[r for r in comparison if r['scope'] in (G,U)],statistics=metrics),ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__': main()
