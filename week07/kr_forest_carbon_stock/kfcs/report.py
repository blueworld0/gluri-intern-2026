"""Auditable outputs and narrowly specified comparisons, never parameter search."""
import copy
import importlib.metadata
import json
import platform
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from .volume.base import file_record

def markdown_table(df):
    if df.empty: return '(자료 없음)'
    df=df.copy().fillna('—')
    def fmt(x):
        return f'{x:,.4f}' if isinstance(x,(float,np.floating)) else str(x).replace('|','/').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |',*['| '+' | '.join(fmt(x) for x in row)+' |' for row in df.itertuples(index=False,name=None)]])

def summary(g, group=None):
    def one(x):
        return {'polygon_count':len(x),'mapped_area_ha':x.area_ha.sum(),'forest_area_ha':x.loc[x.eligible,'area_ha'].sum(),
                'excluded_area_ha':x.loc[~x.eligible,'area_ha'].sum(),'missing_carbon_area_ha':x.loc[x.carbon_tC.isna(),'area_ha'].sum(),
                'unclassified_area_ha':x.loc[x.classification_flag.eq('미분류'),'area_ha'].sum(),
                'volume_m3':x.volume_m3.sum(min_count=1),'carbon_tC':x.carbon_tC.sum(min_count=1),'CO2_tCO2':x.CO2_tCO2.sum(min_count=1)}
    if group is None: return pd.DataFrame([one(g)])
    return pd.DataFrame([{group:key,**one(x)} for key,x in g.groupby(group,dropna=False,sort=True)])

def external_comparison(g,c):
    reference=pd.read_csv(Path(c['reference_dir'])/'external_estimates_shin2023.csv',dtype={'sido_code':str})
    # Explicit full-sido declaration must be accompanied by an official coverage polygon.
    b=c['boundary']; full=set(b.get('full_sido_codes',[])) if b.get('scope')=='full_sido' else set()
    coverage_path=b.get('coverage_path'); coverage=None
    if full:
        if not coverage_path: raise ValueError('full_sido requires official coverage_path for geometric coverage check')
        coverage=gpd.read_file(coverage_path).to_crs(c['spatial']['crs'])
        from .aoi import read_boundary
        aoi=read_boundary(c)
        if aoi is None: raise ValueError('full_sido requires explicit AOI')
    out=[]
    for sc,x in g.groupby('sido_code'):
        ref=reference[reference.sido_code==str(sc)]
        if len(ref)!=1: raise ValueError('External reference sido not unique')
        r=ref.iloc[0]; comparable=False
        if sc in full:
            area=coverage[coverage[b['coverage_code_field']].astype(str)==str(sc)]
            if area.empty: raise ValueError('Coverage polygon missing sido')
            official=area.geometry.union_all()
            comparable=official.difference(aoi).area/official.area<=b['coverage_missing_fraction_tolerance']
        out.append({'sido_name':r.sido_name,'status':'비교 가능' if comparable else '비교 불가(대상 지역이 시도 일부)',
                    'year_reference':int(r.year),'reference_carbon_tC':r.carbon_tC,'reference_volume_m3':r.volume_m3,'reference_area_ha':r.forest_area_ha,
                    'carbon_difference_pct':(x.carbon_tC.sum()/r.carbon_tC-1)*100 if comparable else np.nan})
    return pd.DataFrame(out)

def week4_comparison(g,c,out):
    w=c['validation']['week4']; area=w['fallback_area_ha']; carbon=w['fallback_carbon_tC']; origin='요청문 근삿값'
    if Path(w['summary_path']).exists():
        old=pd.read_csv(w['summary_path']); r=old[(old.year==w['year'])&(old.scope==w['scope'])&(old.group_dimension=='total')]
        if len(r)!=1: raise ValueError('Week4 summary row not unique')
        area=float(r.iloc[0].forest_area_ha); carbon=float(r.iloc[0].carbon_total_tC); origin=w['summary_path']
    df=pd.DataFrame([{'metric':'forest_area_ha','reference':area,'current':g.loc[g.eligible,'area_ha'].sum()},
                     {'metric':'carbon_tC','reference':carbon,'current':g.carbon_tC.sum(min_count=1)}])
    df['difference_pct']=(df.current/df.reference-1)*100
    df['within_tolerance']=df.difference_pct.abs()<=c['validation']['relative_tolerance']*100
    df.to_csv(out/'comparison_week4.csv',index=False,encoding='utf-8-sig')
    message='기준 출처: '+origin
    if not df.within_tolerance.all(): message+='\n1% 초과 원인 후보: 원본 SHP 전지역 처리와 4주차 bbox 처리 차이, 중복 우선순위, 경계 교차의 수치 오차, 수종명·영급 코드, 원본 재적표 판 차이. 구획 비교 CSV에서 차이 집중 구획을 확인함.'
    # Week 4 polygon IDs were reassigned after bbox clipping; compare by source geometry
    # membership/attribute grouping instead of pretending those IDs identify raw records.
    if Path(w['polygon_path']).exists():
        old=pd.read_csv(w['polygon_path']); old=old[old.geumjeong_area_ha>0].copy()
        old['FRTP_CD']=old.forest_code.astype(str); old['species_group']=old.species
        old['age_label']=old.age_class
        keys=['FRTP_CD','species_group','age_label']
        a=old.groupby(keys)[['geumjeong_area_ha','geumjeong_carbon_tC']].sum()
        b=g.groupby(keys)[['area_ha','carbon_tC']].sum()
        a.join(b,how='outer').to_csv(out/'comparison_week4_attribute_groups.csv',encoding='utf-8-sig')
        # Each legacy density is constant within the above keys. Area-overlap matching
        # against optional week4 merged GPKG/SHP is supported without ID assumptions.
        path=w.get('geometry_path')
        if path and Path(path).exists():
            legacy=gpd.read_file(path,encoding=c['inputs']['encoding']).to_crs(g.crs)
            legacy['polygon_id']=legacy[w.get('geometry_id_field','FID_U')].astype(int)
            legacy=legacy.merge(old[['polygon_id','geumjeong_area_ha','geumjeong_carbon_tC']],on='polygon_id')
            pairs=gpd.overlay(g[['polygon_id','carbon_density_tC_ha','geometry']],legacy[['polygon_id','geometry']],how='intersection',keep_geom_type=True)
            pairs['overlap_ha']=pairs.area/c['units']['m2_per_ha']
            pairs['current_carbon_tC']=pairs.overlap_ha*pairs.carbon_density_tC_ha
            pc=pairs.groupby('polygon_id_2')[['overlap_ha','current_carbon_tC']].sum().join(old.set_index('polygon_id')[['geumjeong_area_ha','geumjeong_carbon_tC']],how='outer')
            pc['area_difference_pct']=(pc.overlap_ha/pc.geumjeong_area_ha-1)*100
            pc['carbon_difference_pct']=(pc.current_carbon_tC/pc.geumjeong_carbon_tC-1)*100
            pc.to_csv(out/'comparison_week4_polygons.csv',encoding='utf-8-sig')
        else:
            message+='\n4주차 CSV의 polygon_id는 bbox·중복 처리 후 재부여되어 원본 ID와 같지 않음. 구획별 검산에는 validation.week4.geometry_path의 4주차 merged SHP/GPKG가 필요함. 현재는 같은 임상·수종·영급 그룹별로 비교함.'
    return df,message

def admin_comparison(g,c,out):
    from .volume.tables import administrative_reference
    ref=administrative_reference(c); rows=[]
    for r in ref.itertuples():
        if r.group=='전체': x=g
        elif r.group=='무립목지': x=g[g.FRTP_CD.isin(c['classification']['unstocked_codes'])]
        elif r.group=='죽림': x=g[g.FRTP_CD.isin(c['classification']['bamboo_codes'])]
        elif r.dimension=='forest_type': x=g[g.forest_type==r.group]
        else: x=g[g.eligible & g.age_class.clip(upper=c['classification']['age_max']).eq(float(r.group))]
        area=x.area_ha.sum(); volume=x.volume_m3.sum(min_count=1) if len(x) else 0.
        rows.append({'dimension':r.dimension,'group':r.group,'area_ha':area,'reference_area_ha':r.reference_area_ha,
                     'area_difference_pct':(area/r.reference_area_ha-1)*100 if r.reference_area_ha else np.nan,
                     'volume_m3':volume,'reference_volume_m3':r.reference_volume_m3,
                     'volume_difference_pct':(volume/r.reference_volume_m3-1)*100 if pd.notna(r.reference_volume_m3) and r.reference_volume_m3 else np.nan})
    df=pd.DataFrame(rows); df.to_csv(out/'comparison_administrative.csv',index=False,encoding='utf-8-sig')
    return df,'차이 원인 후보: EPSG:5179 도형 면적과 공식 산림면적의 산출 방식 차이, 임상도 갱신년도와 통계 판 차이, 시군구 실측 재적 대신 시도·임상·영급 평균을 사용함. FRTP=0은 무립목지와 비산림이 합쳐져 있어 공식 무립목지와 직접 일치하지 않음.'

def decomposition(g,c,out):
    from .config import load_config
    from .volume import get_provider
    from .factors import FactorSet
    from .carbon import calculate
    spec=c['validation']['decomposition']; base=load_config(spec['baseline_config'],drive_root=c['drive_root'])
    base['reference_dir']=c['reference_dir']
    def total(cc):
        p=get_provider(cc); d=p.provide(g)
        x=g.drop(columns=[k for k in d.columns if k in g],errors='ignore').copy()
        for col in d: x[col]=d[col]
        # Remove previously attached factors before shared recalculation.
        x=x.drop(columns=[k for k in x if k in ['WD','BEF','R','CF','factor_key','factor_rule','factor_note','factor_targets']],errors='ignore')
        x=calculate(FactorSet(cc).attach(x))
        return float(x.carbon_tC.sum(min_count=1))
    previous=total(base); rows=[{'step':'4주차 기준','carbon_tC':previous,'delta_previous_tC':0.,'comparison':'기준'}]
    for step in spec['steps']:
        for dotted,value in step['changes'].items():
            section,key=dotted.split('.'); base[section][key]=value
        value=total(base); rows.append({'step':step['label'],'carbon_tC':value,'delta_previous_tC':value-previous,'comparison':'앞 단계 대비; 순서 의존'}); previous=value
    df=pd.DataFrame(rows); df.to_csv(out/'comparison_decomposition.csv',index=False,encoding='utf-8-sig')
    return df,'지정한 네 단계만 순서대로 실행했음. 첫 세 변화의 합은 4주차→2025 설정 차이이며, 마지막은 2025 설정에서 임상 3종 계수로 변경한 별도 비교임. 상호작용은 순서에 따라 배분됨.'

def uncertainty_comparison(g,c,factors,mc,out):
    from .uncertainty import monte_carlo
    scenarios=[];audits=[]
    configured=c['uncertainty']['overrides'] or str(Path(c['reference_dir'])/'uncertainty_overrides.csv')
    for label,path in [('미적용',None),('적용',configured)]:
        cc=copy.deepcopy(c);cc['uncertainty']['overrides']=path
        result=mc if path==c['uncertainty']['overrides'] else monte_carlo(g,factors,cc)
        scenarios.append({'scenario':label,'overrides':'none' if path is None else Path(path).name,
                          'point_carbon_tC':float(g.carbon_tC.sum()),**{k:result[k] for k in ['median_tC','p05_tC','p95_tC','n','seed']}})
        audits.append({'scenario':label,**result})
    df=pd.DataFrame(scenarios);df.to_csv(out/'comparison_uncertainty_overrides.csv',index=False,encoding='utf-8-sig')
    return df,audits

def inventory_comparison(g,c,out):
    from .config import load_config
    from .factors import FactorSet
    from .carbon import calculate
    path=c['validation']['national_inventory_comparison']['baseline_config']
    base=load_config(path,drive_root=c['drive_root']);base['reference_dir']=c['reference_dir']
    rows=[]
    for label,fs,rule in [('2025 기준','nifos2014_species','factor_average_50_50'),
                          ('계수 세트만','national_inventory_2025','factor_average_50_50'),
                          ('혼효림 규칙만','nifos2014_species','volume_split_50_50'),
                          ('둘 다(국가 인벤토리)','national_inventory_2025','volume_split_50_50')]:
        cc=copy.deepcopy(base);cc['factors'].update(set=fs,mixed_rule=rule,cf_mode='from_factor_set' if fs=='national_inventory_2025' else base['factors']['cf_mode'])
        result=calculate(FactorSet(cc).attach(g))
        rows.append({'scenario':label,'factor_set':fs,'mixed_rule':rule,'forest_area_ha':float(result.loc[result.eligible,'area_ha'].sum()),
                     'volume_m3':float(result.volume_m3.sum()),'carbon_tC':float(result.carbon_tC.sum()),'CO2_tCO2':float(result.CO2_tCO2.sum())})
    df=pd.DataFrame(rows);baseline=df.iloc[0].carbon_tC
    df['delta_vs_2025_tC']=df.carbon_tC-baseline;df['delta_vs_2025_pct']=df.delta_vs_2025_tC/baseline*100
    df.to_csv(out/'comparison_national_inventory.csv',index=False,encoding='utf-8-sig')
    return df

def write_outputs(g,c,provider,factors,mc,audit,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out)
    g.to_file(out/'polygons.gpkg',driver='GPKG',layer='polygons',engine='pyogrio')
    total=summary(g)
    for name,col in [('total',None),('sido','sido_name'),('forest_type','forest_type'),('species','species_group'),('age','age_label')]:
        summary(g,col).to_csv(out/f'summary_{"total" if col is None else "by_"+name}.csv',index=False,encoding='utf-8-sig')
    fallback=g.groupby(['volume_flag'],dropna=False).agg(count=('area_ha','size'),area_ha=('area_ha','sum')).reset_index()
    fallback.to_csv(out/'volume_flags.csv',index=False,encoding='utf-8-sig')
    unclassified=g[g.classification_flag.eq('미분류')].groupby(['FRTP_CD','KOFTR_GROU','AGCLS_CD'],dropna=False).area_ha.agg(['size','sum']).reset_index()
    unclassified.to_csv(out/'unclassified_codes.csv',index=False,encoding='utf-8-sig')
    fig,ax=plt.subplots(figsize=(9,9)); positive=g[g.area_ha>0]
    if not positive.empty:
        positive.plot(column='carbon_density_tC_ha',ax=ax,legend=True,cmap='YlGn',missing_kwds={'color':'lightgrey','label':'Missing'},legend_kwds={'label':'Carbon density (tC/ha)'})
    ax.set_title('Forest carbon stock density'); ax.set_xlabel('Easting (m)');ax.set_ylabel('Northing (m)')
    fig.tight_layout();fig.savefig(out/'carbon_density_map.png',dpi=180);plt.close(fig)
    ext=external_comparison(g,c);ext.to_csv(out/'comparison_external.csv',index=False,encoding='utf-8-sig')
    comparison=pd.DataFrame(); note=''
    mode=c['validation']['mode']
    if mode=='week4': comparison,note=week4_comparison(g,c,out)
    elif mode=='admin': comparison,note=admin_comparison(g,c,out)
    elif mode=='decomposition': comparison,note=decomposition(g,c,out)
    overrides_table=pd.DataFrame(mc.get('overrides_applied',[]));high_table=pd.DataFrame(mc.get('high_uncertainties',[]))
    sensitivity=pd.DataFrame();uncertainty_audits=[];inventory=pd.DataFrame()
    if c['validation'].get('uncertainty_comparison',{}).get('enabled'):
        sensitivity,uncertainty_audits=uncertainty_comparison(g,c,factors,mc,out)
    if c['validation'].get('national_inventory_comparison',{}).get('enabled'):
        inventory=inventory_comparison(g,c,out)
        if mode=='national_inventory':comparison=inventory;note='면적·재적은 같은 상태로 계수 세트만·혼효림 규칙만·둘 다 바꾸었음. 국가 계수의 상록활엽수는 활엽수로 처리하며 죽림은 NE로 0 처리함.'
    overrides_table.to_csv(out/'uncertainty_overrides_applied.csv',index=False,encoding='utf-8-sig')
    high_table.to_csv(out/'uncertainty_high_values.csv',index=False,encoding='utf-8-sig')
    inputs=[]
    for path in c['inputs']['forest_maps']:
        p=Path(path)
        inputs.extend(file_record(s,edition=c['editions']['forest_map']) for s in sorted(p.parent.glob(p.stem+'.*')) if s.is_file())
    boundary_source=str(c['boundary']['path'])
    if boundary_source.startswith('zip://'): boundary_source=boundary_source.removeprefix('zip://').split('!')[0]
    for p in [boundary_source,c['inputs'].get('attribute_pdf'),c.get('_config_path'),c['validation']['week4']['summary_path'],c['validation']['week4']['polygon_path'],c['validation']['week4'].get('geometry_path'),c['boundary'].get('coverage_path')]:
        if p and p!='none' and Path(p).is_file(): inputs.append(file_record(p))
    for key,edition in [('sfs_zip','sfs2025'),('kosis_xlsx','kosis2020')]:
        p=c['inputs'].get(key)
        if p and Path(p).is_file(): inputs.append(file_record(p,edition=c['editions'][edition]))
    inputs.extend(file_record(p,edition=c['editions']['factors']) for p in sorted(Path(c['reference_dir']).glob('*.csv')))
    if c['uncertainty'].get('overrides'):inputs.append(file_record(c['uncertainty']['overrides']))
    packages={}
    for package in ['geopandas','pyogrio','pyproj','shapely','pandas','numpy','rasterio','xlrd','openpyxl','pyyaml','matplotlib','pytest']:
        try: packages[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: packages[package]=None
    prov={'run_time_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),'inputs':inputs,'editions':c['editions'],
          'volume':provider.provenance(),'factor_set':c['factors'],'config':c,'packages':packages,'forest_map_audit':audit,'uncertainty':mc,'uncertainty_sensitivity':uncertainty_audits,
          'implementation':[file_record(p) for p in sorted(Path(__file__).parent.rglob('*.py'))]}
    (out/'provenance.json').write_text(json.dumps(prov,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    unknown=float(g.loc[g.classification_flag.eq('미분류'),'area_ha'].sum())
    excluded=summary(g[~g.eligible],'forest_type')
    text=f'''# 국내 산림 탄소 저장량 결과

{markdown_table(total)}

탄소 풀: 입목 지상부 + 뿌리만 포함함. 고사목·낙엽층·토양은 제외함. 흡수량(tC/년)이 아니라 저장량(tC)임.
결측 면적이 있으면 위 합계는 관측 가능한 구획의 부분 합계이며 전체 저장량 추정치로 해석할 수 없음.

## 가정과 적용 규칙
공급자: {c['volume']['provider']}, 계수: {c['factors']['set']}, CF: {c['factors']['cf_mode']}, 1영급: {c['volume']['age1_rule']}.
수종 대응 프로필: {c['factors']['mapping_profile']}. 혼효림 규칙: {c['factors']['mixed_rule']}.
factor_average_50_50은 WD·BEF·R·CF를 각각 평균한 후 곱함. volume_split_50_50은 재적을 반씩 나눠 침엽·활엽의 WD×BEF×(1+R)×CF를 각각 곱한 뒤 합산함. 분할 방식의 WD·BEF·R·CF 열은 요약값이며 계산에는 carbon_conversion_tC_m3와 conifer_*/broadleaf_* 성분을 사용함.
national_inventory_2025는 수종그룹을 사용하지 않으며, 활엽수(상록활엽 포함)에 활엽수 계수를 쓰고 죽림은 국가 인벤토리 NE로 0 처리함. 산정 대상은 저장량이며 지역 GIR 연간 흡수량과 단위·시간 범위를 맞추어 해석해야 함.
경계 공식 여부: {c['boundary']['official']}. 금정산 예제 경계는 4주차 자체 경계이며 공식 국립공원 경계가 아님.
면적은 {c['spatial']['crs']}에서 계산했음. 중복 우선순위: {c['spatial']['overlap_priority']}.
중복 제거 면적 {audit['overlap_removed_ha']:,.6f} ha. 완전 중복 행은 0 면적으로 보존함. 미분류 포함 면적 {unknown:,.6f} ha.

## 재적 대체와 결측
{markdown_table(fallback)}

상세 규칙별 건수·면적은 provenance.json의 volume.audit에 기록함. 래스터는 유효 교차 면적으로만 평균하며, 부분 피복 및 전체 결측 면적을 따로 기록함. 임의 대체하지 않음.

## 산정 제외 면적
{markdown_table(excluded)}

## 불확도
{markdown_table(pd.DataFrame([{k:v for k,v in mc.items() if k not in {'overrides_applied','override_input','high_uncertainties'}}]))}

적용한 불확도 덮어쓰기(원본 계수 CSV는 변경하지 않음):
{markdown_table(overrides_table)}

적용 후 {c['uncertainty']['warning_threshold_pct']}% 초과 계수(경고만 하고 계속함):
{markdown_table(high_table)}

덮어쓰기 적용·미적용 비교(같은 구획·재적, n·seed 동일):
{markdown_table(sensitivity)}

미적용 비교에 남은 고불확도 계수:
{markdown_table(pd.DataFrame([dict(scenario=r['scenario'],**h) for r in uncertainty_audits for h in r.get('high_uncertainties',[])]))}

계수 unc_pct는 95% 신뢰구간 반폭(%)으로 간주하여 sd=값×unc_pct/196으로 표본화함. 같은 수종은 모든 구획·혼효 성분에서 같은 표본값을 공유함. 양수가 나올 때까지 재추출함. 빈 불확도와 CF는 고정함. forest_type_3와 national_inventory_2025는 계수별 불확도 자료가 없어 모든 계수를 고정함. 고정 세트의 5–95%가 같은 값인 것은 실제 오차가 없다는 뜻이 아님. 재적·면적 불확도는 포함하지 않음. 보고 구간은 5–95%임. 곰솔 112.782→11.62% 변경은 원문 오기 추정에 대한 팀의 덮어쓰기이며 공식 정오표로 단정하지 않음.

## 지정 시험 비교
{markdown_table(comparison)}

{note}

## 국가 인벤토리 요인별 비교
{markdown_table(inventory)}

각 변화량은 2025 기준 대비임. ‘둘 다’에서 ‘계수 세트만’과 ‘혼효림 규칙만’의 변화량을 뺀 값은 상호작용임. 두 효과를 단순 합산하지 않음.

## 외부 참고값 비교
{markdown_table(ext)}

신재원 외(2023) 표 3은 2020년 Landsat 지수+SRTM 랜덤포레스트 재적과 토지피복 산림면적을 사용한 다른 방법의 결과임. 정답이나 검증 목표가 아님. 시도 일부 AOI는 비교하지 않음. 전체 시도 비교는 공식 시도 도형의 피복 검증을 통과해야 함. 원문 표의 지역값과 합계는 재계산으로 교정하지 않고 기재 수치를 보존함.

## 한계와 재현 정보
시도 통계는 개별 임분을 직접 측정한 값이 아님. 임상도 갱신년도, 재적 통계 기준연도, 계수 판이 다를 수 있음. Ⅰ영급의 재적 0은 국가통계 정의이며 생물학적으로 탄소가 없다는 뜻이 아님. 모든 입력 SHA-256·config·패키지 버전·실행 시각은 provenance.json에 기록함.
'''
    (out/'report.md').write_text(text,encoding='utf-8')
    return total
