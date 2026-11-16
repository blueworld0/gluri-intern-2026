"""Auditable uncertainty overrides and species-correlated component sampling."""
import warnings
import numpy as np
import pandas as pd
from .volume.base import file_record

def effective_uncertainties(factors,config):
    table=factors.species.copy(deep=True);applied=[];high=[];source=None
    if factors.fixed:return table,applied,high,source
    u=config['uncertainty'];path=u.get('overrides')
    if path not in (None,'none'):
        source=file_record(path);overrides=pd.read_csv(path)
        keys=['factor_set','species','parameter']
        if overrides.duplicated(keys).any():raise ValueError('Duplicate uncertainty overrides')
        for row in overrides[overrides.factor_set.eq(config['factors']['set'])].to_dict('records'):
            s=row['species'];k=row['parameter']
            if s not in table.index or k not in u['parameters']:raise ValueError('Unknown uncertainty override key')
            original=table.loc[s,k+'_unc_pct'];replacement=float(row['override_unc_pct'])
            if not np.isfinite(original) or not np.isclose(original,float(row['original_unc_pct'])):
                raise ValueError('Uncertainty override original value mismatch')
            if not np.isfinite(replacement) or replacement<0:raise ValueError('Invalid uncertainty override value')
            table.loc[s,k+'_unc_pct']=replacement;applied.append(row)
    for s,row in table.sort_index().iterrows():
        for k in u['parameters']:
            value=row[k+'_unc_pct']
            if pd.notna(value) and value>u['warning_threshold_pct']:
                high.append({'factor_set':config['factors']['set'],'species':s,'parameter':k,'effective_unc_pct':float(value)})
    if high:
        warnings.warn(f"unc_pct > {u['warning_threshold_pct']}% 계수 {len(high)}개 남아 있음: "+'; '.join(f"{r['species']} {r['parameter']}={r['effective_unc_pct']}%" for r in high),RuntimeWarning,stacklevel=2)
    return table,applied,high,source

def monte_carlo(g,factors,config):
    u=config['uncertainty']
    if not u['enabled']:return {'enabled':False,'overrides_applied':[],'high_uncertainties':[]}
    table,applied,high,source=effective_uncertainties(factors,config)
    rng=np.random.default_rng(u['seed']);n=u['n'];totals=np.zeros(n)
    if factors.fixed:totals[:]=g.carbon_tC.sum(min_count=1)
    else:
        sampled={}
        for s,row in table.sort_index().iterrows():
            sampled[s]={}
            for k in u['parameters']:
                value=float(row[k]);unc=row[k+'_unc_pct'];a=np.full(n,value)
                if np.isfinite(unc) and unc>0:
                    sd=value*unc/u['ci_halfwidth_divisor'];a=rng.normal(value,sd,n)
                    while (a<=0).any():
                        mask=a<=0;a[mask]=rng.normal(value,sd,int(mask.sum()))
                sampled[s][k]=a
        valid=g[g.eligible & g.carbon_tC.notna()]
        for key,rows in valid.groupby('factor_key',sort=True):
            exposure=(rows.area_ha*rows.m3_per_ha).sum();coeff=np.zeros(n)
            for component in factors.component_specs[key]:
                p={k:sum(w*sampled[s][k] for s,w in component['weights'].items()) for k in u['parameters']}
                coeff+=component['weight']*p['WD']*p['BEF']*(1+p['R'])*component['values']['CF']
            totals+=exposure*coeff
    q=np.quantile(totals,u['quantiles'])
    return {'enabled':True,'n':n,'seed':u['seed'],'median_tC':float(q[1]),'p05_tC':float(q[0]),'p95_tC':float(q[2]),
            'coefficient_fixed':factors.fixed,'fixed_reason':'계수별 불확도 자료 없음' if factors.fixed else None,
            'overrides_applied':applied,'override_input':source,'high_uncertainties':high,
            'conditional_on_observed_volume':True,'missing_area_ha':float(g.loc[g.carbon_tC.isna(),'area_ha'].sum())}
