"""Parameter averaging and explicit volume-split component accounting."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

PARAMETERS=('WD','BEF','R','CF')

def conversion(values):
    return values['WD']*values['BEF']*(1+values['R'])*values['CF']

class FactorSet:
    def __init__(self,config):
        self.config=config;ref=Path(config['reference_dir'])
        self.species=pd.read_csv(ref/'nifos2014_carbon_factors.csv').set_index('species')
        self.cf=pd.read_csv(ref/'nifos2010_carbon_fraction.csv').set_index('species').CF.to_dict()
        self.types=pd.read_csv(ref/'forest_type_factors_shin2023.csv').set_index('forest_type')
        self.inventory=pd.read_csv(ref/'national_inventory_factors_nid2025.csv',dtype={'forest_code':str}).set_index('forest_type')
        mapping=pd.read_csv(ref/'species_group_to_factor.csv',dtype={'forest_code':str,'species_code':str}).fillna('')
        mapping=mapping[mapping.profile==config['factors']['mapping_profile']]
        self.mapping={(str(fc),str(sc)):m for (fc,sc),m in mapping.groupby(['forest_code','species_code'],sort=False)}
        self.cf_rules=pd.read_csv(ref/'cf_rules.csv').set_index('forest_type')
        self.weights={};self.component_specs={}
        self.fixed=config['factors']['set'] in {'forest_type_3','national_inventory_2025'}

    def resolve(self,row):
        m=self.mapping.get((row.FRTP_CD,row.KOFTR_GROU))
        if m is None:m=self.mapping.get((row.FRTP_CD,'*'))
        if m is None:return None,'미분류','No CSV mapping'
        weights=dict(zip(m.target_species,m.weight.astype(float)))
        if not np.isclose(sum(weights.values()),1):raise ValueError('Factor weights do not sum to one')
        if set(weights)-set(self.species.index):raise ValueError('Unknown target species')
        return weights,m.rule.iloc[0],m.note.iloc[0]

    def cf_value(self,forest_type,weights=None,fixed=None):
        f=self.config['factors'];mode=f['cf_mode']
        if mode=='from_factor_set':return float(fixed['CF'])
        if mode=='ipcc_inventory':return float(self.cf_rules.loc[forest_type,'ipcc_inventory'])
        if mode=='uniform_0_5':return f['uniform_cf']
        if mode=='nifos2010':return sum(w*self.cf.get(s,f['missing_cf']) for s,w in weights.items())
        raise ValueError('Unknown CF mode')

    def group_component(self,group):
        c=self.config;ft=c['table_parser']['conifer_type' if group=='conifer' else 'broadleaf_type']
        spec=c['factors']['split_groups'][group]
        target=self.species[self.species.forest_type==spec['species_type']]
        if not len(target):raise ValueError('Empty species factor group')
        ws={s:1/len(target) for s in target.index}
        if c['factors']['set']=='national_inventory_2025':values=self.inventory.loc[spec['inventory_type'],list(PARAMETERS)].to_dict()
        elif c['factors']['set']=='forest_type_3':values=self.types.loc[ft,list(PARAMETERS)].to_dict()
        else:values={k:sum(w*float(self.species.loc[s,k]) for s,w in ws.items()) for k in PARAMETERS if k!='CF'}
        values['CF']=self.cf_value(ft,ws,values)
        return {'weight':c['factors']['mixed_weights'][group],'weights':ws,'values':values,'group':group}

    def attach(self,g):
        c=self.config;rows=[];mixed_type=c['table_parser']['mixed_type']
        for row in g.itertuples():
            rule='직접';note='';ws={};components=[];values={k:np.nan for k in PARAMETERS};key=''
            if not row.eligible:
                values={k:0. for k in PARAMETERS};rule='제외';note='죽림 NE / 무립목지·비산림 0; 미분류 임상 별도 보고'
            elif row.forest_type==mixed_type and (c['factors']['mixed_rule']=='volume_split_50_50' or c['factors']['set']=='national_inventory_2025'):
                parts=[self.group_component(group) for group in c['factors']['mixed_weights']]
                if not np.isclose(sum(p['weight'] for p in parts),1):raise ValueError('Mixed volume weights must sum to one')
                values={k:sum(p['weight']*p['values'][k] for p in parts) for k in PARAMETERS}
                if c['factors']['mixed_rule']=='volume_split_50_50':
                    components=parts;rule='재적분할';note='침엽·활엽에 재적 50%씩 배분 후 탄소 합산; 수종그룹 미사용'
                else:
                    components=[{'weight':1.,'weights':{},'values':values}];rule='그룹평균';note='침엽·활엽 계수 각각 1:1 평균 후 곱함'
            elif c['factors']['set']=='national_inventory_2025':
                match=self.inventory[self.inventory.forest_code.eq(row.FRTP_CD)]
                if len(match)==1:
                    values={k:float(match.iloc[0][k]) for k in PARAMETERS};note='NID 2025 표 6-9; 활엽수는 상록활엽 포함; 수종그룹 미사용'
                    components=[{'weight':1.,'weights':{},'values':values}]
                else:rule='미분류';note='Inventory forest type missing'
            elif c['factors']['set']=='forest_type_3':
                if row.forest_type in self.types.index:
                    ws,_,_=self.resolve(row)
                    values={k:float(self.types.loc[row.forest_type,k]) for k in PARAMETERS}
                    values['CF']=self.cf_value(row.forest_type,ws,values)
                    components=[{'weight':1.,'weights':{},'values':values}];note='신재원 외 2023 표 2'
                else:rule='미분류';note='Forest type factor missing'
            else:
                ws,rule,note=self.resolve(row)
                if ws is not None:
                    values.update({k:sum(w*float(self.species.loc[s,k]) for s,w in ws.items()) for k in PARAMETERS if k!='CF'})
                    values['CF']=self.cf_value(row.forest_type,ws,values)
                    components=[{'weight':1.,'weights':ws,'values':values}]
            if components:
                key=json.dumps([{'weight':p['weight'],'weights':p['weights'],'CF':p['values']['CF'],'fixed_values':p['values'] if self.fixed else None} for p in components],ensure_ascii=False,sort_keys=True)
                self.component_specs[key]=components;self.weights[key]=ws or {}
                coefficient=sum(p['weight']*conversion(p['values']) for p in components)
            else:coefficient=0. if not row.eligible else np.nan
            output=dict(**values,factor_key=key,factor_rule=rule,factor_note=note,factor_targets=';'.join(ws or {}),
                        mixed_rule=c['factors']['mixed_rule'],carbon_conversion_tC_m3=coefficient)
            if row.eligible and row.forest_type==mixed_type and c['factors']['mixed_rule']=='volume_split_50_50':
                for p in components:
                    for k in PARAMETERS:output[p['group']+'_'+k]=p['values'][k]
            rows.append(output)
        stale=set(PARAMETERS)|{'factor_key','factor_rule','factor_note','factor_targets','mixed_rule','carbon_conversion_tC_m3'}
        stale|={group+'_'+k for group in c['factors']['mixed_weights'] for k in PARAMETERS}
        return pd.concat([g.drop(columns=list(stale),errors='ignore').reset_index(drop=True),pd.DataFrame(rows)],axis=1)
