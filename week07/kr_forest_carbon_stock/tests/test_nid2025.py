"""NID overrides, volume allocation and fixed-factor uncertainty contracts."""
import copy
import hashlib
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
import pytest
from kfcs.config import ROOT,load_config
from kfcs.factors import FactorSet
from kfcs.carbon import calculate
from kfcs.uncertainty import effective_uncertainties,monte_carlo

def config():return load_config(ROOT/'configs/geumjeongsan_2025.yaml')

def mixed():
    return pd.DataFrame({'eligible':[True],'FRTP_CD':['3'],'forest_type':['혼효림'],'KOFTR_GROU':['77'],
                         'area_ha':[1.],'m3_per_ha':[100.]})

def test_override_preserves_raw_file_and_in_memory_source():
    c=config();f=FactorSet(c);path=Path(c['reference_dir'])/'nifos2014_carbon_factors.csv'
    before=hashlib.sha256(path.read_bytes()).hexdigest()
    table,applied,high,source=effective_uncertainties(f,c)
    assert table.loc['곰솔','BEF_unc_pct']==11.62
    assert f.species.loc['곰솔','BEF_unc_pct']==112.782
    assert len(applied)==1 and high==[] and source['sha256']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==before

def test_disabled_override_warns_and_keeps_high_value():
    c=config();c['uncertainty']['overrides']=None
    with pytest.warns(RuntimeWarning,match='112.782'):
        table,applied,high,_=effective_uncertainties(FactorSet(c),c)
    assert table.loc['곰솔','BEF_unc_pct']==112.782
    assert not applied and high[0]['species']=='곰솔'

def test_stale_override_fails_closed(tmp_path):
    c=config();p=tmp_path/'override.csv';df=pd.read_csv(c['uncertainty']['overrides'])
    df['original_unc_pct']=111.;df.to_csv(p,index=False)
    c['uncertainty']['overrides']=str(p)
    with pytest.raises(ValueError,match='original value mismatch'):effective_uncertainties(FactorSet(c),c)

def test_inventory_two_mixed_rules_have_distinct_expected_carbon():
    c=config();c['factors'].update(set='national_inventory_2025',cf_mode='from_factor_set',mixed_rule='factor_average_50_50')
    avg=calculate(FactorSet(c).attach(mixed())).carbon_tC.iloc[0]
    c['factors']['mixed_rule']='volume_split_50_50'
    split=calculate(FactorSet(c).attach(mixed())).carbon_tC.iloc[0]
    assert avg==pytest.approx(100*.57*1.47*(1+.315)*.495)
    assert split==pytest.approx(50*.46*1.43*(1+.27)*.51+50*.68*1.51*(1+.36)*.48)
    assert split!=pytest.approx(avg)

def test_inventory_ignores_species_codes_and_uncertainty_is_fixed():
    c=load_config(ROOT/'configs/geumjeongsan_national_inventory.yaml')
    g=mixed().drop(columns='KOFTR_GROU');f=FactorSet(c);g=calculate(f.attach(g));m=monte_carlo(g,f,c)
    assert m['coefficient_fixed'] and m['overrides_applied']==[]
    assert m['p05_tC']==m['median_tC']==m['p95_tC']==pytest.approx(g.carbon_tC.sum())

def test_species_volume_split_preserves_shared_samples():
    c=config();c['factors']['mixed_rule']='volume_split_50_50'
    f=FactorSet(c);whole=monte_carlo(calculate(f.attach(mixed())),f,c)
    halves=pd.concat([mixed(),mixed()],ignore_index=True);halves['area_ha']/=2
    f2=FactorSet(c);split=monte_carlo(calculate(f2.attach(halves)),f2,c)
    assert whole==split

def test_inventory_defaults_and_cf_guard(tmp_path):
    import yaml
    c=config();c['factors']['set']='national_inventory_2025';c['factors']['cf_mode']='from_factor_set'
    c['factors'].pop('mixed_rule');c['uncertainty'].pop('overrides')
    p=tmp_path/'config.yaml';p.write_text(yaml.safe_dump(c,allow_unicode=True),encoding='utf-8')
    actual=load_config(p)
    assert actual['factors']['mixed_rule']=='volume_split_50_50'
    assert Path(actual['uncertainty']['overrides']).name=='uncertainty_overrides.csv'
    c['factors']['cf_mode']='uniform_0_5';p.write_text(yaml.safe_dump(c,allow_unicode=True),encoding='utf-8')
    with pytest.raises(ValueError,match='requires cf_mode'):load_config(p)

def test_sensitivity_also_runs_applied_case_when_current_config_is_none(tmp_path):
    from kfcs.report import uncertainty_comparison
    c=config();c['uncertainty']['overrides']=None
    f=FactorSet(c);g=calculate(f.attach(mixed()))
    with pytest.warns(RuntimeWarning,match='112.782'):mc=monte_carlo(g,f,c)
    table,audits=uncertainty_comparison(g,c,f,mc,tmp_path)
    assert table.scenario.tolist()==['미적용','적용']
    assert table.overrides.tolist()==['none','uncertainty_overrides.csv']
    assert audits[0]['high_uncertainties'] and len(audits[1]['overrides_applied'])==1
