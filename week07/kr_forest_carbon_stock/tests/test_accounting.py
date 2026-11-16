"""Test failures that can silently bias totals or misstate uncertainty."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
import pytest
from kfcs.config import load_config,ROOT
from kfcs.factors import FactorSet
from kfcs.carbon import calculate
from kfcs.uncertainty import monte_carlo
from kfcs.volume.stat_table import StatTableProvider
from kfcs.volume.tables import parse_sfs_sheet

def config(): return load_config(ROOT/'configs/geumjeongsan_week4_repro.yaml')

def rows():
    return pd.DataFrame({'sido_code':['26']*3,'forest_type':['침엽수림','혼효림','활엽수림'],'FRTP_CD':['1','3','2'],
        'KOFTR_GROU':['11','77','60'],'eligible':[True]*3,'age_class':[1,7,3],'area_ha':[1.,2.,3.],
        'species_group':['소나무','침활혼효림','기타상록활엽수'],'m3_per_ha':[100.,200.,300.]})

def test_missing_key_fallback_order_and_age_rule():
    c=config();p=StatTableProvider.__new__(StatTableProvider);p.config=c;p.v=c['volume']
    p.lookup={('26','전체',2):10.,('00','침엽수림',2):99.,('00','혼효림',6):30.}
    r=p.provide(rows())
    assert r.m3_per_ha.iloc[:2].tolist()==[10.,30.]
    assert np.isnan(r.m3_per_ha.iloc[2])
    assert r.volume_flag.tolist()==['대체','대체','결측']
    c['volume']['age1_rule']='zero'
    assert p.provide(rows()).m3_per_ha.iloc[0]==0

def test_parameter_means_are_multiplied_after_averaging():
    c=config(); f=FactorSet(c);g=calculate(f.attach(rows()))
    s=f.species
    con=s[s.forest_type=='침엽수'];dec=s[s.forest_type=='낙엽활엽수']
    expected={k:(con[k].mean()+dec[k].mean())/2 for k in ['WD','BEF','R']}
    for k,value in expected.items(): assert g.loc[1,k]==pytest.approx(value)
    assert g.loc[1,'carbon_tC']==pytest.approx(2*200*expected['WD']*expected['BEF']*(1+expected['R'])*.5)
    # Week4 legacy uses the deciduous mean for unspecified evergreen groups.
    assert g.loc[2,'WD']==pytest.approx(dec.WD.mean())

def test_species_samples_shared_across_split_polygons():
    c=config();c['uncertainty']['enabled']=True
    x=rows().iloc[:1].copy();f=FactorSet(c);whole=calculate(f.attach(x))
    a=monte_carlo(whole,f,c)
    halves=pd.concat([x,x],ignore_index=True);halves['area_ha']/=2
    f2=FactorSet(c);b=monte_carlo(calculate(f2.attach(halves)),f2,c)
    assert a==b

def test_forest_type_factors_have_fixed_monte_carlo():
    c=config();c['uncertainty']['enabled']=True;c['factors'].update(set='forest_type_3',cf_mode='from_factor_set')
    f=FactorSet(c);g=calculate(f.attach(rows()));m=monte_carlo(g,f,c)
    assert m['median_tC']==m['p05_tC']==m['p95_tC']==pytest.approx(g.carbon_tC.sum())

def test_sfs_labels_not_absolute_cell_offsets():
    c=config(); p=c['table_parser']; d=pd.DataFrame('',index=range(14),columns=range(22),dtype=object)
    # Shift both the row and column layout away from the real workbook.
    for j,roman in enumerate(p['roman_ages']):d.iat[8,5+j]=roman;d.iat[8,14+j]=roman
    d.iat[8,4]='합계';d.iat[8,13]='합계'
    d.iat[9,1]='총계';d.iat[9,2]='총계';d.iat[9,4]=21;d.iat[9,13]=210
    for j in range(6):d.iat[9,5+j]=j+1;d.iat[9,14+j]=(j+1)*10
    a,v=parse_sfs_sheet(d,c)
    assert a.sum()==21 and v.sum()==210
    d.iat[9,4]=22
    with pytest.raises(AssertionError):parse_sfs_sheet(d,c)

def test_existing_output_is_untouched(tmp_path):
    from kfcs.run import run
    c=config(); c['output_dir']=str(tmp_path)
    sentinel=tmp_path/'original.txt';sentinel.write_text('keep',encoding='utf-8')
    with pytest.raises(FileExistsError):run(c)
    assert list(tmp_path.iterdir())==[sentinel]
    assert sentinel.read_text(encoding='utf-8')=='keep'

def test_zipped_boundary_selects_by_attribute(tmp_path):
    import zipfile
    import geopandas as gpd
    from shapely.geometry import box
    from kfcs.aoi import read_boundary
    folder=tmp_path/'source';folder.mkdir()
    g=gpd.GeoDataFrame({'NPK_NM':['공원A','공원B']},geometry=[box(900000,1900000,900010,1900010),box(900020,1900000,900040,1900010)],crs='EPSG:5179')
    g.to_file(folder/'parks.shp',encoding='cp949')
    path=tmp_path/'parks.zip'
    with zipfile.ZipFile(path,'w') as z:
        for p in folder.iterdir():z.write(p,'공원경계/'+p.name)
    c=config();c['boundary'].update(path=str(path),member='parks.shp',select={'NPK_NM':'공원B'})
    assert read_boundary(c).area==pytest.approx(200)
