"""One calculation shared by every volume supplier."""
import numpy as np

def calculate(g):
    g=g.copy()
    g['volume_m3']=g.area_ha*g.m3_per_ha
    coefficient=g.carbon_conversion_tC_m3 if 'carbon_conversion_tC_m3' in g else g.WD*g.BEF*(1+g.R)*g.CF
    g['carbon_density_tC_ha']=g.m3_per_ha*coefficient
    g['carbon_tC']=g.area_ha*g.carbon_density_tC_ha
    g['CO2_tCO2']=g.carbon_tC*(44/12)
    # Excluded types have zero stock; unidentified eligible records stay missing.
    g.loc[~g.eligible,['m3_per_ha','volume_m3','carbon_density_tC_ha','carbon_tC','CO2_tCO2']]=0.
    g['calculation_flag']=np.where(g.eligible & g.carbon_tC.isna(),'결측','계산됨')
    return g
