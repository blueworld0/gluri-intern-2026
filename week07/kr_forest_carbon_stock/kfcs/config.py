"""Explicit configuration; paths are relative to the experiment, never cwd."""
import copy
import os
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]

def load_config(path, drive_root=None, output_dir=None):
    path = Path(path).resolve()
    c = yaml.safe_load(path.read_text(encoding='utf-8'))
    c = copy.deepcopy(c)
    c['factors'].setdefault('mixed_rule','volume_split_50_50' if c['factors']['set']=='national_inventory_2025' else 'factor_average_50_50')
    if c['factors']['mixed_rule']=='conifer_broadleaf_50_50': c['factors']['mixed_rule']='factor_average_50_50'
    c['uncertainty'].setdefault('overrides','reference/uncertainty_overrides.csv')
    c['drive_root'] = str(Path(drive_root or os.environ.get('KFCS_DRIVE_ROOT') or c['drive_root']).expanduser().resolve())
    def expand(value):
        if isinstance(value, dict): return {k: expand(v) for k,v in value.items()}
        if isinstance(value, list): return [expand(v) for v in value]
        if isinstance(value, str):
            return value.replace('${drive_root}', c['drive_root']).replace('${experiment_root}', str(ROOT))
        return value
    c = expand(c)
    if output_dir: c['output_dir'] = str(Path(output_dir).resolve())
    for key, allowed in [
        (('volume','provider'), {'stat_table','raster'}),
        (('volume','age1_rule'), {'zero','use_age2'}),
        (('factors','set'), {'nifos2014_species','forest_type_3','national_inventory_2025'}),
        (('factors','cf_mode'), {'ipcc_inventory','uniform_0_5','nifos2010','from_factor_set'}),
        (('factors','mixed_rule'), {'factor_average_50_50','volume_split_50_50'}),
        (('factors','mapping_profile'), {'legacy_week4','national_groups'}),
    ]:
        value=c[key[0]][key[1]]
        if value not in allowed: raise ValueError(f'{key}: {value}; allowed={allowed}')
    if c['volume']['provider']=='stat_table' and c['volume']['table'] not in {'sfs2025','kosis2020'}:
        raise ValueError('Unknown volume table')
    if c['volume']['provider']=='raster':
        if c['volume']['raster_unit'] not in {'m3_per_ha','m3_per_pixel'}: raise ValueError('Unknown raster unit')
        if not c['volume'].get('raster_path'): raise ValueError('raster_path required')
    if c['factors']['cf_mode']=='from_factor_set' and c['factors']['set'] not in {'forest_type_3','national_inventory_2025'}:
        raise ValueError('from_factor_set requires a forest-type factor set')
    if c['factors']['set']=='national_inventory_2025' and c['factors']['cf_mode']!='from_factor_set':
        raise ValueError('national_inventory_2025 requires cf_mode=from_factor_set')
    if not c['inputs']['forest_maps']: raise ValueError('forest_maps must not be empty')
    if c['uncertainty']['n'] <= 0: raise ValueError('n must be positive')
    for key in ('reference_dir','output_dir'):
        p=Path(c[key]); c[key]=str(p if p.is_absolute() else ROOT/p)
    overrides=c['uncertainty']['overrides']
    if overrides in (None,'none'): c['uncertainty']['overrides']=None
    else:
        p=Path(overrides);c['uncertainty']['overrides']=str(p if p.is_absolute() else ROOT/p)
    c['_config_path']=str(path)
    return c
