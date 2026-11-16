"""Label-driven, fail-closed source workbook readers. No fixed cell positions."""
import io
import json
import re
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from .base import file_record

COLUMNS=['sido_name','sido_code','forest_type','age_class','area_ha','volume_m3','m3_per_ha','source','edition']

def norm(x): return re.sub(r'\s+','',str(x)) if not pd.isna(x) else ''

def zip_members(path):
    with zipfile.ZipFile(path) as z:
        for item in z.infolist():
            # ZipFile decodes non-UTF8 names as CP437. Restore original CP949 bytes.
            name=item.filename if item.flag_bits & 0x800 else item.filename.encode('cp437').decode('cp949')
            if not item.is_dir(): yield name,z.read(item)

def numeric(value):
    s=norm(value).replace(',','')
    if s in {'-','―','–'}: return 0.
    return float(s)

def sido_lookup(config):
    df=pd.read_csv(Path(config['reference_dir'])/'sido_codes.csv',dtype=str).fillna('')
    lookup={}
    for row in df.itertuples():
        for alias in (row.sido_name,row.table_name,*row.aliases.split('|')):
            if alias: lookup[norm(alias)]=(row.sido_name,row.canonical_code)
    return lookup

def parse_sfs_sheet(df, config):
    p=config['table_parser']; romans=p['roman_ages']
    text=df.map(norm)
    headers=[i for i,row in text.iterrows() if all((row==r).sum()==2 for r in romans)]
    candidates=[]
    for n,h in enumerate(headers):
        first_age=min(j for j,x in enumerate(text.iloc[h]) if x==romans[0])
        end=headers[n+1] if n+1<len(headers) else len(df)
        totals=[i for i in range(h+1,end) if text.iloc[i,:first_age].eq(p['total_label']).sum()>=p['total_label_count']]
        candidates.extend((h,i) for i in totals)
    if len(candidates)!=1: raise ValueError(f'Expected unique header + overall ownership total, got {candidates}')
    h,r=candidates[0]; age_cols=[[j for j,x in enumerate(text.iloc[h]) if x==r] for r in romans]
    first_age=min(x[0] for x in age_cols)
    sumcols=[j for j,x in enumerate(text.iloc[h]) if x==p['sum_label']]
    area_total=[j for j in sumcols if j<first_age]
    if len(area_total)!=1: raise ValueError('Stocked area total header missing')
    areas=[numeric(df.iat[r,cs[0]]) for cs in age_cols]
    volumes=[numeric(df.iat[r,cs[1]]) for cs in age_cols]
    if not np.isclose(sum(areas),numeric(df.iat[r,area_total[0]]),atol=p['area_assert_tolerance_ha'],rtol=0):
        raise AssertionError('Age area sum differs from stocked forest area')
    volume_total=[j for j in sumcols if first_age<j<age_cols[0][1]]
    if len(volume_total)!=1 or not np.isclose(sum(volumes),numeric(df.iat[r,volume_total[0]]),atol=p['volume_assert_tolerance_m3'],rtol=0):
        raise AssertionError('Age volume sum differs from total')
    return np.array(areas),np.array(volumes)

def parse_sfs2025(config):
    path=config['inputs']['sfs_zip']; lookup=sido_lookup(config); rows=[]; files=0
    p=config['table_parser']
    for name,blob in zip_members(path):
        if p['sfs_member_label'] not in name: continue
        match=re.search(r'\(\d+\.([^)]*)\)',name)
        if not match or norm(match[1]) not in lookup: raise ValueError(f'Unknown sido workbook {name}')
        sido,sc=lookup[norm(match[1])]; book=pd.ExcelFile(io.BytesIO(blob)); sets={}
        for sheet in book.sheet_names:
            df=pd.read_excel(book,sheet_name=sheet,header=None)
            # Source sheet names are 시도-1/2/3, not 총계/침엽수/혼효림.
            content=' '.join(norm(x) for x in df.to_numpy().flat if not pd.isna(x))
            for title,ft in p['sfs_sheet_titles'].items():
                if f'({title})' in content:
                    if ft in sets: raise ValueError(f'Duplicate forest type {name}:{ft}')
                    t=df.map(norm)
                    title_rows=[i for i,row in t.iterrows() if row.str.match(p['sfs_title_pattern']).any()]
                    starts=[i for i in title_rows if t.iloc[i].str.contains(f'({title})',regex=False).any()]
                    if len(starts)!=1: raise ValueError('Source title missing or ambiguous')
                    start=starts[0]; end=next((i for i in title_rows if i>start),len(df))
                    try: sets[ft]=parse_sfs_sheet(df.iloc[start:end].reset_index(drop=True),config)
                    except (ValueError,AssertionError) as exc: raise ValueError(f'{name}:{sheet}: {exc}') from exc
        if set(sets)!=set(p['sfs_sheet_titles'].values()): raise ValueError(f'Missing type sheets {name}')
        total,con,mix=(sets[k] for k in [p['all_type'],p['conifer_type'],p['mixed_type']])
        sets[p['broadleaf_type']]=(total[0]-con[0]-mix[0],total[1]-con[1]-mix[1])
        for ft,(areas,volumes) in sets.items():
            if (areas<0).any() or (volumes<0).any(): raise ValueError('Negative broadleaf residual')
            for age,(a,v) in enumerate(zip(areas,volumes),start=config['classification']['age_min']):
                if a==0 and v!=0: raise ValueError('Positive volume with zero area')
                rows.append([sido,sc,ft,age,a,v,v/a if a>0 else np.nan,f'{name}',config['editions']['sfs2025']])
        files+=1
    if files!=p['expected_region_count']: raise ValueError(f'Expected national+17 region books, got {files}')
    return pd.DataFrame(rows,columns=COLUMNS)

def parse_kosis2020(config):
    path=config['inputs']['kosis_xlsx']; p=config['table_parser']; lookup=sido_lookup(config)
    df=pd.read_excel(path,header=None); text=df.map(norm)
    headers=[i for i,row in text.iterrows() if all(label in row.values for label in p['kosis_age_headers'].values())]
    if len(headers)!=1: raise ValueError('KOSIS age headers missing or ambiguous')
    h=headers[0]; agecols={int(age):list(text.iloc[h]).index(label) for age,label in p['kosis_age_headers'].items()}
    labelcols=[j for j,x in enumerate(text.iloc[h]) if p['kosis_group_header'] in x]
    if len(labelcols)!=2: raise ValueError('KOSIS region/type columns missing')
    rows=[]; region=None
    for i in range(h+1,len(df)):
        rn=text.iat[i,labelcols[0]]; ft=text.iat[i,labelcols[1]]
        if rn:
            if rn not in lookup: raise ValueError(f'Unknown KOSIS sido {rn}')
            region=lookup[rn]
        if not ft: continue
        if ft not in p['forest_types'] or region is None: raise ValueError(f'Unknown KOSIS row {i}: {ft}')
        for age,col in agecols.items():
            raw=df.iat[i,col]; value=np.nan if norm(raw) in p['missing_labels'] else numeric(raw)
            rows.append([*region,ft,age,np.nan,np.nan,value,str(path),config['editions']['kosis2020']])
        rows.append([*region,ft,config['classification']['age_min'],np.nan,np.nan,np.nan,str(path),config['editions']['kosis2020']])
    result=pd.DataFrame(rows,columns=COLUMNS)
    if result.sido_code.nunique()!=p['expected_region_count']: raise ValueError('Incomplete KOSIS regions')
    return result

def ensure_table(config, table):
    ref=Path(config['reference_dir']); path=ref/config['volume']['cache_files'][table]
    source=config['inputs']['sfs_zip' if table=='sfs2025' else 'kosis_xlsx']
    record=file_record(source,edition=config['editions'][table]); meta=path.with_suffix('.provenance.json')
    # A source hash mismatch cannot silently reuse a previous extraction.
    if path.exists():
        if not meta.exists(): raise ValueError(f'Unverifiable volume cache: {path}')
        old=json.loads(meta.read_text(encoding='utf-8'))
        if old['sha256']!=record['sha256'] or old['edition']!=record['edition']:
            raise ValueError(f'Volume cache source changed; use a new reference_dir: {path}')
        if old.get('parser_sha256')!=file_record(__file__)['sha256']:
            raise ValueError(f'Parser changed; regenerate cache in a new reference_dir: {path}')
        return path
    df=parse_sfs2025(config) if table=='sfs2025' else parse_kosis2020(config)
    if df.duplicated(['sido_code','forest_type','age_class']).any(): raise ValueError('Duplicate volume keys')
    ref.mkdir(parents=True,exist_ok=True)
    df.to_csv(path,index=False,encoding='utf-8-sig')
    meta.write_text(json.dumps({**record,'parser_sha256':file_record(__file__)['sha256']},ensure_ascii=False,indent=2),encoding='utf-8')
    return path

def administrative_reference(config):
    """Return long reference rows using paired area/volume group headers."""
    p=config['table_parser']; district=config['validation']['admin']['district_name']; rows=[]
    for name,blob in zip_members(config['inputs']['sfs_zip']):
        if p['admin_member_label'] not in name: continue
        book=pd.ExcelFile(io.BytesIO(blob))
        for sheet,dimension in p['admin_sheets'].items():
            df=pd.read_excel(book,sheet_name=sheet,header=None); t=df.map(norm)
            matches=[i for i,row in t.iterrows() if district in row.values]
            if len(matches)!=1: raise ValueError('District not unique in admin reference')
            r=matches[0]
            # Find the major group header and its category row, never absolute addresses.
            hdr=[i for i,row in t.iterrows() if p['admin_area_label'] in row.values and p['admin_volume_label'] in row.values]
            if len(hdr)!=1: raise ValueError('Admin area/volume headers missing')
            h=hdr[0]; ac=list(t.iloc[h]).index(p['admin_area_label']); vc=list(t.iloc[h]).index(p['admin_volume_label'])
            cat=t.iloc[h+1]
            for j in range(ac,vc):
                label=cat.iloc[j]
                if label not in p['admin_categories'][dimension]: continue
                vcols=[k for k in range(vc,len(cat)) if cat.iloc[k]==label]
                rows.append({'dimension':dimension,'group':p['admin_categories'][dimension][label],'reference_area_ha':numeric(df.iat[r,j]),
                             'reference_volume_m3':numeric(df.iat[r,vcols[0]]) if len(vcols)==1 else np.nan,'source':f'{name}:{sheet}'})
    if not rows: raise ValueError('Administrative reference workbook not found')
    return pd.DataFrame(rows)
