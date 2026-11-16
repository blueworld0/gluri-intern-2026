import numpy as np
import pandas as pd
from .base import VolumeProvider, file_record
from .tables import ensure_table

class StatTableProvider(VolumeProvider):
    def __init__(self,config):
        self.config=config; self.v=config['volume']; self.path=ensure_table(config,self.v['table'])
        self.table=pd.read_csv(self.path,dtype={'sido_code':str})
        self.lookup=self.table.set_index(['sido_code','forest_type','age_class']).m3_per_ha.to_dict()
        self.audit={}

    def provide(self,polygons):
        c=self.config; out=[]; reasons=[]
        lo=c['classification']['age_min']; hi=c['classification']['age_max']
        for row in polygons.itertuples():
            source=f'stat_table:{self.v["table"]}:{c["editions"][self.v["table"]]}'
            value=np.nan; flag='결측'; reason='missing_key'
            if not row.eligible:
                value=0.; flag='직접'; reason='excluded_forest_type'
            elif pd.notna(row.age_class) and row.age_class>=lo:
                age=min(hi,int(row.age_class))
                if age==lo and self.v['age1_rule']=='zero': value=0.; flag='직접'; reason='age1_zero'
                else:
                    sub=age==lo
                    if sub: age=c['classification']['age2']
                    candidates=[(str(row.sido_code),row.forest_type,age),(str(row.sido_code),c['table_parser']['all_type'],age),
                                (c['classification']['national_code'],row.forest_type,age)]
                    for rank,key in enumerate(candidates):
                        v=self.lookup.get(key,np.nan)
                        if pd.notna(v):
                            value=float(v); flag='대체' if rank or sub else '직접'
                            reason=('age1_to_age2;' if sub else '')+['direct','sido_all','national_same_type'][rank]
                            source+=':'+':'.join(map(str,key)); break
            out.append((value,source,flag,reason)); reasons.append(reason)
        result=pd.DataFrame(out,index=polygons.index,columns=['m3_per_ha','volume_source','volume_flag','volume_reason'])
        self.audit={}
        for reason in sorted(set(reasons)):
            mask=result.volume_reason.eq(reason)
            self.audit[reason]={'count':int(mask.sum()),'area_ha':float(polygons.loc[mask,'area_ha'].sum())}
        return result

    def provenance(self):
        return {'provider':'stat_table','cache':file_record(self.path,edition=self.config['editions'][self.v['table']]),
                'source':file_record(self.config['inputs']['sfs_zip' if self.v['table']=='sfs2025' else 'kosis_xlsx']),
                'rules':self.v,'audit':self.audit}
