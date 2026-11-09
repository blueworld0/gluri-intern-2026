from forest_map_common import *
from analyze_forest_maps import components
def read_merged(name):
    reader=shapefile.Reader(str(MERGED/name),encoding='cp949')
    return [{'geometry':shape(s.shape.__geo_interface__),'attrs':s.record.as_dict()} for s in reader.iterShapeRecords()]
if __name__=='__main__':
    fs=read_merged('forest2025_geumjeong.shp'); ts=read_merged('trails_geumjeong.shp')
    fg=np.array([f['geometry'] for f in fs if stocked(f['attrs'])],dtype=object)
    comp=components(fg,0); tree=shapely.STRtree(fg)
    votes=defaultdict(lambda:defaultdict(float))
    for t in ts:
        name=t['attrs'].get('MNTN_NM','')
        for i in tree.query(t['geometry'],predicate='intersects'):
            votes[int(comp[i])][name]+=fg[i].intersection(t['geometry']).length
    print('components by descending area:')
    for c in sorted(set(comp),key=lambda c:sum(shapely.area(fg[comp==c])),reverse=True)[:25]:
        print(c,round(sum(shapely.area(fg[comp==c]))/10000,2),dict(votes[c]))
    for year in (2013,2025):
        f=read_merged(f'forest{year}_geumjeong.shp')
        print(year,'forest attribute frequencies',__import__('collections').Counter((r['attrs']['FRTP_CD'],r['attrs']['KOFTR_NM'],r['attrs']['AGCLS_CD']) for r in f).most_common(50))
    a=read_csv(ROOT/'data/reference/nifos2014_allometry.csv')
    print('allometry counts',__import__('collections').Counter((r['species'],r['form']) for r in a))
