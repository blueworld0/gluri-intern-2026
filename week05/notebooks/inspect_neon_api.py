"""Read-only metadata inspection; Python standard library only."""
import json
import urllib.request
from pathlib import Path

dest=Path(__file__).resolve().parent/'api_evidence'
dest.mkdir(exist_ok=True)
urls={
 'rgb_product':'https://data.neonscience.org/api/v0/products/DP3.30010.001',
 'chm_product':'https://data.neonscience.org/api/v0/products/DP3.30015.001',
 'ornl_location':'https://data.neonscience.org/api/v0/locations/ORNL_001.basePlot.vst',
 'ornl_site':'https://data.neonscience.org/api/v0/sites/ORNL',
 'nu_source':'https://raw.githubusercontent.com/NEONScience/NEON-utilities-python/main/neonutilities/by_tile_aop.py',
}
for label,url in urls.items():
    try:
        raw=urllib.request.urlopen(url,timeout=90).read()
        (dest/(label+('.py' if label=='nu_source' else '.json'))).write_bytes(raw)
        if label.endswith('product'):
            data=json.loads(raw)['data']
            sites=data.get('siteCodes',[])
            for site in sites:
                if site['siteCode'] in ['ORNL','GRSM']: print(label,json.dumps(site))
        elif label=='nu_source': print(label,raw[:180].decode())
        else:
            data=json.loads(raw)['data']
            print(label,json.dumps({k:v for k,v in data.items() if 'utm' in k.lower() or 'zone' in k.lower() or k in ['siteCode','siteLatitude','siteLongitude','siteUtmZone']}))
    except Exception as e:
        print(label,type(e).__name__,str(e))
