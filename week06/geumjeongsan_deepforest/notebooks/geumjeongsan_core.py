"""Pure Python helpers shared by Colab and standard-library-only local checks."""
import csv
import hashlib
import math
import re
import struct
from pathlib import Path

SHEETS = ('35909073','35909082','35909083','35909091','35909092','35909093','35913001','35913003','35913013')
CONTROL_5179 = (1140963,1693231,1143384,1696139)

def sheet_bounds(sheet):
    if not re.fullmatch(r'\d{8}', str(sheet)):
        raise ValueError('8자리 1:5000 도엽 번호 필요')
    sheet = str(sheet); latitude=int(sheet[:2]); longitude=120+int(sheet[2])
    big=int(sheet[3:5])-1; small=int(sheet[5:])-1
    if not (0<=big<16 and 0<=small<100):raise ValueError('잘못된 도엽 번호')
    row,column=divmod(big,4); sr,sc=divmod(small,10)
    west=longitude+column*.25+sc*.025
    north=latitude+1-row*.25-sr*.025
    return west,north-.025,west+.025,north

def quadkey_xy(x,y,z=9):
    return ''.join(str((1 if x&(1<<i) else 0)+(2 if y&(1<<i) else 0)) for i in range(z-1,-1,-1))

def quadkey_bounds(quad):
    x=y=0;z=len(quad)
    for c in quad:
        digit=int(c)
        if digit not in (0,1,2,3):raise ValueError('invalid quadkey')
        x=x*2+(digit&1);y=y*2+(digit>>1)
    def latitude(v):return math.degrees(math.atan(math.sinh(math.pi*(1-2*v/(2**z)))))
    return x/(2**z)*360-180,latitude(y+1),(x+1)/(2**z)*360-180,latitude(y)

def quadkeys_for_bounds(bounds,z=9):
    west,south,east,north=bounds
    def xy(lon,lat):
        lat=max(-85.05112878,min(85.05112878,lat)); n=2**z
        return max(0,min(n-1,int((lon+180)/360*n))),max(0,min(n-1,int((1-math.asinh(math.tan(math.radians(lat)))/math.pi)/2*n)))
    x0,y0=xy(west,north); x1,y1=xy(east,south)
    return [quadkey_xy(x,y,z) for y in range(y0,y1+1) for x in range(x0,x1+1)]

def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def read_csv_rows(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def tiff_header(path):
    """Read tiny classic-TIFF IFD metadata only; never read pixel arrays."""
    with Path(path).open('rb') as f:
        order={'II':'<','MM':'>'}[f.read(2).decode()]
        if struct.unpack(order+'H',f.read(2))[0]!=42:raise ValueError('Classic TIFF only')
        f.seek(struct.unpack(order+'I',f.read(4))[0]); count=struct.unpack(order+'H',f.read(2))[0]
        tags={}
        for _ in range(count):
            tag,kind,n,raw=struct.unpack(order+'HHI4s',f.read(12))
            if n==1 and kind in (3,4):tags[tag]=struct.unpack(order+('H' if kind==3 else 'I'),raw[:2 if kind==3 else 4])[0]
            else:tags[tag]=dict(type=kind,count=n)
    return dict(width=tags[256],height=tags[257],compression=tags.get(259),samples=tags.get(277),
                has_geotags=any(k in tags for k in (33550,33922,34264,34735)),tags=tags)

def parse_jucker_r(text):
    match=re.search(r'pa\s*<-\s*c\((.*?)\)\s*\n',text,re.S)
    if not match:raise ValueError('Jucker R coefficient vector not found')
    values=[float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',match.group(1))]
    if len(values)!=100:raise ValueError(f'Expected25×4 coefficients; got{len(values)}')
    table={int(values[i]):dict(a=values[i+1],b=values[i+2],g=values[i+3]) for i in range(0,100,4)}
    if set(table)!=set(range(25)):raise ValueError('Unexpected biome numbering')
    return table

def jucker_diameter(height,crown,forest_code,table):
    def one(k):
        c=table[k];return c['a']*(height*crown)**c['b']*math.exp(c['g']**2/2)
    if height<=0 or crown<=0:return float('nan'),'invalid'
    if str(forest_code)=='1':return one(21),'21'
    if str(forest_code)=='2':return one(20),'20'
    return (one(20)+one(21))/2,'mean(20,21)'

class KoreanAllometry:
    def __init__(self,equations,factors,fractions):
        self.eq={(r['species'],r['form'],r['part']):r for r in equations}
        self.factors={r['species']:r for r in factors}
        self.cf={r['species']:float(r['CF']) for r in fractions}
        self.dec=[s for s,r in self.factors.items() if r['forest_type']=='낙엽활엽수']
        self.con=[s for s,r in self.factors.items() if r['forest_type']=='침엽수']
        self.pines=['강원지방소나무','중부지방소나무']
        self.oaks=['상수리나무','굴참나무','신갈나무','졸참나무']
        assert len(self.factors)==22 and len(self.eq)==176
    def mapping(self,species,dbh,forest_code='',origin='field'):
        aliases={'낙엽송':'일본잎갈나무','편백나무':'편백','개서어나무':'서어나무','소사나무':'서어나무','백합나무':'튤립나무'}
        if species=='등나무':return [],'등나무 제외'
        if species=='소나무':return [('강원지방소나무' if dbh>40 else '중부지방소나무',1.)],'소나무40cm 규칙'
        if species in ('떡갈나무','기타참나무류'):targets=self.oaks;rule='참나무4종 예측량 평균'
        elif aliases.get(species,species) in self.factors:targets=[aliases.get(species,species)];rule='직접/README 별칭'
        elif origin=='field':
            targets=self.con if species in ('노간주나무','기타침엽수') else self.dec
            rule='기타 침엽수 평균' if targets is self.con else '기타 낙엽활엽수 평균'
        elif str(forest_code)=='1':targets=self.pines;rule='수종 미대응→소나무류 평균'
        elif str(forest_code)=='2':targets=self.dec;rule='수종 미대응→낙엽활엽수 평균'
        else:
            return [(s,.5/len(self.pines)) for s in self.pines]+[(s,.5/len(self.dec)) for s in self.dec], '혼효/임상미상→소나무류·낙엽활엽수 그룹 평균'
        return [(s,1/len(targets)) for s in targets],rule
    def predict(self,species,dbh,forest_code='',origin='field'):
        ws,rule=self.mapping(species,dbh,forest_code,origin)
        if not ws:return dict(AGB_kg=0.,carbon_kg=0.,carbon_cf_kg=0.,formula='excluded',mapping=rule,range_flag='',CF=.5,CF_fallback_weight=0.)
        totals={}; flags=set(); descriptions=[]; species_agb={s:0. for s,w in ws}
        for part in ('줄기','가지','잎'):
            totals[part]=0.
            for s,w in ws:
                r=self.eq[s,'a*D^b',part]
                predicted=float(r['a'])*dbh**float(r['b'])
                totals[part]+=w*predicted;species_agb[s]+=predicted
                descriptions.append(f'{s}:{part}:{r["a"]}*D^{r["b"]}:w={w:.8g}:p{r["page"]}')
                if not float(r['DBH_min_cm'])<=dbh<=float(r['DBH_max_cm']):flags.add(s)
        agb=sum(totals.values());carbon_cf=sum(w*species_agb[s]*self.cf.get(s,.5) for s,w in ws)
        cf=carbon_cf/agb if agb>0 else .5
        return dict(AGB_kg=agb,stem_kg=totals['줄기'],branch_kg=totals['가지'],leaf_kg=totals['잎'],
                    carbon_kg=agb*.5,carbon_cf_kg=carbon_cf,CF=cf,
                    CF_fallback_weight=sum(w for s,w in ws if s not in self.cf),mapping=rule,
                    range_flag=';'.join(sorted(flags)),formula=' | '.join(descriptions),
                    target_species=';'.join(f'{s}:{w:.8g}' for s,w in ws))

def training_ratio(reference,estimated):
    good=[(float(a),float(b)) for a,b in zip(reference,estimated) if math.isfinite(float(a)) and math.isfinite(float(b))]
    denom=sum(b for a,b in good)
    return sum(a for a,b in good)/denom if denom>0 else float('nan')

def dimensional_checks():
    assert 1000/1000*10000/400*.5==12.5
    assert math.isclose(training_ratio([10,30],[5,15]),2)
    assert math.isnan(training_ratio([10],[0]))
    assert all(math.isclose(a,b) for a,b in zip(sheet_bounds('35913003'),(129.05,35.225,129.075,35.25)))
