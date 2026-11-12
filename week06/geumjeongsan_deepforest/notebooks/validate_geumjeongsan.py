"""Standard library only. TIFF headers/XML/math, then notebook syntax/order."""
import ast,csv,json,math,struct,sys,xml.etree.ElementTree as ET
from pathlib import Path
from geumjeongsan_core import *
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent

def tm(lon,lat,epsg):
    # Independent Snyder series for local control verification only.
    a=6378137.;f=1/298.257222101;e2=f*(2-f);ep=e2/(1-e2)
    lon0,k,x0,y0=(129,1,200000,600000) if epsg==5187 else (127.5,.9996,1000000,2000000)
    p=math.radians(lat);p0=math.radians(38);lam=math.radians(lon-lon0)
    def meridian(p):return a*((1-e2/4-3*e2**2/64-5*e2**3/256)*p-(3*e2/8+3*e2**2/32+45*e2**3/1024)*math.sin(2*p)+(15*e2**2/256+45*e2**3/1024)*math.sin(4*p)-35*e2**3/3072*math.sin(6*p))
    n=a/math.sqrt(1-e2*math.sin(p)**2);t=math.tan(p)**2;c=ep*math.cos(p)**2;A=math.cos(p)*lam
    x=x0+k*n*(A+(1-t+c)*A**3/6+(5-18*t+t*t+72*c-58*ep)*A**5/120)
    y=y0+k*(meridian(p)-meridian(p0)+n*math.tan(p)*(A*A/2+(5-t+9*c+4*c*c)*A**4/24+(61-58*t+t*t+600*c-330*ep)*A**6/720))
    return x,y

def inverse_tm(x,y,epsg):
    lon=129.;lat=35.25
    for _ in range(8):
        px,py=tm(lon,lat,epsg);xx,xy=tm(lon+1e-5,lat,epsg);yx,yy=tm(lon,lat+1e-5,epsg)
        a=(xx-px)/1e-5;b=(yx-px)/1e-5;c=(xy-py)/1e-5;d=(yy-py)/1e-5
        det=a*d-b*c;dx=x-px;dy=y-py
        lon+=(d*dx-b*dy)/det;lat+=(-c*dx+a*dy)/det
    return lon,lat

def main():
    ortho=ROOT/'data/ortho/2021';rows=[]
    for sheet in SHEETS:
        path=ortho/f'(B060)정사영상_2021_{sheet}.tif';head=tiff_header(path)
        xml=ET.parse(next(ortho.glob('*'+sheet+'.xml'))).getroot()
        assert float(xml.findtext('.//지상표본거리'))==.25
        assert xml.findtext('.//투영원점')=='동부' and xml.findtext('.//기준타원체')=='GRS80'
        assert not head['has_geotags'] and head['samples']==3 and head['compression']==1
        w,s,e,n=sheet_bounds(sheet);corners=[tm(x,y,5187) for x in (w,e) for y in (s,n)]
        left=min(x for x,y in corners)-50;top=max(y for x,y in corners)+50
        right=left+head['width']*.25;bottom=top-head['height']*.25
        widths=max(x for x,y in corners)-min(x for x,y in corners)+100
        heights=max(y for x,y in corners)-min(y for x,y in corners)+100
        out=dict(sheet=sheet,width=head['width'],height=head['height'],left=left,bottom=bottom,right=right,top=top,
                 width_difference_m=head['width']*.25-widths,height_difference_m=head['height']*.25-heights)
        assert abs(out['width_difference_m'])<5 and abs(out['height_difference_m'])<5
        if sheet=='35913003':
            transformed=[tm(*inverse_tm(x,y,5187),5179) for x in (left,right) for y in (bottom,top)]
            bounds=(min(x for x,y in transformed),min(y for x,y in transformed),max(x for x,y in transformed),max(y for x,y in transformed))
            errors=[a-b for a,b in zip(bounds,CONTROL_5179)]
            assert max(map(abs,errors))<=3,errors
            out['control_5179_bbox']=bounds;out['control_5179_errors_m']=errors
        rows.append(out)
    with (HERE/'georef_header_check_local.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));writer.writeheader();writer.writerows(rows)
    plots=read_csv_rows(ROOT/'data/field/seo2011_plots.csv');trees=read_csv_rows(ROOT/'data/field/seo2011_trees.csv');layers=read_csv_rows(ROOT/'data/field/seo2011_layers.csv')
    assert len(plots)==50 and len(trees)==1750
    uncovered=[]
    for p in plots:
        x,y=tm(float(p['경도']),float(p['위도']),5187)
        if not all(any(b['left']<=xx<=b['right'] and b['bottom']<=yy<=b['top'] for b in rows)
                   for xx,yy in ((x-60,y-60),(x-60,y+60),(x+60,y-60),(x+60,y+60))):uncovered.append(p['plot_id'])
    print('60m crop corner coverage failures:',uncovered)
    assert not uncovered
    assert len({r['plot_id'] for r in layers if r['기호']=='T1'})==50
    calc=KoreanAllometry(read_csv_rows(ROOT/'data/reference/nifos2014_allometry.csv'),read_csv_rows(ROOT/'data/reference/nifos2014_carbon_factors.csv'),read_csv_rows(ROOT/'data/reference/nifos2010_carbon_fraction.csv'))
    assert calc.mapping('소나무',40)[0]==[('중부지방소나무',1.)]
    assert calc.mapping('소나무',40.1)[0]==[('강원지방소나무',1.)]
    assert not calc.mapping('등나무',20)[0]
    assert len(calc.mapping('떡갈나무',20)[0])==4
    mass=calc.predict('떡갈나무',20)['AGB_kg']
    assert math.isclose(mass,sum(calc.predict(s,20)['AGB_kg'] for s in calc.oaks)/4)
    cf_mass=sum(calc.predict(s,20)['carbon_cf_kg'] for s in calc.oaks)/4
    assert math.isclose(calc.predict('떡갈나무',20)['carbon_cf_kg'],cf_mass)
    assert calc.predict('곰솔',48)['range_flag']=='곰솔'
    assert math.isclose(training_ratio([10,30],[5,15]),2)
    dimensional_checks()
    outputs=[]
    for threshold in (6,10):
        used=[t for t in trees if t['수종']!='등나무' and float(t['흉고직경_cm'])>=threshold]
        carbon=sum(calc.predict(t['수종'],float(t['흉고직경_cm']))['carbon_kg'] for t in used)
        outputs.append(dict(threshold_cm=threshold,n_trees=len(used),aboveground_carbon_tC=carbon/1000,mean_tC_ha=carbon/1000/2))
    notebook=HERE/'GEUMJEONGSAN_DeepForest_2021.ipynb'
    if notebook.exists():
        nb=json.loads(notebook.read_text(encoding='utf-8'));assert nb['nbformat']==4
        for i,c in enumerate(nb['cells']):
            if c['cell_type']!='code':continue
            assert nb['cells'][i-1]['cell_type']=='markdown' and not c['outputs'] and c['execution_count'] is None
            text=''.join(c['source'])
            if not text.startswith('%'):ast.parse(text)
        assert "files.download('/content/geumjeongsan_results.zip')" in ''.join(nb['cells'][-1]['source'])
        assert len({c['id'] for c in nb['cells']})==len(nb['cells'])
        sources=[''.join(c['source']) for c in nb['cells'] if c['cell_type']=='code']
        combined='\n'.join(s for s in sources if not s.startswith('%'))
        module=ast.parse(combined)
        embedded=next(n.value.value for n in module.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='NEON_PIPELINE_SOURCE' for t in n.targets))
        assert embedded==(HERE.parents[1]/'neon_grsm_ornl/notebooks/neon_pipeline.py').read_text(encoding='utf-8')
        model=next(ast.literal_eval(n.value) for n in module.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='HF_MODEL_PROVENANCE' for t in n.targets))
        assert model['resolved_revision']=='cc21436bc5d572dde8ff5f93c1e71a32f563cace'
        assert combined.index('GEOREF_CHECK=georeference()')<combined.index('DETECTED=estimate_boxes()')
        for required in ['GroupKFold(n_splits=5)','Resampling.bilinear','NEON.predict_crowns','C_calibrated_OOF_tC_ha','reference_proxy_for_circle','forest_fallback']:
            assert required in combined
        rt=ast.parse((HERE/'geumjeongsan_runtime.py').read_text(encoding='utf-8'))
        assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='NEON' and n.func.attr in ('estimate_dbh','estimate_agb','match_one_to_one') for n in ast.walk(rt))
    assert not any(p in sys.modules for p in ('pandas','geopandas','rasterio','pyogrio','numpy','pyproj'))
    result=dict(status='PASS',TIFF_XML_headers=9,control_5179=next(r for r in rows if r['sheet']=='35913003'),
                field_reference=outputs,scientific_packages_imported=False,local_only='stdlib; 원본픽셀/추론/원격다운로드 미실행')
    (HERE/'GEUMJEONGSAN_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
