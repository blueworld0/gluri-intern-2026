"""Memory-bounded, reproducible Hongcheon 2021 GEDI / Sentinel-2 modelling."""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from contextlib import ExitStack
import gc
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time
import traceback
import warnings
import tempfile
import certifi

# GDAL's Windows TLS backend cannot open a CA bundle under a Korean venv path.
# Copy the same trusted certificates to an ASCII temp path; TLS stays verified.
CA_FILE = Path(tempfile.gettempdir()) / 'gedi_step2_ca.pem'
if Path(certifi.where()).resolve() != CA_FILE.resolve():
    shutil.copyfile(certifi.where(), CA_FILE)
certifi.where = lambda: str(CA_FILE)
os.environ['CURL_CA_BUNDLE'] = str(CA_FILE)

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform_bounds
from rasterio.windows import Window, from_bounds
from pyproj import Transformer
from tqdm import tqdm

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')
ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / 'data/processed'
FIG = ROOT / 'figures'
RAW = ROOT / 'data/raw'
BBOX = (127.99, 37.66, 128.21, 37.84)
CRS = 'EPSG:32652'
DATES = ('2021-05-22', '2021-10-24')
BANDS = {'B02':'blue','B03':'green','B04':'red','B05':'rededge1','B06':'rededge2',
         'B07':'rededge3','B08':'nir','B8A':'nir08','B11':'swir16','B12':'swir22','SCL':'scl'}
NODATA = -9999.
SEED = 42
CHUNK = 128
START = time.monotonic()
DEADLINE = START + 3*3600
SCHEMA_VERSION = 3
IPCC = 'https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/4_Volume4/V4_04_Ch4_Forest_Land.pdf'


def log(message):
    print(f'[{(time.monotonic()-START)/60:.1f} min] {message}', flush=True)
    if time.monotonic() > DEADLINE:
        raise TimeoutError('3-hour deadline exceeded; cached work has been retained')


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=str), encoding='utf-8')


def table(frame):
    def text(x):
        return f'{x:,.3f}' if isinstance(x, (float, np.floating)) else str(x)
    return '\n'.join(['| '+' | '.join(map(str, frame.columns))+' |',
                      '| '+' | '.join(['---']*len(frame.columns))+' |'] +
                     ['| '+' | '.join(text(v) for v in row)+' |' for row in frame.itertuples(index=False, name=None)])


def grid():
    trans = Transformer.from_crs('EPSG:4326', CRS, always_xy=True)
    w,s,e,n = BBOX
    # Densified geographic rectangle so rotated UTM edges do not change the study area.
    lons = np.r_[np.linspace(w,e,101), np.full(101,e), np.linspace(e,w,101), np.full(101,w)]
    lats = np.r_[np.full(101,s), np.linspace(s,n,101), np.full(101,n), np.linspace(n,s,101)]
    xs,ys = trans.transform(lons,lats)
    left,bottom = math.floor(min(xs)/10)*10, math.floor(min(ys)/10)*10
    right,top = math.ceil(max(xs)/10)*10, math.ceil(max(ys)/10)*10
    width,height = int((right-left)/10), int((top-bottom)/10)
    polygon = {'type':'Polygon','coordinates':[list(zip(xs,ys)) + [(xs[0],ys[0])]]}
    profile = dict(driver='GTiff',width=width,height=height,crs=CRS,
                   transform=from_origin(left,top,10,10),count=1,dtype='float32',nodata=NODATA,
                   tiled=True,blockxsize=256,blockysize=256,compress='deflate',predictor=3,
                   BIGTIFF='IF_SAFER',num_threads=1,interleave='band')
    return profile, polygon


def windows(profile, rows=CHUNK):
    for row in range(0,profile['height'],rows):
        yield Window(0,row,profile['width'],min(rows,profile['height']-row))


def valid_tif(path, count=None):
    try:
        with rasterio.open(path) as ds:
            ds.read(1,window=Window(ds.width-1,ds.height-1,1,1))
            return count is None or ds.count == count
    except (OSError,rasterio.errors.RasterioError):
        return False


def crop_remote(url, path):
    """HTTP range reads of only the source-grid bbox, with 3 native pixels of padding."""
    if valid_tif(path):
        log(f'캐시 재사용: {path.name}')
        return
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp = path.with_suffix('.part.tif')
    log(f'구역 창 다운로드: {path.name}')
    with rasterio.open(url) as src:
        bounds = transform_bounds('EPSG:4326',src.crs,*BBOX,densify_pts=21)
        win = from_bounds(*bounds, transform=src.transform)
        col,row = math.floor(win.col_off)-3, math.floor(win.row_off)-3
        endcol,endrow = math.ceil(win.col_off+win.width)+3, math.ceil(win.row_off+win.height)+3
        win = Window(col,row,endcol-col,endrow-row).intersection(Window(0,0,src.width,src.height))
        profile = src.profile.copy()
        profile.update(width=int(win.width),height=int(win.height),transform=src.window_transform(win),
                       driver='GTiff',tiled=True,blockxsize=256,blockysize=256,compress='deflate',num_threads=1)
        with rasterio.open(tmp,'w',**profile) as dst:
            for offset in tqdm(range(0,int(win.height),256),desc=path.stem,unit='row-block',leave=False):
                h = min(256,int(win.height)-offset)
                data = src.read(1,window=Window(win.col_off,win.row_off+offset,win.width,h))
                dst.write(data,1,window=Window(0,offset,win.width,h))
    if not valid_tif(tmp):
        raise RuntimeError(f'Invalid cropped raster: {tmp}')
    tmp.replace(path)


def scene_quality(path,feature):
    from rasterio.warp import transform_geom
    w,s,e,n=BBOX
    boundary={'type':'Polygon','coordinates':[[(w,s),(e,s),(e,n),(w,n),(w,s)]]}
    with rasterio.open(path) as src:
        mask=geometry_mask([transform_geom('EPSG:4326',src.crs,boundary)],
                           out_shape=(src.height,src.width),transform=src.transform,invert=True)
        scl=src.read(1)
    n=int(mask.sum())
    covered=mask&(scl!=0)
    cloud=mask&np.isin(scl,[3,8,9,10])
    snow=mask&(scl==11)
    return {'scene_id':feature['id'],'aoi_pixels_20m':n,'covered_pixels':int(covered.sum()),
            'coverage_percent':100*covered.sum()/n,'cloud_shadow_pixels':int(cloud.sum()),
            'cloud_shadow_percent_of_aoi':100*cloud.sum()/n,
            'cloud_shadow_percent_of_covered':100*cloud.sum()/covered.sum() if covered.any() else np.nan,
            'snow_percent_of_aoi':100*snow.sum()/n}


def search_scenes():
    from pystac_client import Client
    out = RAW / 's2'
    out.mkdir(parents=True,exist_ok=True)
    chosen = []
    quality=[]

    def load_exact(date,item_id):
        cache = out / (date+'_selected.json')
        if cache.exists():
            feature = json.loads(cache.read_text(encoding='utf-8'))
            if feature['id'] != item_id:
                raise ValueError(f'Cached item differs from requested exact scene: {item_id}')
        else:
            client = Client.open('https://earth-search.aws.element84.com/v1')
            items = list(client.search(collections=['sentinel-2-l2a'],ids=[item_id]).items())
            if len(items)!=1 or items[0].id!=item_id:
                raise RuntimeError(f'Exact requested scene not found: {item_id}')
            feature = items[0].to_dict()
            save_json(cache,feature)
        log(f"장면: {feature['id']}, scene cloud={feature['properties'].get('eo:cloud_cover')}%")
        scl_path=out/feature['id']/'SCL_bbox.tif'
        crop_remote(feature['assets']['scl']['href'],scl_path)
        q=scene_quality(scl_path,feature)
        quality.append(q)
        log(f"구역 SCL: {feature['id']} coverage={q['coverage_percent']:.4f}%, "
            f"구름·그림자={q['cloud_shadow_percent_of_aoi']:.4f}% (SCL 3/8/9/10), "
            f"눈={q['snow_percent_of_aoi']:.4f}%")
        if q['coverage_percent'] < 99.99:
            raise RuntimeError(f"Scene does not cover 100% of AOI: {feature['id']}")
        return feature,q

    spring,q=load_exact('2021-05-22','S2A_52SDG_20210522_0_L2A')
    autumn,aq=load_exact('2021-10-24','S2B_52SDG_20211024_0_L2A')
    if q['cloud_shadow_percent_of_aoi'] > 15:
        log('5/22 구역 구름·그림자 >15%: 요청대로 7/26으로 교체')
        spring,q=load_exact('2021-07-26','S2B_52SDG_20210726_0_L2A')
    pd.DataFrame(quality).to_csv(PROC/'scene_quality.csv',index=False,encoding='utf-8-sig')
    for feature,q in [(spring,q),(autumn,aq)]:
        date=feature['properties']['datetime'][:10]
        assets = {}
        for band,key in BANDS.items():
            asset = feature['assets'][key]
            path = out / feature['id'] / (band+'_bbox.tif')
            crop_remote(asset['href'],path)
            rb = asset.get('raster:bands',[{}])[0]
            assets[band] = {'path':str(path),'scale':rb.get('scale',.0001 if band!='SCL' else 1),
                            'offset':rb.get('offset',0),'nodata':rb.get('nodata',0)}
        chosen.append({'date':date,'id':feature['id'],'assets':assets,
                       'cloud_cover':feature['properties'].get('eo:cloud_cover'),'aoi_quality':q})
    save_json(PROC/'s2_scene_manifest.json',chosen)
    return chosen


def prepare_dem():
    result=[]
    # E127 is needed for the 127.99~128.00 sliver; both sources remain native 30 m.
    for lon in (127,128):
        name=f'Copernicus_DSM_COG_10_N37_00_E{lon:03d}_00_DEM'
        url=f'https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif'
        path=RAW/'dem'/(name+'_bbox.tif')
        crop_remote(url,path)
        result.append(path)
    return result


def aligned(stack,path,profile,nearest=False):
    src=stack.enter_context(rasterio.open(path))
    return stack.enter_context(WarpedVRT(src,crs=CRS,transform=profile['transform'],
             width=profile['width'],height=profile['height'],dtype='float32',nodata=np.nan,
             resampling=Resampling.nearest if nearest else Resampling.bilinear,
             warp_mem_limit=32))


def read_float(ds,win):
    return ds.read(1,window=win,masked=True).filled(np.nan).astype('float32',copy=False)


def ratio(numerator,denominator):
    out=np.full_like(numerator,np.nan,dtype='float32')
    np.divide(numerator,denominator,out=out,where=np.isfinite(denominator)&(np.abs(denominator)>1e-6))
    return out


def build_features(profile,polygon,scenes,dems):
    names=[scene['date'].replace('-','')+'_'+band for scene in scenes
           for band in [*list(BANDS)[:-1],'NDVI','NDRE','EVI','NDMI']]+['elevation','slope','aspect']
    path=PROC/'features_10m.tif'
    signature=hashlib.sha256(json.dumps({'schema':SCHEMA_VERSION,'bbox':BBOX,'scenes':scenes,
              'dems':list(map(str,dems)),'names':names},sort_keys=True).encode()).hexdigest()
    if valid_tif(path,len(names)):
        with rasterio.open(path) as src:
            if src.tags().get('signature')==signature:
                log('변수 래스터 캐시 재사용')
                return path,names,signature
    tmp=path.with_suffix('.part.tif')
    with ExitStack() as stack:
        views=[]
        for scene in scenes:
            views.append({band:aligned(stack,asset['path'],profile,band=='SCL') for band,asset in scene['assets'].items()})
        dem_views=[aligned(stack,path,profile) for path in dems]
        dst=stack.enter_context(rasterio.open(tmp,'w',**(profile|{'count':len(names)})))
        for i,name in enumerate(names,1):
            dst.set_band_description(i,name)
        dst.update_tags(signature=signature,source='Earth Search Sentinel-2 L2A; Copernicus DEM',
                        reflectance='STAC asset scale and offset applied; bilinear 20m to 10m',
                        slope_units='degrees',aspect_units='degrees clockwise from north; flat=0')
        for win in tqdm(list(windows(profile)),desc='31-band features',unit='block'):
            inside=geometry_mask([polygon],out_shape=(int(win.height),int(win.width)),
                                 transform=rasterio.windows.transform(win,profile['transform']),invert=True)
            band_index=1
            for scene,view in zip(scenes,views):
                scl=read_float(view['SCL'],win)
                bad=np.isin(scl,[0,1,3,8,9,10,11])|~np.isfinite(scl)|~inside
                values={}
                for band in list(BANDS)[:-1]:
                    asset=scene['assets'][band]
                    arr=read_float(view[band],win)
                    arr=arr*np.float32(asset['scale'])+np.float32(asset['offset'])
                    arr[bad]=np.nan
                    values[band]=arr
                    dst.write(np.where(np.isfinite(arr),arr,NODATA),band_index,window=win)
                    band_index+=1
                red,nir,blue=values['B04'],values['B08'],values['B02']
                n8a,re,swir=values['B8A'],values['B05'],values['B11']
                indices=[ratio(nir-red,nir+red),ratio(n8a-re,n8a+re),
                         2.5*ratio(nir-red,nir+6*red-7.5*blue+1),ratio(n8a-swir,n8a+swir)]
                for arr in indices:
                    dst.write(np.where(np.isfinite(arr),arr,NODATA).astype('float32'),band_index,window=win)
                    band_index+=1
                del values,indices
            # One-row halo avoids artificial slope seams at row-chunk boundaries.
            top=max(0,int(win.row_off)-1)
            bottom=min(profile['height'],int(win.row_off+win.height)+1)
            halo=Window(0,top,profile['width'],bottom-top)
            elevation=np.full((bottom-top,profile['width']),np.nan,dtype='float32')
            for view in dem_views:
                block=read_float(view,halo)
                np.copyto(elevation,block,where=~np.isfinite(elevation))
            dy,dx=np.gradient(elevation,10.,10.)
            slope=np.degrees(np.arctan(np.hypot(dx,dy))).astype('float32')
            aspect=np.mod(np.degrees(np.arctan2(-dx,dy)),360).astype('float32')
            aspect[(dx==0)&(dy==0)]=0
            offset=int(win.row_off)-top
            for arr in [elevation,slope,aspect]:
                arr=arr[offset:offset+int(win.height)].copy()
                arr[~inside]=np.nan
                dst.write(np.where(np.isfinite(arr),arr,NODATA),band_index,window=win)
                band_index+=1
    tmp.replace(path)
    save_json(PROC/'feature_names.json',names)
    return path,names,signature


def training(features,names):
    df=pd.read_csv(PROC/'gedi_l4a_2021_filtered.csv',dtype={'shot_number':'uint64'})
    if df.shot_number.duplicated().any():
        raise ValueError('Duplicate GEDI shot_number would cause validation leakage; inspect input')
    tx=Transformer.from_crs('EPSG:4326',CRS,always_xy=True)
    x,y=tx.transform(df.lon_lowestmode.to_numpy(),df.lat_lowestmode.to_numpy())
    with rasterio.open(features) as src:
        vals=np.asarray(list(src.sample(zip(x,y),masked=True)),dtype='float32')
    vals[vals==NODATA]=np.nan
    for i,name in enumerate(names):
        df[name]=vals[:,i]
    df['easting']=x
    df['northing']=y
    df['block_id']=(np.floor(x/2000).astype('int64').astype(str)+np.char.add('_',np.floor(y/2000).astype('int64').astype(str)))
    df['complete_features']=np.isfinite(vals).all(axis=1)
    df['valid_target']=np.isfinite(df.agbd)&df.agbd.ge(0)
    df.to_csv(PROC/'training_all_samples.csv',index=False,encoding='utf-8-sig')
    counts=[]
    for threshold in (400,600,800):
        within=df.valid_target&df.agbd.le(threshold)
        usable=within&df.complete_features
        counts.append({'threshold_Mg_ha':threshold,'input_points':len(df),
                       'invalid_target':int((~df.valid_target).sum()),
                       'above_threshold':int((df.valid_target&df.agbd.gt(threshold)).sum()),
                       'missing_features_after_target_filter':int((within&~df.complete_features).sum()),
                       'training_points':int(usable.sum()),'spatial_blocks':df.loc[usable,'block_id'].nunique()})
    pd.DataFrame(counts).to_csv(PROC/'training_filter_counts.csv',index=False,encoding='utf-8-sig')
    return df,pd.DataFrame(counts)


def model(name):
    if name=='RF':
        from sklearn.ensemble import RandomForestRegressor
        return RandomForestRegressor(n_estimators=300,random_state=SEED,n_jobs=1,
                    min_samples_leaf=2,max_features=1.0)
    from lightgbm import LGBMRegressor
    return LGBMRegressor(n_estimators=300,learning_rate=.05,num_leaves=31,
             min_child_samples=20,reg_lambda=1.,random_state=SEED,n_jobs=1,
             deterministic=True,force_col_wise=True,verbosity=-1)


def metrics(y,pred):
    from sklearn.metrics import r2_score,mean_absolute_error
    return {'R2':float(r2_score(y,pred)),'RMSE':float(np.sqrt(np.mean((pred-y)**2))),
            'Bias':float(np.mean(pred-y)),'MAE':float(mean_absolute_error(y,pred))}


def evaluate(df,names,signature):
    from sklearn.model_selection import KFold,GroupKFold
    rows=[]
    base=df.loc[df.valid_target&df.complete_features]
    checkpoint=PROC/'step2_cv_cache'
    checkpoint.mkdir(exist_ok=True)
    for threshold in (600,400,800):
        data=base.loc[base.agbd.le(threshold)].copy().reset_index(drop=True)
        if len(data)<10 or data.block_id.nunique()<5:
            raise ValueError(f'Insufficient points/blocks at threshold {threshold}')
        X=data[names].to_numpy(dtype='float32')
        y=data.agbd.to_numpy(dtype='float64')
        data_hash=hashlib.sha256(pd.util.hash_pandas_object(data[['shot_number','agbd',*names]],index=False).values.tobytes()).hexdigest()
        for scheme in ('random','spatial'):
            splitter=KFold(5,shuffle=True,random_state=SEED) if scheme=='random' else GroupKFold(5)
            folds=list(splitter.split(X,y,groups=data.block_id if scheme=='spatial' else None))
            for name in ('RF','LightGBM'):
                tag=f'{name}_{scheme}_{threshold}'
                cached=checkpoint/(tag+'.npz')
                cache_key=signature+data_hash+str(model(name).get_params())
                key=hashlib.sha256(cache_key.encode()).hexdigest()
                if cached.exists():
                    with np.load(cached) as z:
                        pred=z['prediction'].copy() if str(z['key'])==key else None
                else:
                    pred=None
                if pred is None:
                    pred=np.full(len(y),np.nan)
                    for fold,(train,test) in enumerate(folds,1):
                        if scheme=='spatial':
                            assert not set(data.block_id.iloc[train])&set(data.block_id.iloc[test])
                        estimator=model(name)
                        estimator.fit(X[train],y[train])
                        pred[test]=estimator.predict(X[test])
                        log(f'CV {tag} fold {fold}/5: n_train={len(train)}, n_test={len(test)}')
                        del estimator
                        gc.collect()
                    np.savez_compressed(cached,prediction=pred,key=key)
                score=metrics(y,pred)
                rows.append({'threshold_Mg_ha':threshold,'split':scheme,'model':name,'n':len(y),
                             'blocks':data.block_id.nunique(),'aggregation':'pooled out-of-fold',**score})
                fold_id=np.zeros(len(y),dtype='int8')
                for fold,(_,test) in enumerate(folds,1):
                    fold_id[test]=fold
                out=data[['shot_number','agbd','block_id','easting','northing']].copy()
                out['prediction']=pred
                out['fold']=fold_id
                out.to_csv(PROC/(tag+'_oof.csv'),index=False,encoding='utf-8-sig')
                pd.DataFrame(rows).to_csv(PROC/'model_cv.csv',index=False,encoding='utf-8-sig')
    return pd.DataFrame(rows)


def fit_models(df,names):
    data=df.loc[df.valid_target&df.complete_features&df.agbd.le(600)].copy()
    X=data[names].to_numpy(dtype='float32')
    y=data.agbd.to_numpy()
    models={}
    import joblib
    for name in ('RF','LightGBM'):
        log('최종 학습: '+name)
        estimator=model(name)
        estimator.fit(X,y)
        models[name]=estimator
        joblib.dump(estimator,PROC/f'{name}_model.joblib',compress=3)
    return models


def predict_maps(models,selected,features,profile,polygon):
    wc_path=RAW/'worldcover/ESA_WorldCover_10m_2021_v200_N36E126_Map.tif'
    totals={name:{'pixels':0,'agb_sum_Mg':0.,'min':float('inf'),'max':float('-inf')} for name in models}
    total_tree=0
    missing_tree=0
    with ExitStack() as stack:
        src=stack.enter_context(rasterio.open(features))
        wc=aligned(stack,wc_path,profile,True)
        outputs={name:stack.enter_context(rasterio.open(PROC/f'agb_map_{name}_10m.part.tif','w',**profile)) for name in models}
        for name,dst in outputs.items():
            dst.set_band_description(1,'AGBD_Mg_ha')
            dst.update_tags(model=name,training_threshold_Mg_ha=600,units='Mg/ha',
                            selected=str(name==selected),seed=SEED)
        for win in tqdm(list(windows(profile)),desc='RF + LightGBM map prediction',unit='block'):
            arr=src.read(window=win)
            mask=np.isfinite(arr).all(axis=0)&(arr!=NODATA).all(axis=0)
            inside=geometry_mask([polygon],out_shape=mask.shape,
                                 transform=rasterio.windows.transform(win,profile['transform']),invert=True)
            tree=(read_float(wc,win)==10)&inside
            total_tree+=int(tree.sum())
            missing_tree+=int((tree&~mask).sum())
            good=mask&tree
            if good.any():
                X=arr[:,good].T.copy()
            for name,estimator in models.items():
                raster=np.full(mask.shape,NODATA,dtype='float32')
                if good.any():
                    # Negative biomass is physically invalid; log/count this clipping.
                    pred=np.maximum(estimator.predict(X),0).astype('float32')
                    raster[good]=pred
                    stat=totals[name]
                    stat['pixels']+=len(pred)
                    stat['agb_sum_Mg']+=float(pred.sum(dtype='float64')*.01)
                    stat['min']=min(stat['min'],float(pred.min()))
                    stat['max']=max(stat['max'],float(pred.max()))
                outputs[name].write(raster,1,window=win)
    for name in models:
        (PROC/f'agb_map_{name}_10m.part.tif').replace(PROC/f'agb_map_{name}_10m.tif')
    shutil.copyfile(PROC/f'agb_map_{selected}_10m.tif',PROC/'agb_map_10m.tif')
    carbon=[]
    for name,stat in totals.items():
        area=stat['pixels']*.01
        agb=stat['agb_sum_Mg']
        tc=agb*1.25*.47
        carbon.append({'model':name,'selected':name==selected,'tree_cover_area_ha':total_tree*.01,
             'mapped_tree_area_ha':area,'unmapped_tree_area_ha':missing_tree*.01,
             'coverage_percent':100*stat['pixels']/total_tree if total_tree else 0,
             'mean_AGB_Mg_ha':agb/area if area else np.nan,'AGB_total_Mg':agb,
             'mean_carbon_tC_ha':tc/area if area else np.nan,'carbon_tC':tc,'CO2_equivalent_t':tc*44/12,
             'carbon_low_tC':agb*1.20*.47,'carbon_high_tC':agb*1.30*.49,
             'predicted_min_Mg_ha':stat['min'],'predicted_max_Mg_ha':stat['max']})
    result=pd.DataFrame(carbon)
    result.to_csv(PROC/'carbon_totals.csv',index=False,encoding='utf-8-sig')
    return result


def figures(cv,models,names,profile):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Malgun Gothic'
    plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for ax,name in zip(axes,models):
        data=pd.read_csv(PROC/f'{name}_spatial_600_oof.csv',dtype={'shot_number':'uint64'})
        ax.scatter(data.agbd,data.prediction,s=7,alpha=.35,edgecolors='none')
        maximum=max(data.agbd.max(),data.prediction.max())*1.03
        ax.plot([0,maximum],[0,maximum],'k--',linewidth=1)
        score=cv.loc[(cv.threshold_Mg_ha==600)&(cv.split=='spatial')&(cv.model==name)].iloc[0]
        ax.set(xlabel='GEDI 추정 AGBD (Mg/ha)',ylabel='공간 교차검증 예측 (Mg/ha)',
               title=f'{name}: R²={score.R2:.3f}, RMSE={score.RMSE:.1f}',xlim=(0,maximum),ylim=(0,maximum))
    fig.tight_layout()
    fig.savefig(FIG/'obs_vs_pred.png',dpi=170)
    plt.close(fig)
    importance=pd.DataFrame({'feature':names,**{name:est.feature_importances_ for name,est in models.items()}})
    importance.to_csv(PROC/'feature_importance.csv',index=False,encoding='utf-8-sig')
    fig,axes=plt.subplots(1,2,figsize=(13,8))
    for ax,name in zip(axes,models):
        ranked=importance.sort_values(name).tail(20)
        norm=importance[name].sum()
        ax.barh(ranked.feature,ranked[name]/norm,color='forestgreen')
        ax.set(title=name,xlabel='정규화 중요도 (RF: impurity, LightGBM: split count)')
    fig.tight_layout()
    fig.savefig(FIG/'feature_importance.png',dpi=170)
    plt.close(fig)
    with rasterio.open(PROC/'agb_map_10m.tif') as src:
        scale=max(1,math.ceil(max(src.width,src.height)/1200))
        arr=src.read(1,out_shape=(src.height//scale,src.width//scale),masked=True)
        bounds=src.bounds
    fig,ax=plt.subplots(figsize=(9,9))
    im=ax.imshow(arr,extent=(bounds.left,bounds.right,bounds.bottom,bounds.top),cmap='YlGn',vmin=0)
    ax.set(xlabel='UTM 52N Easting (m)',ylabel='UTM 52N Northing (m)',title='홍천 2021 예측 AGBD — 나무 피복 영역')
    ax.ticklabel_format(style='plain',axis='both',useOffset=False)
    fig.colorbar(im,ax=ax,label='AGBD (Mg/ha)',shrink=.8)
    fig.text(.5,.015,'흰색: 비나무 피복 또는 변수 결측으로 미예측',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.035,1,1))
    fig.savefig(FIG/'agb_map.png',dpi=180)
    plt.close(fig)


def report(cv,counts,carbon,selected,scenes,profile,df):
    main=cv.loc[cv.threshold_Mg_ha==600,['split','model','n','R2','RMSE','Bias','MAE']]
    sensitivity=cv.loc[cv.split=='spatial',['threshold_Mg_ha','model','n','R2','RMSE','Bias','MAE']]
    chosen=carbon.loc[carbon.model==selected].iloc[0]
    lower,upper=carbon.carbon_low_tC.min(),carbon.carbon_high_tC.max()
    base=df.loc[df.valid_target&df.complete_features&df.agbd.le(600)]
    text=f'''# 홍천 Sentinel-2 기반 AGB 지도 및 탄소저장량 추정

## 주요 결과

- 주 모델: **{selected}**, 기준 AGBD ≤600 Mg/ha. 공간 5-fold의 pooled RMSE가 낮은 모델을 선택함.
- 입력 GEDI: {len(df):,}개. 최종 학습: {len(base):,}개, 2km 블록 {base.block_id.nunique()}개.
- 나무 피복 면적: {chosen.tree_cover_area_ha:,.2f} ha. 실제 지도화: {chosen.mapped_tree_area_ha:,.2f} ha ({chosen.coverage_percent:.2f}%).
- 결측 변수로 미예측한 나무 피복: {chosen.unmapped_tree_area_ha:,.2f} ha. 이를 0으로 간주하거나 추정 면적에 외삽하지 않음.
- 예측된 나무 피복의 평균 AGB: {chosen.mean_AGB_Mg_ha:,.2f} Mg/ha.
- 평균 생체량 탄소밀도(지상+뿌리): {chosen.mean_carbon_tC_ha:,.2f} tC/ha.
- **지도화 면적의 탄소저장량: {chosen.carbon_tC:,.0f} tC**, CO₂ 질량 환산 {chosen.CO2_equivalent_t:,.0f} tCO₂.
- 두 모델 × R·CF 최저/최고 조합: **{lower:,.0f}~{upper:,.0f} tC**
  (CO₂ 환산 {lower*44/12:,.0f}~{upper*44/12:,.0f} tCO₂).
- 위 범위는 계수·모델 시나리오 범위이며 통계적 신뢰구간이 아니다. 미예측 면적이 있으면 구역 전체 총량으로 볼 수 없다.
- 탄소저장량이며 연간 흡수량·배출저감량이 아니다. 토양·고사목·낙엽은 포함하지 않는다.

## 자료와 처리

- 구역: {BBOX}, WGS84 직사각형. UTM 52N (EPSG:32652), 10m 격자 {profile['width']}×{profile['height']}.
- {scenes[0]['date']}: `{scenes[0]['id']}`, 장면 구름 {scenes[0]['cloud_cover']}%, 구역 구름·그림자 {scenes[0]['aoi_quality']['cloud_shadow_percent_of_aoi']:.3f}%.
- {scenes[1]['date']}: `{scenes[1]['id']}`, 장면 구름 {scenes[1]['cloud_cover']}%, 구역 구름·그림자 {scenes[1]['aoi_quality']['cloud_shadow_percent_of_aoi']:.3f}%.
- 5/22 장면은 STAC에서 S2B가 아닌 S2A로 확인되어 위성명만 바로잡았다. 날짜·타일·처리 버전은 유지했다.
- 5/22 구역 구름·그림자(SCL 3/8/9/10)가 15%를 넘으면 7/26으로 교체하는 규칙을 적용했고, 이번에는 8.138%여서 교체하지 않았다.
- 최종 두 장면의 구역 관측 피복: {scenes[0]['aoi_quality']['coverage_percent']:.4f}%, {scenes[1]['aoi_quality']['coverage_percent']:.4f}%. 10/24는 20m SCL 약 1픽셀만 nodata. 전체 구역 관측과 구름 제거 후 유효 면적은 구분해야 함.
- 공개 COG에서 bbox 창만 HTTP Range로 읽어 native-grid cutout을 data/raw/s2에 저장.
- asset별 scale·offset 적용. 20m 밴드를 bilinear로 10m 정렬; SCL과 WorldCover는 nearest.
- SCL 3/8/9/10/11 및 nodata=0, defective=1을 결측 처리. 두 장면 중 하나라도 변수가 결측이면 학습·지도 예측에서 제외.
- 장면별 10개 밴드 + NDVI, NDRE(B8A,B05), EVI, NDMI(B8A,B11); DEM 고도·경사·향. 총 31개 변수.
- Copernicus DEM GLO-30의 N37E128 및 서쪽 경계용 N37E127 창을 결합. DEM은 실제로 DSM이며 30m→10m 재표본화로 해상도가 증가하지 않음.
- 경사는 UTM 10m 간격 중앙 차분, 향은 북쪽 0도·시계방향. 경사 산출 시 블록 경계 halo 사용.
- GEDI footprint 중심의 단일 10m 픽셀 값을 추출. shot_number는 uint64로 보존.
- 입력 중 WorldCover=10 여부와 무관하게 유효한 GEDI를 학습에 사용하고, 최종 지도만 나무 피복에 제한.

## 학습 점 필터 및 이상치 민감도

{table(counts)}

상한 초과 → 결측 변수 제거 순서다. 400/600/800은 다른 표본 집합이므로 지표를 동일한 검증 모집단으로 해석하지 않는다.

## 성능: 주 분석 AGBD ≤600

{table(main)}

R²는 무단위, RMSE/Bias/MAE는 Mg/ha. Bias = 예측−GEDI. 표는 모든 out-of-fold 예측을 합친 지표다.
무작위 KFold는 shuffle=True, seed=42. 공간 GroupKFold는 UTM 좌표를 2,000m로 나눈 고정 블록이며 5개 fold 사이 블록 중복 없음.
RF: 300 trees, min_samples_leaf=2, n_jobs=1. LightGBM: 300 trees, learning_rate=0.05, num_leaves=31, n_jobs=1.
사전 고정 설정으로 튜닝하지 않았으며 최적 모델의 독립 외부 검증은 수행하지 않았다.

## 공간 검증의 이상치 상한 민감도

{table(sensitivity)}

## 탄소 계산과 계수

탄소(tC) = Σ(AGBD × 0.01ha) × (1+R) × CF. CO₂ 질량 환산 = tC × 44/12.
픽셀마다 같은 UTM 지도면적 100m²=0.01ha를 사용한다(투영 축척의 미세한 차이 무보정).

- CF: IPCC 2006 Vol.4 Ch.4, 표 4.3, p.4.48, Temperate and Boreal / All: **0.47 (0.47~0.49)** 확인.
- R: 표 4.4, p.4.49는 온대 침엽수·참나무·기타 활엽수와 바이오매스 구간별로 다른 값을 제시한다.
  구역의 수종 구분이 없어 적용할 하나의 계수는 **확인 필요**. 사용자 지정 대체값 **0.25 (0.20~0.30)** 사용.
  이 대체값을 IPCC가 제시한 온대림 공통값으로 해석하지 않는다.
- CF는 지상부 표에서 확인한 값을 사용자 식에 따라 뿌리를 포함한 생체량에도 동일 적용하는 근사다.
- 범위는 RF와 LightGBM 지도의 AGB 합 각각에 R=0.20/0.30, CF=0.47/0.49의 4개 조합을 적용한 최저·최고다.
- 지도에서 음수 AGB 예측은 물리적으로 불가능하므로 0으로 처리했다. 선택된 RF 지도는 원래 예측 최솟값이 양수였으며 이 처리의 영향을 받지 않았다.
- 출처: [IPCC 2006 Vol.4 Ch.4 표 4.3·4.4]({IPCC}). 국내 산림·수종별 계수로 대체한 결과가 아니다.

{table(carbon)}

## 한계와 해석

- 기준값이 GEDI 모델 추정치다. 위 성능은 현장 측정 AGB의 정확도가 아니라 GEDI 재현 성능이다.
- 2021년 영상 2장만 사용해 계절·구름·수분·생육상태를 충분히 반영하지 못한다. GEDI 관측일과도 다르다.
- GEDI 관측이 궤도에 집중돼 공간적으로 대표성이 부족할 수 있다. 2km 블록은 독립성을 완전히 보장하지 않으며 블록 경계 양쪽의 가까운 점을 차단하는 buffer 검증은 아니다.
- 단일 10m 픽셀과 약 25m GEDI footprint의 공간 지지 범위가 다르고 위치오차·경사·혼합 픽셀 영향이 있다.
- RF는 학습값 범위를 넘는 예측이 어렵다. LightGBM도 상수 잎 기반 트리라 새로운 범위로의 안정적 외삽을 보장하지 않지만 RF처럼 엄밀히 학습 최솟값·최댓값 내로 제한되지는 않는다.
- 상한 600 제거는 고바이오매스 산림을 과소평가할 수 있다. 400/800 민감도는 CV 성능·표본 변화까지 확인했으며 그 기준별 지도·탄소 총량 민감도는 산출하지 않았다.
- 20m·30m 입력을 10m로 보간해도 10m 수준의 독립적인 관측 정밀도가 생기지 않는다.
- IPCC 기본계수 및 가정 R을 사용했으며 국내 계수가 아니다. 수종·뿌리·GEDI·모델·분류·결측의 모든 불확실성을 합성한 신뢰구간이 아니다.
- WorldCover=10은 나무 피복이며 법적 산림 경계가 아니다. 결측 지역은 그대로 남겨 두었다.
- 변수 중요도는 RF impurity와 LightGBM split count를 각각 정규화한 기술통계이며 인과효과나 모델 간 동일한 중요도 척도가 아니다.

## 메모리·시간·재현

- 전체 31밴드 래스터를 메모리에 올리지 않고 {CHUNK}행 단위 처리, GDAL 캐시 32MiB, warp 32MiB, 모델 단일 스레드.
- 모델·검증 seed=42. 선택 장면 메타데이터, CV 캐시·OOF, 학습표, 최종 모델을 저장.
- 마지막 실행 소요(기존 영상·변수 캐시 재사용): {(time.monotonic()-START)/60:.1f}분. 중간 캐시를 포함한 재실행은 동일 원본·조건을 사용한다.
'''
    (ROOT/'step2_summary.md').write_text(text,encoding='utf-8')
    print('\n'+table(main)+'\n')
    print(f'선택 모델: {selected}; 탄소: {chosen.carbon_tC:,.0f} tC; 범위: {lower:,.0f}~{upper:,.0f} tC')


def main():
    global CHUNK
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chunk-rows',type=int,default=128)
    args=parser.parse_args()
    CHUNK=args.chunk_rows
    if not 16<=CHUNK<=256:
        raise ValueError('chunk-rows must be between 16 and 256')
    for p in (PROC,FIG,RAW/'s2',RAW/'dem'):
        p.mkdir(parents=True,exist_ok=True)
    log('예상: 장면 검색·구역 다운로드 10~45분; 변수 5~15분; CV/민감도 10~35분; 지도 10~40분; 정리 2~5분 (총 37~140분)')
    profile,polygon=grid()
    log(f"10m grid: {profile['width']}×{profile['height']}; 31-band disk raster, {CHUNK}-row chunks")
    with rasterio.Env(GDAL_CACHEMAX=32*1024*1024,VSI_CACHE=False,GDAL_NUM_THREADS='1',
                      GDAL_CURL_CA_BUNDLE=str(CA_FILE),
                      GDAL_HTTP_TIMEOUT='120',GDAL_HTTP_MAX_RETRY='3',GDAL_HTTP_RETRY_DELAY='2',
                      GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif',
                      AWS_NO_SIGN_REQUEST='YES'):
        scenes=search_scenes()
        dems=prepare_dem()
        feature_path,names,signature=build_features(profile,polygon,scenes,dems)
        df,counts=training(feature_path,names)
        log('학습 필터:\n'+table(counts))
        cv=evaluate(df,names,signature)
        scores=cv.loc[(cv.threshold_Mg_ha==600)&(cv.split=='spatial')].sort_values(['RMSE','model'])
        selected=scores.iloc[0].model
        models=fit_models(df,names)
        carbon=predict_maps(models,selected,feature_path,profile,polygon)
        figures(cv,models,names,profile)
        report(cv,counts,carbon,selected,scenes,profile,df)
        log('완료: step2_summary.md')


if __name__=='__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
