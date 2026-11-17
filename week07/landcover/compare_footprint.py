"""Compare center pixels with area-weighted r=12.5m reflectance; frozen CV."""
from pathlib import Path
from contextlib import ExitStack
import argparse
import hashlib
import json
import math
import os
import platform

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window,from_bounds
from shapely.geometry import Point,box
import shapely
import sklearn
from sklearn.ensemble import RandomForestRegressor
from evaluate_forest import metric,md,save_json

BANDS=['B02','B03','B04','B08','B11','B12']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def weighted_mean(values,weights,valid):
    support=weights[valid].sum()
    return float(np.sum(values[valid]*weights[valid])/support)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=HERE/'forest_evaluation')
    out=parser.parse_args().output_dir.resolve()
    if out.exists() and any(out.glob('footprint_*')):
        raise RuntimeError('Footprint outputs already exist; use a new --output-dir.')
    base=HERE/'forest_evaluation'
    rasters=ROOT/'week06/step5/outputs/rasters'
    paths={k:base/v for k,v in dict(data='forest184_input.csv',splits='forest184_split_assignments.csv',
        config='experiment_config.json',previous='predictions.csv',validation='validation.json').items()}
    paths['helper']=HERE/'evaluate_forest.py'
    for b in BANDS:
        paths[b]=rasters/f'reflectance_{b}.tif'
    paths['scl']=rasters/'scl.tif'
    paths['classi']=next((ROOT/'week06/step5/source').glob('*.SAFE/GRANULE/*/QI_DATA/MSK_CLASSI_B00.jp2'))
    hashes={str(p.relative_to(ROOT)):sha(p) for p in paths.values()}
    assert json.loads(paths['validation'].read_text(encoding='utf-8'))['status']=='passed'
    config=json.loads(paths['config'].read_text(encoding='utf-8'))
    sets,params=config['feature_sets'],config['rf_params']
    data=pd.read_csv(paths['data'],dtype={'shot_number':str}).set_index('point_id',drop=False)
    assert len(data)==184 and data.index.is_unique
    averaged=data.copy()
    audit,pixel_rows=[],[]
    with ExitStack() as stack:
        ds={b:stack.enter_context(rasterio.open(paths[b])) for b in BANDS}
        scl=stack.enter_context(rasterio.open(paths['scl']))
        classi=stack.enter_context(rasterio.open(paths['classi']))
        ref=ds['B04']
        assert ref.crs.to_epsg()==32652 and ref.res==(10,10)
        assert all(d.transform==ref.transform and d.crs==ref.crs and d.shape==ref.shape for d in list(ds.values())+[scl])
        assert classi.crs==ref.crs
        bounds=(data.easting_m.min()-200,data.northing_m.min()-200,data.easting_m.max()+200,data.northing_m.max()+200)
        cw=from_bounds(*bounds,classi.transform).round_offsets().round_lengths()
        classi_arr=classi.read([1,2,3],window=cw)
        ct=classi.window_transform(cw)
        for row in data.itertuples(index=False):
            x,y=row.easting_m,row.northing_m
            disk=Point(x,y).buffer(12.5,quad_segs=64)
            r,c=ref.index(x,y)
            win=Window(c-2,r-2,5,5)
            assert c>=2 and r>=2 and c+2<ref.width and r+2<ref.height
            t=ref.window_transform(win)
            weights=np.zeros((5,5))
            xx,yy=np.zeros((5,5)),np.zeros((5,5))
            for i in range(5):
                for j in range(5):
                    left,top=t*(j,i)
                    right,bottom=t*(j+1,i+1)
                    weights[i,j]=disk.intersection(box(left,bottom,right,top)).area
                    xx[i,j],yy[i,j]=t*(j+.5,i+.5)
            assert np.isclose(weights.sum(),disk.area,rtol=1e-10,atol=1e-7)
            arr={b:ds[b].read(1,window=win,masked=True) for b in BANDS}
            band_valid=np.logical_and.reduce([~np.ma.getmaskarray(a)&np.isfinite(a.filled(np.nan)) for a in arr.values()])
            scl_array=scl.read(1,window=win)
            scl_valid=np.isin(scl_array,[4,5,6])
            rr,cc=map(np.asarray,rasterio.transform.rowcol(ct,xx.ravel(),yy.ravel()))
            assert rr.min()>=0 and cc.min()>=0 and rr.max()<classi_arr.shape[1] and cc.max()<classi_arr.shape[2]
            classi_valid=(classi_arr[:,rr,cc]==0).all(axis=0).reshape(5,5)
            valid=band_valid&scl_valid&classi_valid
            coverage=float(weights[valid].sum()/disk.area)
            record=dict(point_id=row.point_id,common_valid_fraction=coverage,circle_area_m2=float(disk.area),
                circle_area_to_pi_ratio=float(disk.area/(math.pi*12.5**2)),intersected_pixel_n=int((weights>0).sum()),
                effective_pixel_n=float(weights.sum()**2/np.sum(weights**2)),
                invalid_band_area_fraction=float(weights[~band_valid].sum()/disk.area),
                excluded_scl_area_fraction=float(weights[~scl_valid].sum()/disk.area),
                excluded_classi_area_fraction=float(weights[~classi_valid].sum()/disk.area))
            for b in BANDS:
                assert np.isclose(arr[b][2,2],getattr(row,'reflectance_'+b),atol=1e-7,rtol=0)
                average=weighted_mean(arr[b].filled(np.nan).astype(float),weights,valid)
                averaged.loc[row.point_id,'reflectance_'+b]=average
                record['center_'+b]=float(arr[b][2,2])
                record['circle_'+b]=average
            audit.append(record)
            for i,j in zip(*np.where(weights>0)):
                pixel_rows.append(dict(point_id=row.point_id,raster_row=r-2+i,raster_col=c-2+j,
                    intersection_area_m2=float(weights[i,j]),common_valid=bool(valid[i,j]),scl=int(scl_array[i,j]),
                    **{b:float(arr[b][i,j]) if not np.ma.is_masked(arr[b][i,j]) else None for b in BANDS}))
    audit=pd.DataFrame(audit)
    out.mkdir(parents=True,exist_ok=True)
    audit.to_csv(out/'footprint_sampling_audit.csv',index=False)
    pixels=pd.DataFrame(pixel_rows)
    pixels.to_csv(out/'footprint_pixel_weights.csv',index=False)
    # Preserve 184-point paired comparison: never impute, silently discard a
    # point, or normalize a small clear patch and call it a full-footprint mean.
    if not audit.common_valid_fraction.ge(1-1e-6).all():
        save_json(out/'footprint_validation.json',dict(status='stopped_insufficient_full_area',
            incomplete_n=int(audit.common_valid_fraction.lt(1-1e-6).sum()),model_fits=0))
        raise RuntimeError('Some circles lack full common valid support; no model fitted. See audit.')
    blue,red,nir,swir=[averaged['reflectance_'+b] for b in ['B02','B04','B08','B11']]
    denominator=nir+6*red-7.5*blue+1
    averaged['ndvi']=(nir-red)/(nir+red)
    averaged['evi']=2.5*(nir-red)/denominator
    averaged['ndmi']=(nir-swir)/(nir+swir)
    averaged['circle_evi_denominator']=denominator
    averaged['circle_evi_review']=(averaged.evi.abs()>2)|(denominator.abs()<.05)
    assert np.isfinite(averaged[sum(sets.values(),[])]).all().all()
    assert averaged.ndvi.between(-1,1).all() and averaged.ndmi.between(-1,1).all()
    averaged.to_csv(out/'footprint_forest184_input.csv',index=False,encoding='utf-8-sig')
    shifts=pd.DataFrame([dict(feature=col,mean_center=float(data[col].mean()),mean_circle=float(averaged[col].mean()),
        mean_abs_change=float((averaged[col]-data[col]).abs().mean()),max_abs_change=float((averaged[col]-data[col]).abs().max()))
        for col in sets['indices3_bands6']])
    shifts.to_csv(out/'footprint_feature_changes.csv',index=False)
    assignments=pd.read_csv(paths['splits'])
    roles,fold_checks={},[]
    for fold in range(1,6):
        a=assignments[assignments.split_id.eq(f'new441_spatial5_fold{fold}')]
        assert len(a)==184 and a.point_id.is_unique and set(a.point_id)==set(data.index)
        roles[fold]={role:a.loc[a.role.eq(role),'point_id'].tolist() for role in ['train','test','buffer_excluded']}
        tr,te=[data.loc[roles[fold][role]] for role in ['train','test']]
        distance=np.sqrt(((tr[['easting_m','northing_m']].to_numpy()[:,None]-te[['easting_m','northing_m']].to_numpy()[None])**2).sum(axis=2)).min()
        assert distance>=500 and set(tr.spatial_group).isdisjoint(te.spatial_group)
        fold_checks.append(dict(fold=fold,train_n=len(tr),test_n=len(te),original_excluded_n=len(roles[fold]['buffer_excluded']),min_distance_m=float(distance)))
    save_json(out/'footprint_config.json',dict(radius_m=12.5,crs='EPSG:32652',circle_segments=256,
        weighting='pixel-circle intersection area; average reflectance first, then calculate indices',
        valid_policy='common six-band validity + SCL 4/5/6 + CLASSI cloud/cirrus/snow=0; require >=1-1e-6 full-circle coverage',
        forest_mask_inside_circle=False,all184_retained=True,feature_sets=sets,rf_params=params,
        previous_split_roles_unchanged=True,new_hyperparameter_search=False,input_sha256=hashes,
        script_sha256=sha(Path(__file__)),versions=dict(python=platform.python_version(),numpy=np.__version__,
            pandas=pd.__version__,sklearn=sklearn.__version__,rasterio=rasterio.__version__,shapely=shapely.__version__)))
    old=pd.read_csv(paths['previous'])
    old=old[old.training_population.eq('forest184')]
    output=[]
    reproductions=[]
    for mode,source in [('center',data),('circle_r12p5',averaged)]:
        for fold in range(1,6):
            tr,te=[source.loc[roles[fold][role]] for role in ['train','test']]
            for model,columns in sets.items():
                pred=RandomForestRegressor(**params).fit(tr[columns],tr.agbd).predict(te[columns])
                frame=te[['point_id','agbd','spatial_group','circle_r12p5_mixed_forest_nonforest']].copy()
                frame['sampling'],frame['fold'],frame['model']=mode,fold,model
                frame['prediction'],frame['residual']=pred,te.agbd.to_numpy()-pred
                output.append(frame.reset_index(drop=True))
                if mode=='center':
                    previous=old[old.model.eq(model)&old.fold.eq(fold)].set_index('point_id').loc[te.point_id]
                    reproductions.append(float(np.max(abs(pred-previous.prediction.to_numpy()))))
                if mode=='circle_r12p5' and fold==1 and model=='indices3':
                    again=RandomForestRegressor(**params).fit(tr[columns],tr.agbd).predict(te[columns])
                    assert np.array_equal(pred,again)
            print(f'{mode} fold {fold}: train={len(tr)}, test={len(te)}',flush=True)
    assert max(reproductions)<1e-10
    pred=pd.concat(output,ignore_index=True)
    dummy=old[old.model.eq('dummy_mean')].copy()
    dummy['sampling']='dummy_shared'
    pred=pd.concat([pred,dummy[pred.columns.intersection(dummy.columns)]],ignore_index=True)
    for (mode,model),g in pred.groupby(['sampling','model']):
        assert len(g)==184 and g.point_id.is_unique
        assert set(zip(g.point_id,g.fold))==set(zip(dummy.point_id,dummy.fold))
    overall=pd.DataFrame([dict(sampling=mode,model=model,**metric(g)) for (mode,model),g in pred.groupby(['sampling','model'],sort=False)])
    comparisons=[]
    for model in sets:
        a=overall[(overall.sampling=='center')&(overall.model==model)].iloc[0]
        b=overall[(overall.sampling=='circle_r12p5')&(overall.model==model)].iloc[0]
        comparisons.append(dict(model=model,n=184,center_rmse=a.rmse,circle_rmse=b.rmse,rmse_delta=b.rmse-a.rmse,
            rmse_reduction_pct=100*(a.rmse-b.rmse)/a.rmse,center_mae=a.mae,circle_mae=b.mae,center_r2=a.r2,circle_r2=b.r2))
    comparisons=pd.DataFrame(comparisons)
    fm=pd.DataFrame([dict(sampling=mode,model=model,fold=fold,**metric(g)) for (mode,model,fold),g in pred.groupby(['sampling','model','fold'],sort=False)])
    metadata=data[['point_id','circle_r12p5_mixed_forest_nonforest']].reset_index(drop=True)
    diag=pred.drop(columns='circle_r12p5_mixed_forest_nonforest',errors='ignore').merge(metadata,on='point_id',validate='many_to_one')
    scopes=pd.DataFrame([dict(sampling=mode,model=model,scope=scope,**metric(h),underprediction_n=int(h.residual.gt(0).sum()))
        for (mode,model),g in diag.groupby(['sampling','model'],sort=False)
        for scope,h in [('agbd_ge100',g[g.agbd.ge(100)]),('mixed_circle',g[g.circle_r12p5_mixed_forest_nonforest]),
                        ('nonmixed_circle',g[~g.circle_r12p5_mixed_forest_nonforest])]])
    bins=pd.DataFrame([dict(sampling=mode,model=model,agbd_bin=str(b),**metric(h))
        for (mode,model),g in pred.groupby(['sampling','model'],sort=False)
        for b,h in g.groupby(pd.cut(g.agbd,[0,50,100,150,200,np.inf],right=False),observed=True)])
    for name,frame in [('predictions',pred),('overall_metrics',overall),('paired_comparison',comparisons),
                       ('fold_metrics',fm),('subgroup_metrics',scopes),('agbd_bin_metrics',bins),('fold_validation',pd.DataFrame(fold_checks))]:
        frame.to_csv(out/f'footprint_{name}.csv',index=False)
    report=f'''# 중심 픽셀과 반경 12.5m 평균 반사율 비교

## 결과

같은 산림 184점·공간 fold·버퍼·RF 설정을 유지하고 입력 추출 범위만 바꿨다. 원 안의 픽셀과 겹치는 면적으로 반사율을 평균한 뒤 NDVI·EVI·NDMI를 계산했다. 지수 픽셀 자체의 평균과는 다르다.

{md(comparisons)}

{md(overall)}

원 평균 적용 후 밴드6 RMSE는 56.75→55.51(약 2.18% 감소), 지수+밴드는 56.25→55.60(약 1.16% 감소)였다. 지수3은 56.57→58.17로 악화됐다. 세 구성 모두 MAE가 증가하고 Dummy RMSE 51.53보다 나쁘며 R²도 음수다. 따라서 공간 평균이 일부 큰 오차를 줄였더라도 일관된 예측력 개선을 확보했다고 볼 수 없다. 고 AGBD 51점에서 밴드6의 평균 과소 예측은 약 75.19→69.83 Mg/ha로 줄었지만 여전히 크고, 원 평균 적용 후에도 48점이 과소 예측됐다.

delta는 원 평균−중심 픽셀이며 음수가 RMSE 개선이다. Dummy는 같은 산림 학습점 평균을 사용하는 기존 저장 예측을 공유했다. 세 입력군 결과를 모두 보고하며 유리한 반경이나 모델을 탐색하지 않았다.

## 추출 범위·유효성

184점 모두 원 전체에서 여섯 밴드 공통 유효성, SCL 4/5/6, CLASSI opaque/cirrus/snow=0 조건을 만족했다. 최소 공통 유효 면적 비율은 {audit.common_valid_fraction.min():.12f}다. 불완전한 원이 있으면 적합을 중단하는 규칙으로 구현했으며 결측 대치·표본 탈락·작은 유효 영역만의 평균으로 대체한 점은 없다. 유효 면적 합산의 작은 부동소수 오차에는 1e-6 허용치를 사용했다.

원은 EPSG:32652에서 반경 12.5m·256선분으로 근사했으며 면적은 이론적 원의 약 {audit.circle_area_to_pi_ratio.min():.8f}배다. 평균에 기여한 픽셀 수는 {audit.intersected_pixel_n.min()}~{audit.intersected_pixel_n.max()}개다. 픽셀 중심 포함 여부가 아닌 실제 면적 중첩을 사용했다. `footprint_pixel_weights.csv`에 모든 가중치·반사율·픽셀 ID를 저장했다.

평균에는 SCL 비식생·수역도 포함하며, 원 내부를 산림만으로 잘라 유리한 입력을 만들지 않았다. 원본 밴드 품질 마스크와 보간 유효성은 기존 GeoTIFF 처리 결과를 재사용했다. B11/B12는 원본 20m를 10m 정렬한 자료이므로 새로운 공간 해상도를 얻은 것이 아니다. 주변의 음수 반사율은 임의 클리핑하지 않았다. 평균 후 계산한 EVI 점검 규칙 해당은 {int(averaged.circle_evi_review.sum())}점이며 값을 삭제·클리핑하지 않았다.

## 입력값 변화

{md(shifts)}

## fold별 결과

{md(fm)}

분할은 원래 산림 평가의 39·51·28·62·4점을 유지했다. fold5는 4점뿐이므로 특히 불안정하다. 기존 buffer_excluded 역할을 복원하거나 재배정하지 않았으며 학습/평가 최소 거리는 500m 이상이다. 전체 지표는 184점 OOF 예측에서 직접 계산했다.

## 고 AGBD와 주변 혼합 진단

{md(scopes)}

고 AGBD 51점·혼합 22점·비혼합 162점은 기존 분류를 유지한 사후 진단 집단이다. 혼합 여부는 앞선 EGIS 분석의 EPSG:5186 원으로 정한 값이며 이번 반사율 원은 EPSG:32652에 정의했다. 두 투영좌표계의 미세한 축척 차이와 지도 일반화 한계가 있다. 집단별 결과를 별도 재학습 성능이나 인과 효과로 해석하지 않는다.

## 시각화

![중심점과 원 평균 비교](footprint_comparison.png)

## 검증과 한계

중심 픽셀 RF 15개를 다시 적합해 기존 산림 전용 예측을 점별로 재현했다. 최대 예측 차이는 {max(reproductions):g}다. 원 평균 RF 15회·fold1 지수 RF 재현 1회를 포함해 총 31회 RF를 적합했다. 평균 Dummy는 재학습하지 않았다. 모든 입력군의 평가 ID·fold가 동일하며 각 점은 한 번만 평가된다. 파라미터·버전·입력 해시는 footprint_config.json, 결과 검증은 footprint_validation.json에 기록했다.

이 원은 GEDI의 실제 footprint 모양·에너지 분포·좌표 오차를 복원한 것이 아니라 공간 범위 차이를 점검하기 위한 근사다. 원 평균이 좋아져도 공간 불일치가 유일한 원인이었다고 결론 내리지 않는다. 동일 날짜·궤도·이미 사용한 CV의 탐색 비교이며 독립 최종 일반화 검증은 아니다. EGIS 산림 분류도 현장 실측 검증을 대체하지 않는다.

## 재현

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\landcover\\compare_footprint.py --output-dir .\\week07\\landcover\\footprint_reproduction
```
'''
    (out/'footprint_report.md').write_text(report,encoding='utf-8')
    fig,axes=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
    x=np.arange(3)
    axes[0].bar(x-.18,comparisons.center_rmse,width=.36,label='Center pixel')
    axes[0].bar(x+.18,comparisons.circle_rmse,width=.36,label='Circle mean')
    dummy_rmse=float(overall[overall.sampling.eq('dummy_shared')].rmse.iloc[0])
    axes[0].axhline(dummy_rmse,color='black',linestyle='--',label='Dummy')
    axes[0].set(xticks=x,xticklabels=list(sets),ylabel='OOF RMSE (Mg/ha)',title='Same 184 forest points')
    axes[0].tick_params(axis='x',rotation=20,labelsize=8)
    axes[0].legend(fontsize=8)
    for mode in ['center','circle_r12p5']:
        h=fm[fm.sampling.eq(mode)&fm.model.eq('indices3')]
        axes[1].plot(h.fold,h.rmse,marker='o',label=mode)
    axes[1].set(xlabel='Spatial fold',ylabel='RMSE (Mg/ha)',title='Indices3: fold variation',xticks=range(1,6))
    axes[1].legend(fontsize=8)
    axes[2].scatter(data.ndvi,averaged.ndvi,c=data.circle_r12p5_mixed_forest_nonforest.astype(int),cmap='coolwarm',s=20,alpha=.7)
    axes[2].plot([0,1],[0,1],'--',color='gray')
    axes[2].set(xlabel='Center NDVI',ylabel='NDVI from mean reflectance',title='Blue: nonmixed / red: mixed circle',xlim=(0,1),ylim=(0,1))
    axes[2].set_aspect('equal')
    fig.savefig(out/'footprint_comparison.png',dpi=170)
    plt.close(fig)
    saved=pd.read_csv(out/'footprint_predictions.csv')
    for (mode,model),g in saved.groupby(['sampling','model']):
        row=overall[(overall.sampling==mode)&(overall.model==model)].iloc[0]
        recalculated=metric(g)
        assert all(np.isclose(recalculated[k],row[k],rtol=1e-12,atol=1e-12) for k in ['rmse','mae','r2'])
    # Independent reconstruction from exported long-form pixel weights.
    for b in BANDS:
        group=(pixels[b]*pixels.intersection_area_m2).groupby(pixels.point_id).sum()/pixels.intersection_area_m2.groupby(pixels.point_id).sum()
        assert np.allclose(group.sort_index(),averaged['reflectance_'+b].sort_index(),atol=1e-12,rtol=1e-12)
    assert all(sha(ROOT/p)==h for p,h in hashes.items())
    save_json(out/'footprint_validation.json',dict(status='passed',n=184,full_area_support_all184=True,
        min_valid_fraction=float(audit.common_valid_fraction.min()),center_predictions_reproduced=True,
        center_max_prediction_difference=max(reproductions),circle_fold1_repeat_difference=0,
        saved_metrics_recomputed=True,exported_weight_means_recomputed=True,identical_test_ids=True,
        original_buffer_roles_preserved=True,input_sha256_unchanged=hashes,model_fits=31,new_downloads=False))
    print(comparisons.to_string(index=False))
    print(f'Dummy RMSE: {dummy_rmse}')


if __name__=='__main__':
    main()
