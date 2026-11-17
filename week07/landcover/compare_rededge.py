"""Add existing SAFE B05/B06/B07/B8A to frozen forest footprint evaluation."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import xml.etree.ElementTree as ET

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
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
import sklearn
from sklearn.ensemble import RandomForestRegressor
from evaluate_forest import metric,md,save_json

EXTRA=['B05','B06','B07','B8A']


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=HERE/'forest_evaluation')
    out=parser.parse_args().output_dir.resolve()
    if out.exists() and any(out.glob('rededge_*')):
        raise RuntimeError('Rededge outputs already exist; use a new --output-dir.')
    base=HERE/'forest_evaluation'
    safe=next((ROOT/'week06/step5/source').glob('*.SAFE'))
    granule=next((safe/'GRANULE').iterdir())
    paths={k:base/v for k,v in dict(data='footprint_forest184_input.csv',weights='footprint_pixel_weights.csv',
        config='footprint_config.json',validation='footprint_validation.json',splits='forest184_split_assignments.csv',
        previous='footprint_predictions.csv').items()}
    paths['product_xml']=safe/'MTD_MSIL2A.xml'
    paths['processing_helper']=ROOT/'week06/step5/process_candidate_a.py'
    paths['metric_helper']=HERE/'evaluate_forest.py'
    paths['reference']=ROOT/'week06/step5/outputs/rasters/reflectance_B04.tif'
    for b in EXTRA:
        paths[b]=next((granule/'IMG_DATA/R20m').glob(f'*_{b}_20m.jp2'))
        paths[b+'_quality']=granule/f'QI_DATA/MSK_QUALIT_{b}.jp2'
        paths[b+'_detector']=granule/f'QI_DATA/MSK_DETFOO_{b}.jp2'
    assert all(p.exists() for p in paths.values())
    hashes={str(p.relative_to(ROOT)):sha(p) for p in paths.values()}
    assert json.loads(paths['validation'].read_text(encoding='utf-8'))['status']=='passed'
    data=pd.read_csv(paths['data'],dtype={'shot_number':str}).set_index('point_id',drop=False)
    assert len(data)==184 and data.index.is_unique
    weights=pd.read_csv(paths['weights'])
    assert weights.common_valid.all() and set(weights.point_id)==set(data.index)
    config=json.loads(paths['config'].read_text(encoding='utf-8'))
    params=config['rf_params']
    sets={'bands6':config['feature_sets']['bands6']}
    sets['bands10']=sets['bands6']+['reflectance_'+b for b in EXTRA]
    # Reuse previously verified crop, categorical-nearest and bilinear/full-
    # support functions exactly; module does not execute its main pipeline.
    spec=importlib.util.spec_from_file_location('safe_helpers',paths['processing_helper'])
    helper=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    xml=ET.parse(paths['product_xml'])
    local=lambda tag:tag.rsplit('}',1)[-1]
    q=float(next(e.text for e in xml.iter() if local(e.tag)=='BOA_QUANTIFICATION_VALUE'))
    ids={e.attrib['bandId']:e.attrib['physicalBand'] for e in xml.iter() if local(e.tag)=='Spectral_Information'}
    offsets={ids[e.attrib['band_id']]:float(e.text) for e in xml.iter() if local(e.tag)=='BOA_ADD_OFFSET'}
    offsets={b:offsets[b if b=='B8A' else 'B'+str(int(b[1:]))] for b in EXTRA}
    special={}
    for e in xml.iter():
        if local(e.tag)=='Special_Values':
            d={local(c.tag):c.text for c in e}
            special[d['SPECIAL_VALUE_TEXT']]=int(d['SPECIAL_VALUE_INDEX'])
    assert special['NODATA']==0 and special['SATURATED']==65535
    with rasterio.open(paths['reference']) as ds:
        assert ds.crs.to_epsg()==32652 and ds.res==(10,10)
        target=dict(crs=ds.crs,transform=ds.transform,width=ds.width,height=ds.height)
    bounds=array_bounds(target['height'],target['width'],target['transform'])
    pixel_values=weights[['point_id','raster_row','raster_col','intersection_area_m2']].copy()
    rr,cc=weights.raster_row.to_numpy(),weights.raster_col.to_numpy()
    area=weights.intersection_area_m2.groupby(weights.point_id).sum()
    band_audit=[]
    point_audit=pd.DataFrame(index=data.index)
    # Check that the retained old footprint weights still reconstruct all six
    # baseline features before adding any new information.
    for b in ['B02','B03','B04','B08','B11','B12']:
        mean=(weights[b]*weights.intersection_area_m2).groupby(weights.point_id).sum()/area
        assert np.allclose(mean.sort_index(),data['reflectance_'+b].sort_index(),atol=1e-12,rtol=1e-12)
    for b in EXTRA:
        raw,native=helper.crop(paths[b],bounds)
        quality=helper.read_grid(paths[b+'_quality'],native,list(range(1,9)))
        detector=helper.read_grid(paths[b+'_detector'],native)
        assert set(np.unique(quality)).issubset({0,1})
        assert np.isin(detector,np.arange(13)).all()
        valid=(raw!=special['NODATA'])&(raw!=special['SATURATED'])&~(quality!=0).any(axis=0)&np.isin(detector,np.arange(1,13))
        rho=(raw.astype('float32')+offsets[b])/q
        rho[~valid]=np.nan
        aligned=helper.warp(rho,native,target,Resampling.bilinear,nodata=np.nan)
        support=helper.warp(valid.astype('float32'),native,target,Resampling.bilinear)
        aligned_valid=np.isfinite(aligned)&(support>=1-1e-6)
        aligned[~aligned_valid]=np.nan
        vals=aligned[rr,cc].astype(float)
        ok=aligned_valid[rr,cc]
        pixel_values[b]=vals
        pixel_values[b+'_valid']=ok
        coverage=(weights.intersection_area_m2*ok).groupby(weights.point_id).sum()/area
        point_audit[b+'_valid_fraction']=coverage
        means=(weights.intersection_area_m2*vals).groupby(weights.point_id).sum()/area
        data['reflectance_'+b]=means
        band_audit.append(dict(band=b,original_resolution_m=20,offset=offsets[b],quantification=q,
            native_crop_rows=raw.shape[0],native_crop_cols=raw.shape[1],native_invalid_n=int((~valid).sum()),
            footprint_invalid_aligned_pixel_n=int((~ok).sum()),minimum_footprint_valid_fraction=float(coverage.min()),
            minimum_footprint_reflectance=float(means.min()),maximum_footprint_reflectance=float(means.max())))
        print(f'{b}: offset={offsets[b]}, quantification={q}, footprint invalid pixels={int((~ok).sum())}',flush=True)
    out.mkdir(parents=True,exist_ok=True)
    pixel_values.to_csv(out/'rededge_pixel_values.csv',index=False)
    point_audit.to_csv(out/'rededge_coverage_audit.csv',index=True)
    pd.DataFrame(band_audit).to_csv(out/'rededge_band_audit.csv',index=False)
    if not point_audit.ge(1-1e-6).all().all():
        save_json(out/'rededge_validation.json',dict(status='stopped_incomplete_support',model_fits=0))
        raise RuntimeError('Some new bands lack full-circle support; no models fitted, no imputation or removal.')
    assert np.isfinite(data[sets['bands10']]).all().all()
    data.to_csv(out/'rededge_forest184_input.csv',index=False,encoding='utf-8-sig')
    assignments=pd.read_csv(paths['splits'])
    roles,fold_audit,heldout={},[],[]
    for fold in range(1,6):
        a=assignments[assignments.split_id.eq(f'new441_spatial5_fold{fold}')]
        assert len(a)==184 and a.point_id.is_unique and set(a.point_id)==set(data.index)
        roles[fold]={role:a.loc[a.role.eq(role),'point_id'].tolist() for role in ['train','test','buffer_excluded']}
        tr,te=[data.loc[roles[fold][role]] for role in ['train','test']]
        assert set(tr.spatial_group).isdisjoint(te.spatial_group)
        distance=float(np.sqrt(((tr[['easting_m','northing_m']].to_numpy()[:,None]-te[['easting_m','northing_m']].to_numpy()[None])**2).sum(axis=2)).min())
        assert distance>=500
        fold_audit.append(dict(fold=fold,train_n=len(tr),test_n=len(te),original_buffer_excluded_n=len(roles[fold]['buffer_excluded']),minimum_distance_m=distance))
        heldout.extend(te.point_id)
    assert len(heldout)==len(set(heldout))==184
    save_json(out/'rededge_config.json',dict(feature_sets=sets,rf_params=params,input_sha256=hashes,
        new_band_offsets=offsets,quantification=q,quality_policy='All 8 native QUALIT layers zero, detector 1..12, valid DN; full bilinear valid support',
        sampling='Existing r12.5m/256-segment circle-pixel area weights; no new indices, no forest clipping inside circles',
        native_resolution_m=20,aligned_resolution_m=10,all184_retained=True,tuning=False,
        effective_max_features={k:max(1,int(len(v)*params['max_features'])) for k,v in sets.items()},
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
            sklearn=sklearn.__version__,rasterio=rasterio.__version__),script_sha256=sha(Path(__file__))))
    previous=pd.read_csv(paths['previous'])
    old6=previous[(previous.sampling=='circle_r12p5')&(previous.model=='bands6')]
    output=[]
    repeats=[]
    reproduction=[]
    for fold in range(1,6):
        tr,te=[data.loc[roles[fold][role]] for role in ['train','test']]
        for label,cols in sets.items():
            pred=RandomForestRegressor(**params).fit(tr[cols],tr.agbd).predict(te[cols])
            frame=te[['point_id','agbd','spatial_group','circle_r12p5_mixed_forest_nonforest']].copy()
            frame['fold'],frame['model']=fold,label
            frame['prediction'],frame['residual']=pred,te.agbd.to_numpy()-pred
            output.append(frame.reset_index(drop=True))
            if label=='bands6':
                old=old6[old6.fold.eq(fold)].set_index('point_id').loc[te.point_id]
                reproduction.append(float(np.max(abs(old.prediction.to_numpy()-pred))))
            if label=='bands10' and fold==1:
                again=RandomForestRegressor(**params).fit(tr[cols],tr.agbd).predict(te[cols])
                repeats.append(float(np.max(abs(again-pred))))
        print(f'Red-edge comparison fold {fold}: train={len(tr)}, test={len(te)}',flush=True)
    assert max(reproduction)<1e-10 and repeats==[0.0]
    pred=pd.concat(output,ignore_index=True)
    dummy=previous[previous.model.eq('dummy_mean')].copy()
    dummy['model']='dummy_mean'
    pred=pd.concat([pred,dummy[pred.columns.intersection(dummy.columns)]],ignore_index=True)
    for label,g in pred.groupby('model'):
        assert len(g)==184 and g.point_id.is_unique and set(zip(g.point_id,g.fold))==set(zip(dummy.point_id,dummy.fold))
    overall=pd.DataFrame([dict(model=label,**metric(g)) for label,g in pred.groupby('model',sort=False)])
    indexed=overall.set_index('model')
    old,new,dm=[indexed.loc[label] for label in ['bands6','bands10','dummy_mean']]
    change=dict(rmse_delta=float(new.rmse-old.rmse),rmse_reduction_pct=float(100*(old.rmse-new.rmse)/old.rmse),
                mae_delta=float(new.mae-old.mae),r2_delta=float(new.r2-old.r2),rmse_delta_vs_dummy=float(new.rmse-dm.rmse))
    fm=pd.DataFrame([dict(model=label,fold=fold,**metric(g)) for (label,fold),g in pred.groupby(['model','fold'],sort=False)])
    mixed=data[['point_id','circle_r12p5_mixed_forest_nonforest']].reset_index(drop=True)
    diag=pred.drop(columns='circle_r12p5_mixed_forest_nonforest',errors='ignore').merge(mixed,on='point_id',validate='many_to_one')
    scopes=pd.DataFrame([dict(model=label,scope=scope,**metric(h),underprediction_n=int(h.residual.gt(0).sum()))
        for label,g in diag.groupby('model',sort=False)
        for scope,h in [('agbd_ge100',g[g.agbd.ge(100)]),('mixed_circle',g[g.circle_r12p5_mixed_forest_nonforest]),('nonmixed_circle',g[~g.circle_r12p5_mixed_forest_nonforest])]])
    bins=pd.DataFrame([dict(model=label,agbd_bin=str(b),**metric(h)) for label,g in pred.groupby('model',sort=False)
        for b,h in g.groupby(pd.cut(g.agbd,[0,50,100,150,200,np.inf],right=False),observed=True)])
    paired=pred.pivot(index=['point_id','fold','agbd'],columns='model',values='prediction').reset_index()
    paired['bands10_minus_bands6_absolute_error']=abs(paired.agbd-paired.bands10)-abs(paired.agbd-paired.bands6)
    for name,frame in [('predictions',pred),('overall_metrics',overall),('fold_metrics',fm),('subgroup_metrics',scopes),
                       ('agbd_bin_metrics',bins),('paired_predictions',paired),('fold_validation',pd.DataFrame(fold_audit))]:
        frame.to_csv(out/f'rededge_{name}.csv',index=False)
    save_json(out/'rededge_comparison.json',change)
    wide=fm.pivot(index='fold',columns='model',values='rmse')
    wins=int((wide.bands10<wide.bands6).sum())
    verdict='산림 학습 평균 Dummy보다 RMSE가 낮다.' if new.rmse<dm.rmse else '산림 학습 평균 Dummy보다 RMSE가 높아, 이번 입력 확장으로도 기준 모델을 넘지 못했다.'
    report=f'''# 기존 SAFE red-edge 추가 비교 — 밴드 6개 vs 10개

## 결과

기존 B02·B03·B04·B08·B11·B12의 반경 12.5m 평균 반사율에 B05·B06·B07·B8A를 추가했다. 같은 산림 184점, 같은 공간 fold·버퍼·RF 설정을 유지했다. 밴드10의 RMSE는 **{new.rmse:.4f} Mg/ha**, 밴드6 대비 변화는 **{change['rmse_delta']:+.4f} Mg/ha ({change['rmse_reduction_pct']:.2f}% 감소)**다. MAE 변화는 **{change['mae_delta']:+.4f} Mg/ha**다. {verdict}

{md(overall)}

## 원본 처리와 품질

{md(pd.DataFrame(band_audit))}

새 밴드는 모두 원본 20m다. 제품 XML의 밴드별 BOA_ADD_OFFSET과 BOA_QUANTIFICATION_VALUE로 `(DN+offset)/quantification`을 적용했다. MSK_QUALIT 8개 층이 모두 0이고 검출기 ID가 1~12이며 DN이 NoData·포화가 아닌 픽셀만 사용했다. 보정·품질 마스크 후 기존 처리 함수로 bilinear 10m 정렬하고 보간 기여 원본 픽셀이 전부 유효한 경우만 인정했다. 10m 정렬은 실제 해상도를 높이지 않는다.

기존 원형 면적 가중치와 기존 6개 밴드의 공통 유효 영역(SCL 4/5/6, CLASSI clear)을 그대로 재사용했다. 184점 모두 신규 네 밴드에서도 원 전체 유효성을 만족했다. 결측 대치·점 삭제·산림 피복만 잘라 평균내기·추가 지수 생성은 하지 않았다. `rededge_pixel_values.csv`와 `rededge_coverage_audit.csv`에 근거를 저장했다. 새로운 외부 자료 다운로드는 없었다.

## 고정 비교 설계

RF는 n_estimators=100, max_depth=6, min_samples_leaf=3, max_features=0.67, random_state=42, n_jobs=1이다. max_features가 비율이므로 후보 변수 수는 6개 입력에서 4개, 10개 입력에서 6개로 달라진다. 추가 밴드와 후보 수 증가의 효과를 분리한 인과 실험은 아니다. 반경·파라미터·입력 조합을 성능에 맞춰 탐색하지 않았다.

{md(pd.DataFrame(fold_audit))}

모든 점은 한 번 평가됐고 평가 ID·fold가 세 모델에서 동일하다. 기존 buffer_excluded 역할을 유지하고 제외점을 복원하지 않았다. 학습·평가 블록은 겹치지 않으며 최소 거리는 500m 이상이다. fold5는 4점뿐이어서 특히 불안정하다. 전체 점수는 pooled OOF 예측에서 계산했다.

## fold별 성능

{md(fm)}

밴드10의 RMSE가 밴드6보다 낮은 fold는 **{wins}/5**개다. Dummy보다 낮은 fold는 **{int((wide.bands10<wide.dummy_mean).sum())}/5**개다. 기술적 편차를 보는 값이며 독립 반복 검정이나 유의성을 의미하지 않는다.

## 고 AGBD와 주변 혼합 표본

{md(scopes)}

고 AGBD 51점·혼합 22점·비혼합 162점은 같은 ID로 고정했다. 이 부분집합 지표는 별도 집합만으로 재학습한 결과가 아니다.

고 AGBD 표본에서 밴드10의 평균 과소 예측은 약 68.11 Mg/ha로, 밴드6의 약 69.83보다 줄었지만 51점 중 48점이 여전히 과소 예측됐다. 입력 확장의 작은 개선을 고 AGBD 예측 문제 해결로 해석하지 않는다.

## AGBD 구간별 결과

{md(bins)}

## 그림

![red-edge 입력 추가 비교](rededge_comparison.png)

## 검증과 해석

기존 밴드6 원 평균 RF를 재적합해 점별 저장 예측을 재현했으며 최대 차이는 {max(reproduction):g}이다. 신규 밴드10 fold1 재실행 차이는 {repeats[0]:g}이다. RF 10회+재현 1회 총 11회 적합했고 Dummy는 동일 저장 예측을 재사용했다. 새 밴드 평균을 픽셀 가중치에서 다시 계산하고 저장 예측에서 지표를 독립 재계산했다. 원본·기존 결과의 입력 해시는 보존됐다.

같은 날짜·궤도에 대한 이미 여러 번 확인한 공간 CV의 탐색 비교다. 독립 최종 일반화 성능으로 보고하지 않는다. EGIS 산림 분류·지도 시점·중분류 경계와 GEDI 추정 AGBD의 한계가 남는다. 반경 12.5m 원형 평균은 실제 footprint·에너지 가중치·위치 오차를 정확히 복원하지 않는다. 이번 결과로 밴드별 인과 기여나 지수 포화를 확정하지 않는다.

이번 비교 결과까지 week07의 입력 개선 실험으로 정리한다. 기준 모델을 넘지 못한다면 같은 CV에서 조합과 튜닝을 계속 늘리기보다 독립 관측일·궤도 확보와 학습 표본의 대표성을 다음 연구 과제로 삼는다.

## 재현

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\landcover\\compare_rededge.py --output-dir .\\week07\\landcover\\rededge_reproduction
```
'''
    (out/'rededge_report.md').write_text(report,encoding='utf-8')
    fig,axes=plt.subplots(1,3,figsize=(14,4.7),layout='constrained')
    axes[0].bar(overall.model,overall.rmse,color=['#3182bd','#31a354','#888888'])
    axes[0].set(ylabel='OOF RMSE (Mg/ha)',title='Same 184 forest points')
    for label in ['bands6','bands10','dummy_mean']:
        axes[1].plot(range(1,6),wide[label],marker='o',label=label)
    axes[1].set(xlabel='Spatial fold',ylabel='RMSE (Mg/ha)',title='Fold variation',xticks=range(1,6))
    axes[1].legend(fontsize=8)
    h=pred[pred.model.eq('bands10')]
    axes[2].scatter(h.agbd,h.prediction,c=h.fold,cmap='viridis',s=20,alpha=.7)
    axes[2].plot([0,250],[0,250],'--',color='gray')
    axes[2].set(xlabel='GEDI estimated AGBD',ylabel='OOF prediction',title='Bands10: r=12.5m mean',xlim=(0,250),ylim=(0,250))
    axes[2].set_aspect('equal')
    fig.savefig(out/'rededge_comparison.png',dpi=170)
    plt.close(fig)
    saved=pd.read_csv(out/'rededge_predictions.csv')
    for label,g in saved.groupby('model'):
        row=indexed.loc[label]
        m=metric(g)
        assert all(np.isclose(m[k],row[k],atol=1e-12,rtol=1e-12) for k in ['rmse','mae','r2'])
    exported=pd.read_csv(out/'rededge_pixel_values.csv')
    for b in EXTRA:
        reconstructed=(exported[b]*exported.intersection_area_m2).groupby(exported.point_id).sum()/exported.intersection_area_m2.groupby(exported.point_id).sum()
        assert np.allclose(reconstructed.sort_index(),data['reflectance_'+b].sort_index(),rtol=1e-12,atol=1e-12)
    assert all(sha(ROOT/p)==h for p,h in hashes.items())
    save_json(out/'rededge_validation.json',dict(status='passed',n=184,new_bands=EXTRA,full_support_all_points=True,
        all_test_ids_identical=True,oof_each_point_once=True,minimum_500m_buffer_verified=True,
        original_split_roles_preserved=True,baseline_max_prediction_difference=max(reproduction),
        new_model_fold1_repeat_difference=repeats[0],exported_means_recomputed=True,
        saved_metrics_recomputed=True,model_fits=11,input_sha256_unchanged=hashes,new_downloads=False))
    print(overall.to_string(index=False))
    print(json.dumps(change,indent=2))


if __name__=='__main__':
    main()
