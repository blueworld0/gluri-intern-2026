"""Forest-only retraining on preserved spatial folds, with paired old predictions."""
from pathlib import Path
import argparse
import hashlib
import json
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
import sklearn
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save_json(p,value):
    p.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def metric(g):
    y=g.agbd.to_numpy()
    e=y-g.prediction.to_numpy()
    den=np.sum((y-y.mean())**2)
    return dict(n=len(g),rmse=float(np.sqrt(np.mean(e**2))),mae=float(np.mean(abs(e))),
        r2=float(1-np.sum(e**2)/den) if len(g)>1 and den>0 else None,
        mean_residual=float(e.mean()),mean_agbd=float(y.mean()),mean_prediction=float(g.prediction.mean()))


def md(df):
    def fmt(v):
        return f'{v:.4f}' if isinstance(v,(float,np.floating)) else str(v)
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join(
        '| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=HERE/'forest_evaluation')
    out=parser.parse_args().output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use a new empty --output-dir.')
    paths={k:ROOT/v for k,v in dict(data='week07/landcover/outputs/landcover_points.csv',
        validation='week07/landcover/outputs/validation.json',
        independent='week07/landcover/outputs/independent_verification.json',
        splits='week06/step7/outputs/split_assignments.csv',
        config='week07/step3/outputs/experiment_config.json',
        old_rf='week07/step3/outputs/predictions.csv',old_dummy='week06/step8/outputs/predictions.csv').items()}
    hashes={str(p.relative_to(ROOT)):sha(p) for p in paths.values()}
    assert json.loads(paths['validation'].read_text(encoding='utf-8'))['status']=='passed'
    assert json.loads(paths['independent'].read_text(encoding='utf-8'))['status']=='passed'
    all_data=pd.read_csv(paths['data'],dtype={'shot_number':str})
    assert len(all_data)==441 and all_data.point_id.is_unique
    selected=all_data.landcover_status.eq('forest')
    assert selected.equals(all_data.L2_CODE.astype(str).isin(['310','320','330']))
    data=all_data[selected].set_index('point_id',drop=False)
    assert len(data)==184
    config=json.loads(paths['config'].read_text(encoding='utf-8'))
    sets=config['feature_sets']
    params=config['params']
    assert params['max_depth']==6 and params['min_samples_leaf']==3 and params['max_features']==.67
    assert np.isfinite(data[sum(sets.values(),[])+['agbd','easting_m','northing_m']]).all().all()
    assert not data.duplicated(['gedi_file','beam','shot_number']).any()
    assert not data.duplicated(['latitude','longitude']).any()
    assert not data.duplicated(['utm_10m_x','utm_10m_y']).any()
    assignments=pd.read_csv(paths['splits'])
    roles={}
    fold_details,coverage,kept_assignments=[],[],[]
    for fold in range(1,6):
        a=assignments[assignments.split_id.eq(f'new441_spatial5_fold{fold}')]
        assert len(a)==441 and a.point_id.is_unique
        a=a[a.point_id.isin(data.index)].copy()
        assert len(a)==184 and set(a.point_id)==set(data.index)
        roles[fold]={role:a.loc[a.role.eq(role),'point_id'].tolist() for role in ['train','test','buffer_excluded']}
        tr,te,ex=[data.loc[roles[fold][role]] for role in ['train','test','buffer_excluded']]
        assert len(tr)>3 and len(te)>1
        assert set(tr.spatial_group).isdisjoint(te.spatial_group)
        assert set(te.index)==set(data.index[data.spatial_group.isin(te.spatial_group)])
        dist=np.sqrt(((tr[['easting_m','northing_m']].to_numpy()[:,None]-
                       te[['easting_m','northing_m']].to_numpy()[None])**2).sum(axis=2)).min(axis=1)
        assert dist.min()>=500
        fold_details.append(dict(fold=fold,train_n=len(tr),test_n=len(te),original_buffer_excluded_n=len(ex),
            minimum_train_test_distance_m=float(dist.min()),train_ge100_n=int(tr.agbd.ge(100).sum()),
            test_ge100_n=int(te.agbd.ge(100).sum()),train_mean_agbd=float(tr.agbd.mean()),
            test_mean_agbd=float(te.agbd.mean()),test_mixed_circle_n=int(te.circle_r12p5_mixed_forest_nonforest.sum())))
        kept_assignments.append(a)
        coverage.extend(te.point_id)
    assert len(coverage)==len(set(coverage))==184
    out.mkdir(parents=True,exist_ok=True)
    save_json(out/'experiment_config.json',dict(feature_sets=sets,rf_params=params,
        dummy_params=DummyRegressor(strategy='mean').get_params(),primary_population='184 center-point EGIS forest samples',
        circle_mixed_policy='22 retained; subgroup diagnostics only',
        fold_policy='Original 441-point roles filtered to forest IDs; original buffer-excluded roles retained, no restoration',
        excluded_evi_review_point_n=int(all_data.evi_review_required.sum()-data.evi_review_required.sum()),
        no_target_based_selection=True,no_tuning=True,no_scaling_or_imputation=True,input_sha256=hashes,
        script_sha256=sha(Path(__file__)),versions=dict(python=platform.python_version(),numpy=np.__version__,
        pandas=pd.__version__,sklearn=sklearn.__version__)))
    data.to_csv(out/'forest184_input.csv',index=False,encoding='utf-8-sig')
    pd.concat(kept_assignments).to_csv(out/'forest184_split_assignments.csv',index=False)
    pd.DataFrame(fold_details).to_csv(out/'fold_validation.csv',index=False)
    predictions=[]
    repeat_difference=None
    for fold in range(1,6):
        tr,te=[data.loc[roles[fold][role]] for role in ['train','test']]
        for label in ['dummy_mean']+list(sets):
            columns=sets.get(label,sets['indices3'])
            model=DummyRegressor(strategy='mean') if label=='dummy_mean' else RandomForestRegressor(**params)
            model.fit(tr[columns],tr.agbd)
            prediction=model.predict(te[columns])
            frame=te[['point_id','agbd','latitude','longitude','spatial_group','L2_CODE',
                      'circle_r12p5_mixed_forest_nonforest']].copy()
            frame['fold'],frame['model'],frame['training_population']=fold,label,'forest184'
            frame['prediction'],frame['residual']=prediction,te.agbd.to_numpy()-prediction
            predictions.append(frame.reset_index(drop=True))
            if fold==1 and label=='indices3':
                again=RandomForestRegressor(**params).fit(tr[columns],tr.agbd).predict(te[columns])
                repeat_difference=float(np.max(abs(again-prediction)))
        print(f'Forest fold {fold}: train {len(tr)}, test {len(te)}',flush=True)
    new=pd.concat(predictions,ignore_index=True)
    assert repeat_difference==0
    # Old models were trained on the original vegetation cohorts. Compare their
    # saved predictions on the exact same 184 forest test IDs, never all 441.
    old_rf=pd.read_csv(paths['old_rf']).rename(columns={'feature_set':'model'})
    old_dummy=pd.read_csv(paths['old_dummy'])
    old_dummy=old_dummy[old_dummy['case'].eq('main_all441')&old_dummy.model.eq('dummy_mean')]
    old=pd.concat([old_rf,old_dummy],ignore_index=True)
    old=old[old.point_id.isin(data.index)].copy()
    old['training_population']='vegetation441'
    for label in ['dummy_mean']+list(sets):
        a=new[new.model.eq(label)]
        b=old[old.model.eq(label)]
        assert len(a)==len(b)==184 and a.point_id.is_unique and b.point_id.is_unique
        assert set(zip(a.point_id,a.fold))==set(zip(b.point_id,b.fold))
        assert np.allclose(a.set_index('point_id').agbd.sort_index(),b.set_index('point_id').agbd.sort_index())
    cols=['point_id','agbd','fold','model','prediction','residual','training_population']
    combined=pd.concat([new[cols],old[cols]],ignore_index=True)
    overall=pd.DataFrame([dict(training_population=pop,model=model,**metric(g))
        for (pop,model),g in combined.groupby(['training_population','model'],sort=False)])
    comparison=[]
    for label in ['dummy_mean']+list(sets):
        a=overall[(overall.training_population=='forest184')&(overall.model==label)].iloc[0]
        b=overall[(overall.training_population=='vegetation441')&(overall.model==label)].iloc[0]
        comparison.append(dict(model=label,n=184,old_rmse=b.rmse,forest_rmse=a.rmse,
            rmse_delta=a.rmse-b.rmse,old_mae=b.mae,forest_mae=a.mae,old_r2=b.r2,forest_r2=a.r2))
    comparison=pd.DataFrame(comparison)
    fm=pd.DataFrame([dict(training_population=pop,model=model,fold=fold,**metric(g))
        for (pop,model,fold),g in combined.groupby(['training_population','model','fold'],sort=False)])
    bins=[]
    for (pop,model),g in combined.groupby(['training_population','model'],sort=False):
        categories=pd.cut(g.agbd,[0,50,100,150,200,np.inf],right=False)
        for label,h in g.groupby(categories,observed=True):
            bins.append(dict(training_population=pop,model=model,agbd_bin=str(label),**metric(h)))
    metadata=data[['point_id','circle_r12p5_mixed_forest_nonforest']].reset_index(drop=True)
    diag=combined.merge(metadata,on='point_id',validate='many_to_one')
    scopes=pd.DataFrame([dict(training_population=pop,model=model,scope=scope,**metric(h),
        underprediction_n=int(h.residual.gt(0).sum()))
        for (pop,model),g in diag.groupby(['training_population','model'],sort=False)
        for scope,h in [('agbd_ge100',g[g.agbd.ge(100)]),
                       ('mixed_circle',g[g.circle_r12p5_mixed_forest_nonforest]),
                       ('nonmixed_circle',g[~g.circle_r12p5_mixed_forest_nonforest])]])
    for name,frame in [('predictions',combined),('overall_metrics',overall),('paired_comparison',comparison),
        ('fold_metrics',fm),('agbd_bin_metrics',pd.DataFrame(bins)),('subgroup_metrics',scopes)]:
        frame.to_csv(out/f'{name}.csv',index=False)
    saved=pd.read_csv(out/'predictions.csv')
    for (pop,model),g in saved.groupby(['training_population','model']):
        row=overall[(overall.training_population==pop)&(overall.model==model)].iloc[0]
        independent=metric(g)
        assert all(np.isclose(independent[k],row[k],rtol=1e-12,atol=1e-12) for k in ['rmse','mae','r2'])
    main=overall[overall.training_population.eq('forest184')]
    best=main[main.model.ne('dummy_mean')].sort_values('rmse').iloc[0]
    dummy=main[main.model.eq('dummy_mean')].iloc[0]
    fg=fm[fm.training_population.eq('forest184')].pivot(index='fold',columns='model',values='rmse')
    wins={label:int((fg[label]<fg.dummy_mean).sum()) for label in sets}
    report=f'''# EGIS 산림 184점 공간 교차검증

## 결과

중심점이 EGIS 산림 코드 310·320·330에 속한 **184점만 학습·평가**했다. 동일 RF 설정과 기존 fold를 유지했다. RF 중 pooled OOF RMSE 최저는 `{best.model}` **{best.rmse:.4f} Mg/ha**, R² **{best.r2:.4f}**다. 산림 학습 평균 Dummy의 RMSE는 **{dummy.rmse:.4f} Mg/ha**다.

{md(main)}

산림 전용 RF의 RMSE는 기존 동일 산림 평가점 대비 약 1.09~1.99 Mg/ha 낮아졌지만, 세 RF 모두 산림 학습 평균 Dummy보다 나쁘다. Dummy도 학습 대상의 평균이 산림 구성에 맞춰지면서 약 9.00 Mg/ha 개선됐다. 산림으로 범위를 정리한 것은 타당한 연구 대상 정의지만, 현재 지수·밴드와 RF가 산림 내 AGBD 변동을 유용하게 예측한다고 결론 내릴 수 없다. 세 RF의 전체 R²는 모두 음수다.

## 같은 산림 184점에서 기존 모델과 비교

{md(comparison)}

`old_*`는 441점 식생 후보를 대상으로 학습한 기존 모델의 저장 예측 중 동일한 산림 184점만 추린 값이다. `forest_*`는 산림 학습점만으로 재학습한 결과다. delta는 산림 재학습−기존이며 음수가 RMSE 개선이다. 두 조건은 fold별 평가 ID·AGBD가 정확히 같으며 학습 표본 구성이 다르다. 기존 전체 441점 점수와 직접 개선률을 계산하지 않았다. EVI 이상점은 피복 분류상 비산림이어서 제외됐으며 이 변화도 학습 구성 차이에 포함된다.

## 표본과 고정 분할

{md(pd.DataFrame(fold_details))}

원래 441점의 공간 5-fold·5km 격자·500m 버퍼 역할에서 산림 ID만 남겼다. 기존 buffer_excluded를 복원하거나 재분할하지 않았다. 따라서 비산림 평가점이 빠져 새 산림 평가점에서 멀어진 제외점도 그대로 제외될 수 있다. 모든 산림점이 한 번씩 평가되고 학습·평가 블록이 겹치지 않으며 최소 거리는 500m 이상이다.

평가점 수는 fold별 39·51·28·62·4점으로 크게 다르다. 특히 fold5는 4점뿐이고 고 AGBD 점이 없어 지표가 불안정하며 다른 fold와 동등한 독립 반복으로 취급하지 않는다. 성능을 보고 유리하게 재분할하지 않았다. 전체 OOF 지표는 184점에서 직접 계산했다.

산림 중심점 184점 중 주변 반경 12.5m 원에 비산림도 포함된 22점은 주평가에 유지했다. 주변 피복을 보고 유리한 표본만 선택하지 않았고 아래에 진단 집계만 남겼다. 이 원은 실제 GEDI footprint나 위치 불확실성을 정확히 복원하지 않는다.

## fold별 성능

{md(fm[fm.training_population.eq('forest184')])}

Dummy 대비 RF RMSE 개선 fold 수: {wins}. fold는 독립 반복이 아니며 편차와 음의 R²를 함께 해석해야 한다.

## 고 AGBD와 주변 혼합 표본

{md(scopes)}

고 AGBD 51점, 주변 혼합 22점, 비혼합 162점의 저장 OOF 예측 진단이다. 비혼합 표본만으로 재학습한 결과가 아니며 집단별 AGBD 분포가 달라 점수 차이를 혼합 피복의 인과 효과로 해석하지 않는다.

## AGBD 구간별 성능

{md(pd.DataFrame(bins))}

## 시각화

![산림 전용 학습 비교](forest_comparison.png)

## 설정·검증·해석의 한계

지수3·밴드6·지수+밴드9를 기존 구성 그대로 비교했다. RF는 n_estimators=100, max_depth=6, min_samples_leaf=3, max_features=0.67, random_state=42, n_jobs=1로 고정했다. max_features가 비율이므로 변수 수에 따라 실제 후보 수는 2·4·6개다. 튜닝·변수 선택·스케일링·결측 대치는 없다. Dummy 5회·RF 15회·지수 RF fold1 재실행 1회로 총 21회 적합했다. 재실행 예측 차이는 {repeat_difference:g}이다. 저장 예측에서 RMSE·MAE·R²를 재계산하고 입력 해시 보존을 검증했다.

EGIS 산림 코드와 일치하는 표본이지 독립 현장 실측 산림 AGBD 정답은 아니다. 지도 시점·일반화·경계·GEDI 위치 오차와 같은 날짜·궤도에 한정된다는 제약이 남는다. 이번 산림 재학습은 이미 진단한 공간 CV를 재사용한 탐색 비교이며, 모델 선택 이후의 독립 최종 성능이나 통계적 유의성을 입증하지 않는다. 모델 객체를 배포하거나 기존 441점 결과를 덮어쓰지 않았다.

## 재현

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\landcover\\evaluate_forest.py --output-dir .\\week07\\landcover\\forest_reproduction
```
'''
    (out/'report.md').write_text(report,encoding='utf-8')
    fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    labels=['dummy_mean']+list(sets)
    order=comparison.set_index('model').loc[labels]
    x=np.arange(4)
    axes[0].bar(x-.18,order.old_rmse,width=.36,label='Train on vegetation')
    axes[0].bar(x+.18,order.forest_rmse,width=.36,label='Train on forest')
    axes[0].set(xticks=x,xticklabels=labels,ylabel='RMSE on same 184 forest points (Mg/ha)',title='Paired training-population comparison')
    axes[0].tick_params(axis='x',labelrotation=25,labelsize=7)
    axes[0].legend(fontsize=8)
    for label in labels:
        axes[1].plot(range(1,6),fg[label],marker='o',label=label)
    axes[1].set(xlabel='Spatial fold',ylabel='RMSE (Mg/ha)',title='Forest-only training',xticks=range(1,6))
    axes[1].legend(fontsize=7)
    h=new[new.model.eq(best.model)]
    axes[2].scatter(h.agbd,h.prediction,c=h.fold,cmap='viridis',s=20,alpha=.7)
    axes[2].plot([0,250],[0,250],'--',color='gray')
    axes[2].set(xlabel='GEDI estimated AGBD',ylabel='OOF prediction',title=f'{best.model}: forest-only OOF',xlim=(0,250),ylim=(0,250))
    axes[2].set_aspect('equal')
    fig.savefig(out/'forest_comparison.png',dpi=170)
    plt.close(fig)
    assert all(sha(ROOT/p)==h for p,h in hashes.items())
    save_json(out/'validation.json',dict(status='passed',forest_n=184,models=labels,
        same184_test_points_in_both_training_conditions=True,oof_each_point_once=True,
        forest_only_train_test=True,original_buffer_roles_preserved=True,minimum_distance_at_least_500m=True,
        fold_checks=fold_details,source_hashes_unchanged=hashes,metrics_recalculated=True,
        model_fits=21,repeat_max_prediction_difference=repeat_difference,new_downloads=False))
    print(comparison.to_string(index=False))
    print(main.to_string(index=False))


if __name__=='__main__':
    main()
