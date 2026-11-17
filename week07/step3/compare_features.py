"""Frozen spatial-CV feature comparison using existing local reflectance only."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform

ROOT = Path(__file__).resolve().parents[2]
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[key] = '1'
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestRegressor

INDICES = ['ndvi', 'evi', 'ndmi']
BANDS = ['reflectance_' + b for b in ['B02', 'B03', 'B04', 'B08', 'B11', 'B12']]
SETS = {'indices3': INDICES, 'bands6': BANDS, 'indices3_bands6': INDICES + BANDS}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def metrics(g):
    y, p = g.agbd.to_numpy(), g.prediction.to_numpy()
    e = y-p
    return dict(n=len(g), rmse=float(np.sqrt(np.mean(e**2))), mae=float(np.mean(abs(e))),
                r2=float(1-np.sum(e**2)/np.sum((y-y.mean())**2)), mean_residual=float(e.mean()))


def md(df):
    def fmt(v):
        return f'{v:.4f}' if isinstance(v, (float, np.floating)) else str(v)
    return '| ' + ' | '.join(df.columns) + ' |\n| ' + ' | '.join(['---']*len(df.columns)) + ' |\n' + '\n'.join(
        '| ' + ' | '.join(fmt(v) for v in row) + ' |' for row in df.itertuples(index=False, name=None))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).parent/'outputs')
    out = parser.parse_args().output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use a new empty --output-dir.')
    paths = {k: ROOT/v for k, v in dict(
        data='week06/step7/outputs/new441_evaluation_input.csv',
        bands='week06/step5/outputs/vegetation_candidates.csv',
        splits='week06/step7/outputs/split_assignments.csv',
        config='week06/step8/outputs/experiment_config.json',
        previous='week06/step8/outputs/predictions.csv',
        audit='week07/step2/outputs/validation.json',
        errors='week07/step1/outputs/top20_absolute_errors.csv').items()}
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths.values()}
    assert json.loads(paths['audit'].read_text(encoding='utf-8'))['status'] == 'passed'
    data = pd.read_csv(paths['data'], dtype={'shot_number': str})
    bands = pd.read_csv(paths['bands'], dtype={'shot_number': str})
    assert len(data) == len(bands) == 441 and data.point_id.is_unique and bands.point_id.is_unique
    assert set(data.point_id) == set(bands.point_id)
    bands = bands.set_index('point_id').loc[data.point_id]
    for a,b in [('agbd','agbd_Mg_ha'), ('ndvi','ndvi'), ('evi','evi'), ('ndmi','ndmi')]:
        assert np.array_equal(data[a].to_numpy(), bands[b].to_numpy())
    for b in BANDS:
        data[b] = bands[b].to_numpy()
    data = data.set_index('point_id', drop=False)
    assert np.isfinite(data[INDICES+BANDS+['agbd','easting_m','northing_m']]).all().all()
    assert not data.duplicated(['latitude','longitude']).any()
    params = json.loads(paths['config'].read_text(encoding='utf-8'))['models']['week05_spatial_B_fixed']
    assert params['max_depth'] == 6 and params['min_samples_leaf'] == 3 and params['max_features'] == .67
    assignments = pd.read_csv(paths['splits'])
    previous = pd.read_csv(paths['previous'])
    previous = previous[previous['case'].eq('main_all441')]
    roles, fold_checks, coverage = {}, [], []
    for fold in range(1,6):
        a = assignments[assignments.split_id.eq(f'new441_spatial5_fold{fold}')]
        assert len(a) == 441 and a.point_id.is_unique and set(a.point_id) == set(data.index)
        assert set(a.role) <= {'train','test','buffer_excluded'}
        roles[fold] = {r: a.loc[a.role.eq(r),'point_id'].tolist() for r in ['train','test','buffer_excluded']}
        tr, te, ex = [data.loc[roles[fold][r]] for r in ['train','test','buffer_excluded']]
        assert set(tr.spatial_group).isdisjoint(te.spatial_group)
        assert set(te.index) == set(data.index[data.spatial_group.isin(te.spatial_group)])
        xy = te[['easting_m','northing_m']].to_numpy()
        distance = lambda g: np.sqrt(((g[['easting_m','northing_m']].to_numpy()[:,None]-xy[None])**2).sum(axis=2)).min(axis=1)
        minimum = float(distance(tr).min())
        assert minimum >= 500 and (len(ex)==0 or (distance(ex)<500).all())
        fold_checks.append(dict(fold=fold, train_n=len(tr), test_n=len(te), excluded_n=len(ex), min_distance_m=minimum))
        coverage.extend(te.point_id)
    assert len(coverage) == len(set(coverage)) == 441
    out.mkdir(parents=True, exist_ok=True)
    config = dict(feature_sets=SETS, target='agbd', params=params, cases=['main_all441'],
                  note='All 441 points retained, including preflagged EVI point; no tuning or new data.',
                  effective_max_features={k:max(1,int(.67*len(v))) for k,v in SETS.items()},
                  input_sha256=hashes, versions=dict(python=platform.python_version(), numpy=np.__version__,
                  pandas=pd.__version__, sklearn=sklearn.__version__), script_sha256=sha(Path(__file__)))
    save_json(out/'experiment_config.json', config)
    data.to_csv(out/'evaluation_input.csv', index=False)
    predictions, fold_metrics, repeat = [], [], []
    for fold in range(1,6):
        tr, te = [data.loc[roles[fold][r]] for r in ['train','test']]
        for label, features in SETS.items():
            model = RandomForestRegressor(**params).fit(tr[features], tr.agbd)
            p = model.predict(te[features])
            frame = te[['point_id','agbd','latitude','longitude','spatial_group']].copy()
            frame['fold'], frame['feature_set'] = fold, label
            frame['prediction'], frame['residual'] = p, te.agbd.to_numpy()-p
            predictions.append(frame.reset_index(drop=True))
            fold_metrics.append(dict(fold=fold, feature_set=label, train_n=len(tr), **metrics(frame)))
            if fold == 1 and label == 'indices3_bands6':
                again = RandomForestRegressor(**params).fit(tr[features], tr.agbd).predict(te[features])
                repeat.append(float(np.max(abs(p-again))))
        print(f'Fold {fold}: three feature sets evaluated', flush=True)
    pred = pd.concat(predictions, ignore_index=True)
    reference = pred[pred.feature_set.eq('indices3')]
    for label in SETS:
        g = pred[pred.feature_set.eq(label)]
        assert len(g)==441 and g.point_id.is_unique
        assert set(zip(g.point_id,g.fold)) == set(zip(reference.point_id,reference.fold))
    # Baseline rerun must reproduce every original OOF prediction, not just metrics.
    baseline = pred[pred.feature_set.eq('indices3')].merge(
        previous[previous.model.eq('week05_spatial_B_fixed')][['point_id','fold','prediction']],
        on=['point_id','fold'], suffixes=('','_week06'), validate='one_to_one')
    baseline_difference = float(abs(baseline.prediction-baseline.prediction_week06).max())
    assert len(baseline)==441 and baseline_difference < 1e-10 and repeat == [0.0]
    overall = pd.DataFrame([dict(feature_set=k, **metrics(pred[pred.feature_set.eq(k)])) for k in SETS])
    dummy = previous[previous.model.eq('dummy_mean')]
    dummy_metrics = metrics(dummy)
    baseline_rmse = float(overall.iloc[0].rmse)
    overall['rmse_delta_vs_indices3'] = overall.rmse-baseline_rmse
    overall['rmse_reduction_pct_vs_indices3'] = 100*(baseline_rmse-overall.rmse)/baseline_rmse
    overall['rmse_delta_vs_dummy'] = overall.rmse-dummy_metrics['rmse']
    fm = pd.DataFrame(fold_metrics)
    pred['agbd_bin'] = pd.cut(pred.agbd, [0,50,100,150,200,np.inf], right=False,
                             labels=['[0,50)','[50,100)','[100,150)','[150,200)','[200,inf)']).astype(str)
    bins = pd.DataFrame([dict(feature_set=k, agbd_bin=b, **metrics(g), mean_prediction=float(g.prediction.mean()))
        for (k,b),g in pred.groupby(['feature_set','agbd_bin'], sort=False)])
    top_ids = set(pd.read_csv(paths['errors']).point_id)
    high = pd.DataFrame([dict(feature_set=k, scope=scope, **metrics(g), underprediction_n=int((g.residual>0).sum()))
        for k in SETS for scope,mask in [('agbd_ge100',pred.agbd.ge(100)), ('previous_top20',pred.point_id.isin(top_ids))]
        for g in [pred[pred.feature_set.eq(k)&mask]]])
    paired = pred.pivot(index=['point_id','fold','agbd','spatial_group'], columns='feature_set', values='prediction').reset_index()
    for label in ['bands6','indices3_bands6']:
        paired[label+'_absolute_error_delta'] = abs(paired.agbd-paired[label])-abs(paired.agbd-paired.indices3)
    for name,frame in [('predictions',pred),('overall_metrics',overall),('fold_metrics',fm),
                       ('agbd_bin_metrics',bins),('high_agbd_and_previous_top20',high),
                       ('paired_predictions',paired),('fold_validation',pd.DataFrame(fold_checks))]:
        frame.to_csv(out/f'{name}.csv', index=False)
    saved = pd.read_csv(out/'predictions.csv')
    for label in SETS:
        recalculated = metrics(saved[saved.feature_set.eq(label)])
        row = overall[overall.feature_set.eq(label)].iloc[0]
        assert all(np.isclose(recalculated[k],row[k],rtol=1e-12,atol=1e-12) for k in ['rmse','mae','r2'])
    assert all(sha(ROOT/p)==h for p,h in hashes.items())
    wide = fm.pivot(index='fold',columns='feature_set',values='rmse')
    best = overall.loc[overall.rmse.idxmin()]
    report = f'''# Week07 Step 3 — 기존 밴드 추가 효과 비교

## 결과

새 자료 없이 기존 Sentinel-2 SAFE에서 추출·검증한 반사율 6개만 사용했다. 동일한 441점, 공간 5-fold·5km 격자·500m 버퍼, 동일 RF 설정을 유지하고 입력변수만 비교했다. 최저 pooled OOF RMSE는 `{best.feature_set}`의 **{best.rmse:.4f} Mg/ha**다. 이는 이번 고정 비교 내 관측 결과이며 외부 일반화 검증은 아니다.

{md(overall)}

밴드만 사용하면 RMSE는 기존보다 약 0.30 Mg/ha(0.68%) 낮지만, 지수에 밴드를 추가하면 약 0.19 Mg/ha 높다. 입력 추가만으로 뚜렷하고 일관된 개선을 확보하지 못했다. bands6의 fold5 RMSE는 약 5.41 낮아진 반면 fold1은 약 2.42 높아져 지역별 결과가 다르다. 고 AGBD 56점은 세 입력군 모두 과소 예측하며 평균 과소 예측은 약 90.74→91.22(밴드만) 또는 92.22(지수+밴드) Mg/ha로 줄지 않았다. 따라서 이 결과만으로 기존 모델을 대체할 근거는 부족하다.

Dummy는 week06의 동일 평가점별 저장 예측을 재사용했다: RMSE {dummy_metrics['rmse']:.4f}, MAE {dummy_metrics['mae']:.4f}, R² {dummy_metrics['r2']:.4f}. delta는 신규−기존이며 음수가 RMSE 개선이다. pooled 지표는 전체 OOF 예측에서 직접 계산했으며 fold 평균과 다르다.

## 입력과 통제

- `indices3`: NDVI·EVI·NDMI.
- `bands6`: B02·B03·B04·B08·B11·B12 반사율.
- `indices3_bands6`: 위 지수 3개 + 반사율 6개.
- RF: n_estimators=100, max_depth=6, min_samples_leaf=3, max_features=0.67, random_state=42, n_jobs=1. 전체 파라미터와 라이브러리 버전은 experiment_config.json에 기록했다.
- 전역 스케일링·결측 대치·변수 선택·하이퍼파라미터 탐색은 없다. point_id 기준으로 반사율을 결합하고 기존 지수·AGBD 값 일치를 검증했다.
- 기존 EVI 점검 대상 1점도 모든 입력군의 동일 441점에 유지했다. 이번에는 추가 민감도 실험으로 유리한 조건을 선택하지 않았다.
- B11/B12는 원본 20m에서 기존 방식으로 10m 정렬한 반사율이며 실제 해상도가 높아진 것은 아니다.

`max_features=0.67`은 비율이므로 분할 후보 변수 수는 입력군에 따라 2·4·6개로 달라진다. 같은 설정을 사용한 입력 구성 비교이며, 변수 추가의 순수 효과와 후보 수 증가의 효과를 분리한 인과 실험은 아니다.

## 공간 fold별 결과

{md(fm)}

기존 지수 대비 RMSE가 낮은 fold 수: bands6 **{int((wide.bands6<wide.indices3).sum())}/5**, 지수+밴드 **{int((wide.indices3_bands6<wide.indices3).sum())}/5**. fold는 독립 반복 실험이 아니며 유의성·신뢰구간을 계산하지 않았다.

## 고 AGBD 및 기존 상위 오차점

{md(high)}

`agbd_ge100`은 같은 고 AGBD 56점, `previous_top20`은 step1에서 고정한 동일 20점이다. mean_residual 양수는 과소 예측이다. 기존 오차가 큰 점을 사후 선택한 진단 집합이므로 이 부분집합의 개선을 독립 검증 성능으로 해석하지 않는다.

## AGBD 구간별 결과

{md(bins)}

## 시각화

![입력변수별 성능 비교](feature_comparison.png)

## 검증과 한계

기존 지수 모델 441점의 week06 예측 재현 최대 차이 {baseline_difference:g}. 지수+밴드 fold1 재실행 차이 {repeat[0]:g}. 모든 입력군에서 평가 ID·fold가 동일하고 각 점은 한 번만 평가되며, 학습/평가 블록 불일치와 최소 500m 거리를 확인했다. 저장 예측으로 지표를 독립 재계산했고 입력 해시는 실행 전후 동일하다. 총 16회 RF 적합이며 Dummy는 재학습하지 않았다.

상위 오차 진단에 사용한 기존 공간 CV를 다시 비교한 탐색 실험이다. 새 독립 관측일·궤도·현장 표본을 추가하지 않았으므로 모델 선택 후 최종 일반화 성능으로 보고하지 않는다. 산림 미검증 식생 후보와 GEDI 추정 AGBD라는 기존 한계가 유지된다. 이번 밴드 비교에는 새 자료가 필요 없지만, 산림/농경지 구분에는 별도 피복 자료가, 시간·궤도 일반화 확인에는 독립 관측이 필요하다.

## 재현

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\step3\\compare_features.py --output-dir .\\week07\\step3\\reproduction_outputs
```
'''
    (out/'report.md').write_text(report,encoding='utf-8')
    fig, axes = plt.subplots(1,3, figsize=(15,4.5), layout='constrained')
    labels = list(SETS)
    axes[0].bar(labels,overall.rmse,color=['#737373','#3182bd','#31a354'])
    axes[0].axhline(dummy_metrics['rmse'],color='black',linestyle='--',label='Dummy')
    axes[0].set(ylabel='Pooled OOF RMSE (Mg/ha)',title='Overall error',ylim=(0,max(dummy_metrics['rmse'],overall.rmse.max())*1.15))
    axes[0].legend()
    for label in labels:
        axes[1].plot(range(1,6),wide[label],marker='o',label=label)
    axes[1].set(xlabel='Spatial fold',ylabel='RMSE (Mg/ha)',title='Fold variation',xticks=range(1,6))
    axes[1].legend(fontsize=8)
    h=high[high.scope.eq('agbd_ge100')].set_index('feature_set').loc[labels]
    axes[2].bar(labels,h.mean_residual,color=['#737373','#3182bd','#31a354'])
    axes[2].axhline(0,color='black',linewidth=.7)
    axes[2].set(ylabel='Mean GEDI minus prediction (Mg/ha)', title='AGBD >=100: underprediction bias (n=56)')
    for ax in [axes[0],axes[2]]:
        ax.tick_params(axis='x',labelrotation=20,labelsize=8)
    fig.savefig(out/'feature_comparison.png',dpi=170)
    plt.close(fig)
    save_json(out/'validation.json',dict(status='passed',input_hashes_unchanged=hashes,
        oof_each_point_once=True,identical_test_ids=True,buffer_checks=fold_checks,
        baseline_max_prediction_difference=baseline_difference,repeat_max_prediction_difference=repeat[0],
        saved_metrics_recomputed=True,model_fits=16,new_downloads=False))
    print(overall.to_string(index=False))


if __name__ == '__main__':
    main()
