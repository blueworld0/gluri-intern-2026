"""Evaluate frozen Step 7 spatial folds. No network, tuning, or source mutation."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['MPLCONFIGDIR'] = str(Path(__file__).parent / '.mplconfig')
import numpy as np
import pandas as pd
import sklearn
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

FEATURES = ['ndvi', 'evi', 'ndmi']
MODELS = ['dummy_mean', 'week04_rf_baseline', 'week05_spatial_B_fixed']
LABELS = dict(zip(MODELS, ['Dummy', 'Week04 RF', 'Week05 Spatial RF']))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def metric(y, pred):
    y, pred = np.asarray(y, float), np.asarray(pred, float)
    return dict(n=len(y), rmse=float(np.sqrt(mean_squared_error(y, pred))),
                mae=float(mean_absolute_error(y, pred)), r2=float(r2_score(y, pred)))


def nearest(a, b):
    return np.concatenate([np.sqrt(((a[i:i+64, None] - b[None])**2).sum(axis=2)).min(axis=1)
                           for i in range(0, len(a), 64)])


def table(frame):
    def fmt(x):
        return f'{x:.4f}' if isinstance(x, (float, np.floating)) else str(x)
    return '| ' + ' | '.join(frame.columns) + ' |\n| ' + ' | '.join(['---']*len(frame.columns)) + ' |\n' + '\n'.join(
        '| ' + ' | '.join(fmt(x) for x in row) + ' |' for row in frame.itertuples(index=False, name=None))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).parent / 'outputs')
    out = parser.parse_args().output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise RuntimeError('Output directory must be empty; use a new --output-dir to preserve existing results.')
    src = ROOT / 'week06/step7/outputs'
    paths = [src / x for x in ['new441_evaluation_input.csv', 'split_assignments.csv',
             'split_summary.csv', 'evaluation_plan.json', 'validation.json', 'summary.json', 'report.md']]
    paths += [ROOT / x for x in ['week04/train_baseline.py', 'week05/improvement/experiment_design.json',
              'week05/improvement/model_comparison.md', 'week06/step5/outputs/vegetation_candidates.csv']]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    validation = dict(status='validating', input_sha256_before=hashes, errors=[],
                      block_size_m=5000, buffer_m=500, folds=5, crs='EPSG:32652',
                      source_data_modified=False, external_calls=False, hyperparameter_search=False)
    def check(condition, description):
        if not bool(condition):
            validation['errors'].append(description)
    data = pd.read_csv(paths[0], dtype={'shot_number': str}).set_index('point_id', drop=False)
    assignments = pd.read_csv(paths[1])
    summary = pd.read_csv(paths[2]).set_index('split_id')
    plan = json.loads(paths[3].read_text(encoding='utf-8'))
    check(json.loads(paths[4].read_text(encoding='utf-8'))['status'] == 'passed', 'Step 7 validation not passed')
    check(len(data) == 441 and data.index.is_unique, 'Expected 441 unique input points')
    check(set(data.cohort) == {'new441'}, 'Unexpected cohort')
    check(plan['primary_features'] == FEATURES and plan['target'] == 'agbd', 'Unexpected features/target')
    check(np.isfinite(data[FEATURES + ['agbd', 'easting_m', 'northing_m']]).all().all(), 'Nonfinite required values')
    check(not data.duplicated(['gedi_file','beam','shot_number']).any(), 'Duplicate GEDI shots')
    check(not data.duplicated(['latitude','longitude']).any(), 'Duplicate coordinates')
    check(not data.duplicated(['utm_10m_x','utm_10m_y']).any(), 'Duplicate 10m pixel')
    check((np.floor(data.easting_m / 5000) == data.grid_5km_x).all() and
          (np.floor(data.northing_m / 5000) == data.grid_5km_y).all(), '5km grid inconsistent')
    expected_group = 'X' + data.grid_5km_x.astype(str) + '_Y' + data.grid_5km_y.astype(str)
    check((data.spatial_group == expected_group).all(), 'Spatial group ID inconsistent')
    raw = pd.read_csv(paths[-1], dtype={'shot_number': str})
    check(len(raw) == 441, 'Step5 source count changed')
    for normalized, original in [('agbd','agbd_Mg_ha'), ('ndvi','ndvi'), ('evi','evi'), ('ndmi','ndmi')]:
        check(np.array_equal(data[normalized].to_numpy(), raw[original].to_numpy()), 'Step5 values/order mismatch: ' + normalized)
    expected = {
        'dummy_mean': {'class':'DummyRegressor','strategy':'mean'},
        'week04_rf_baseline': {'class':'RandomForestRegressor','n_estimators':100,'random_state':42,
                            'max_depth':None,'min_samples_leaf':1,'max_features':1.0,'n_jobs':1},
        'week05_spatial_B_fixed': {'class':'RandomForestRegressor','n_estimators':100,'random_state':42,
                                'max_depth':6,'min_samples_leaf':3,'max_features':0.67,'n_jobs':1}}
    for name in MODELS:
        check(plan['models'][name] == expected[name], 'Unexpected frozen parameters: ' + name)
    design = json.loads(paths[8].read_text(encoding='utf-8'))
    check(design['complexity_candidates']['depth6_leaf3_two_thirds'] ==
          {'max_depth':6,'min_samples_leaf':3,'max_features':0.67}, 'Week05 source settings mismatch')
    marked = data.evi_review_required.astype(str).str.lower().eq('true')
    check(marked.sum() == 1, 'Expected exactly one preflagged EVI point')
    outlier_id = data.index[marked][0]
    check(abs(data.loc[outlier_id, 'evi'] + 41.465836) < 1e-5, 'Unexpected EVI flagged value')
    fold_roles, details, heldout = {}, [], []
    for fold in range(1, 6):
        split = f'new441_spatial5_fold{fold}'
        a = assignments.loc[assignments.split_id.eq(split)].copy()
        check(len(a) == 441 and a.point_id.is_unique and set(a.point_id) == set(data.index), split + ': incomplete IDs')
        check(set(a.role) <= {'train','test','buffer_excluded'}, split + ': invalid role')
        check(set(a.cohort) == {'new441'}, split + ': cohort mismatch')
        check((a.spatial_group.to_numpy() == data.loc[a.point_id, 'spatial_group'].to_numpy()).all(), split + ': group mismatch')
        roles = {role:a.loc[a.role.eq(role), 'point_id'].tolist() for role in ['train','test','buffer_excluded']}
        fold_roles[fold] = roles
        tr, te, ex = [data.loc[roles[r]] for r in ['train','test','buffer_excluded']]
        check(set(tr.spatial_group).isdisjoint(te.spatial_group), split + ': train/test block overlap')
        check(set(ex.spatial_group).isdisjoint(te.spatial_group), split + ': excluded point in test block')
        check(set(te.index) == set(data.index[data.spatial_group.isin(te.spatial_group)]), split + ': incomplete test blocks')
        distances = nearest(tr[['easting_m','northing_m']].to_numpy(), te[['easting_m','northing_m']].to_numpy())
        check(distances.min() >= 500, split + ': buffer violation')
        if len(ex):
            check((nearest(ex[['easting_m','northing_m']].to_numpy(), te[['easting_m','northing_m']].to_numpy()) < 500).all(), split + ': incorrect excluded points')
        for role in roles:
            check(len(roles[role]) == int(summary.loc[split, role.replace('train','train_n').replace('test','test_n').replace('buffer_excluded','buffer_excluded_n')]), split + ': count mismatch')
        check(abs(distances.min() - summary.loc[split,'minimum_train_test_distance_m']) < 1e-6, split + ': distance mismatch')
        check(summary.loc[split,'buffer_m'] == 500, split + ': buffer setting mismatch')
        details.append(dict(fold=fold, train_n=len(tr), test_n=len(te), buffer_excluded_n=len(ex),
                            minimum_train_test_distance_m=float(distances.min()),
                            outlier_role=a.set_index('point_id').loc[outlier_id,'role']))
        heldout.extend(roles['test'])
    check(len(heldout) == 441 and len(set(heldout)) == 441, 'OOF coverage not exactly once')
    validation.update(data_n=len(data), spatial_group_n=int(data.spatial_group.nunique()),
                      required_values_finite=True, fold_details=details, outlier_point_id=outlier_id,
                      outlier_evi=float(data.loc[outlier_id,'evi']))
    validation['status'] = 'invalid' if validation['errors'] else 'inputs_passed'
    write_json(out / 'validation.json', validation)
    if validation['errors']:
        (out / 'report.md').write_text('# Step 8 검증 중단\n\n기존 설정을 변경하거나 학습하지 않았습니다.\n\n' + '\n'.join('- '+x for x in validation['errors']), encoding='utf-8')
        raise RuntimeError(validation['errors'])
    def model_for(name):
        params = {k:v for k,v in plan['models'][name].items() if k != 'class'}
        return DummyRegressor(**params) if name == 'dummy_mean' else RandomForestRegressor(**params)
    config = dict(features=FEATURES, target='agbd', models={n:model_for(n).get_params() for n in MODELS},
                  fold_source='week06/step7/outputs/split_assignments.csv',
                  cases=['main_all441','sensitivity_exclude_evi'], no_scaling_imputation_or_feature_selection=True,
                  random_state=42, n_jobs=1, thread_limits=1,
                  versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, sklearn=sklearn.__version__),
                  script_sha256=sha(Path(__file__)), input_sha256=hashes)
    write_json(out / 'experiment_config.json', config)
    predictions, fold_metrics, reproduction = [], [], []
    for case in config['cases']:
        for fold, roles in fold_roles.items():
            def select(role):
                return data.loc[[p for p in roles[role] if case == 'main_all441' or p != outlier_id]]
            train, test = select('train'), select('test')
            for name in MODELS:
                model = model_for(name)
                model.fit(train[FEATURES], train.agbd)
                pred = model.predict(test[FEATURES])
                frame = test[['point_id','agbd','latitude','longitude','spatial_group','evi']].copy()
                frame['case'], frame['fold'], frame['model'] = case, fold, name
                frame['prediction'], frame['residual'] = pred, test.agbd.to_numpy() - pred
                predictions.append(frame.reset_index(drop=True))
                fold_metrics.append(dict(case=case, fold=fold, model=name, train_n=len(train),
                                         buffer_excluded_n=len(select('buffer_excluded')), **metric(test.agbd, pred)))
                if case == 'main_all441' and fold == 1 and name != 'dummy_mean':
                    repeat = model_for(name).fit(train[FEATURES], train.agbd).predict(test[FEATURES])
                    reproduction.append(dict(model=name, fold=1, max_abs_prediction_difference=float(np.max(np.abs(pred-repeat)))))
            print(f'{case} fold {fold}: train={len(train)} test={len(test)}', flush=True)
    pred = pd.concat(predictions, ignore_index=True)
    fm = pd.DataFrame(fold_metrics)
    overall = pd.DataFrame([dict(case=c, model=m, **metric(g.agbd,g.prediction)) for (c,m),g in pred.groupby(['case','model'],sort=False)])
    sensitivity = []
    for name in MODELS:
        for scope in ['overall',1,2,3,4,5]:
            base = pred.loc[pred.case.eq('main_all441') & pred.model.eq(name)]
            exc = pred.loc[pred.case.eq('sensitivity_exclude_evi') & pred.model.eq(name)]
            if scope != 'overall':
                base, exc = base.loc[base.fold.eq(scope)], exc.loc[exc.fold.eq(scope)]
            common = base.loc[base.point_id.ne(outlier_id)]
            paired = common.merge(exc[['point_id','prediction']],on='point_id',suffixes=('_base','_refit'),validate='one_to_one')
            b, c, e = metric(base.agbd, base.prediction), metric(common.agbd,common.prediction), metric(exc.agbd,exc.prediction)
            row = dict(model=name,scope=scope, n_included=b['n'], n_common=c['n'], n_excluded=e['n'],
                       max_abs_common_prediction_change=float(np.abs(paired.prediction_base-paired.prediction_refit).max()))
            for key in ['rmse','mae','r2']:
                row.update({f'included_{key}':b[key], f'base_common_{key}':c[key], f'excluded_{key}':e[key],
                            f'evaluation_removal_delta_{key}':c[key]-b[key],
                            f'refit_common_delta_{key}':e[key]-c[key], f'total_delta_{key}':e[key]-b[key]})
            sensitivity.append(row)
    sens = pd.DataFrame(sensitivity)
    pred['agbd_bin'] = pd.cut(pred.agbd,[0,50,100,150,200,np.inf],right=False,
                              labels=['[0,50)','[50,100)','[100,150)','[150,200)','[200,inf)']).astype(str)
    bins = pd.DataFrame([dict(case=c,model=m,agbd_bin=b,**metric(g.agbd,g.prediction),
                             mean_agbd=float(g.agbd.mean()),mean_prediction=float(g.prediction.mean()),
                             mean_residual=float(g.residual.mean()), overprediction_fraction=float((g.residual < 0).mean()))
                         for (c,m,b),g in pred.groupby(['case','model','agbd_bin'],sort=False)])
    for name, frame in [('fold_metrics',fm),('overall_metrics',overall),('predictions',pred),
                        ('sensitivity_results',sens),('agbd_bin_metrics',bins),('fold_validation',pd.DataFrame(details))]:
        frame.to_csv(out / (name+'.csv'),index=False,encoding='utf-8')
    # Check saved predictions independently, including common test IDs across all models.
    saved = pd.read_csv(out / 'predictions.csv')
    for case, n in [('main_all441',441),('sensitivity_exclude_evi',440)]:
        groups = [saved.loc[saved.case.eq(case) & saved.model.eq(m)] for m in MODELS]
        check(all(len(g) == n and g.point_id.nunique() == n for g in groups), case + ': OOF coverage')
        check(all(set(zip(g.point_id,g.fold)) == set(zip(groups[0].point_id,groups[0].fold)) for g in groups), case + ': model test mismatch')
    for (c,m),g in saved.groupby(['case','model']):
        residual = g.agbd.to_numpy() - g.prediction.to_numpy()
        independent = [np.sqrt(np.mean(residual**2)),np.mean(np.abs(residual)),1-np.sum(residual**2)/np.sum((g.agbd-g.agbd.mean())**2)]
        stored = overall.loc[overall.case.eq(c)&overall.model.eq(m),['rmse','mae','r2']].to_numpy()[0]
        check(np.allclose(independent,stored,rtol=1e-12,atol=1e-12), 'Saved metrics mismatch')
    check((sens.loc[sens.scope.eq(2),'max_abs_common_prediction_change'] == 0).all(), 'Fold2 train unchanged but predictions differ')
    check(all(r['max_abs_prediction_difference'] == 0 for r in reproduction), 'Reproduction failed')
    hashes_after = {str(p.relative_to(ROOT)):sha(p) for p in paths}
    check(hashes == hashes_after, 'Input hashes changed')
    validation.update(status='passed' if not validation['errors'] else 'failed', training_performed=True,
                      oof_all_points_once=True, identical_test_points_across_models=True,
                      buffers_verified=True, block_disjointness_verified=True, duplicate_counts_zero=True,
                      sensitivity_original_fold_roles_preserved=True, fold2_common_predictions_identical=True,
                      metrics_independently_recomputed_from_saved_predictions=True, reproduction_checks=reproduction,
                      input_sha256_after=hashes_after, input_hashes_unchanged=hashes==hashes_after,
                      prediction_rows=len(saved), model_fits=32, created_at_utc=datetime.now(timezone.utc).isoformat())
    write_json(out / 'validation.json',validation)
    if validation['errors']:
        raise RuntimeError(validation['errors'])
    make_report(out,data,fm,overall,sens,bins,details,outlier_id,config)
    print(overall.to_string(index=False))


def make_report(out,data,fm,overall,sens,bins,details,outlier_id,config):
    main = overall.loc[overall.case.eq('main_all441')].sort_values('rmse')
    best = main.iloc[0]
    indexed = main.set_index('model')
    dummy = indexed.loc['dummy_mean','rmse']
    rf = indexed.loc['week04_rf_baseline','rmse']
    improved = indexed.loc['week05_spatial_B_fixed','rmse']
    fold = fm.loc[fm.case.eq('main_all441')]
    ranges = fold.groupby('model').agg(rmse_min=('rmse','min'),rmse_max=('rmse','max'),
                                      rmse_sample_sd=('rmse','std'), r2_min=('r2','min'),r2_max=('r2','max')).reset_index()
    sp = sens.loc[sens.scope.eq('overall'),['model','n_included','n_excluded','included_rmse','base_common_rmse',
                                         'excluded_rmse','evaluation_removal_delta_rmse','refit_common_delta_rmse','total_delta_rmse']]
    wide = fold.pivot(index='fold',columns='model',values='rmse')
    wins = {m:int((wide[m] < wide['dummy_mean']).sum()) for m in MODELS[1:]}
    sb = sens.loc[sens.scope.eq('overall')].set_index('model')
    outlier = pd.read_csv(out/'predictions.csv')
    outlier = outlier.loc[outlier.point_id.eq(outlier_id),['model','fold','agbd','prediction','residual']]
    text = f'''# Week06 Step 8 — 확대 데이터 기반 모델 성능 평가

## 목적과 기본 결과

산림으로 독립 검증되지 않은 Sentinel-2 SCL 식생 분류 통과 GEDI L4A 표본 441점만 평가했다. 기존 90점을 합치지 않았다. 기본 결과는 EVI 이상점도 포함한다. 전체 OOF RMSE가 가장 낮은 모델은 **{LABELS[best['model']]} ({best.rmse:.4f} Mg/ha)**이다. 최선이라는 표현은 이번 고정 비교의 관측 지표에 한정한다.

{table(main[['model','n','rmse','mae','r2']])}

Dummy 대비 RMSE 감소량은 Week04 RF **{dummy-rf:.4f} Mg/ha ({100*(dummy-rf)/dummy:.2f}%)**, Week05 Spatial RF **{dummy-improved:.4f} Mg/ha ({100*(dummy-improved)/dummy:.2f}%)**이다. 음수는 악화다. 기존 RF 대비 개선 RF의 RMSE 변화(개선−기존)는 **{improved-rf:+.4f} Mg/ha**이다. pooled OOF 지표는 각 점의 한 번의 평가 예측을 모아 직접 계산했으며 fold 지표 평균과 다르다.

## 고정 실험 설정과 검증

입력: `week06/step7/outputs/new441_evaluation_input.csv`. NDVI·EVI·NDMI → AGBD(Mg/ha). EPSG:32652의 원점 0 기준 5km 격자, 기존 5-fold 배정, 평가점 중심 500m 미만 학습점 제외를 그대로 재사용했다. 차단된 점을 다른 fold에서 임의로 재배정하지 않았다. 점·shot·좌표·공통 10m 셀 중복은 없고, 모든 441점이 정확히 한 번 평가된다. 학습/평가 블록은 겹치지 않고 모든 학습/평가점 거리가 500m 이상이다. 세 모델의 평가 ID는 fold 및 실험 조건마다 같다.

{table(pd.DataFrame(details))}

- Dummy: 각 fold 학습 AGBD의 평균만 예측한다.
- Week04 RF: n_estimators=100, max_depth=None, min_samples_leaf=1, max_features=1.0.
- Week05 개선 RF: **Spatial B**의 max_depth=6, min_samples_leaf=3, max_features=0.67을 고정했다. Week05 Random B(depth=3, leaf=4, max_features=1.0)는 적용하지 않았다. 근거는 `week05/improvement/model_comparison.md`, `experiment_design.json`, Step7 `evaluation_plan.json`이다.
- RF 공통 random_state=42, n_jobs=1. 나머지 sklearn 기본값은 `experiment_config.json`의 전체 get_params로 기록했다. 전체 데이터에 기반한 스케일링·결측 대치·변수 선택·신규 탐색은 없다. 저장된 기존 모델을 불러온 평가가 아니라 고정 train ID에서의 재학습이다.
- Python {config['versions']['python']}, numpy {config['versions']['numpy']}, pandas {config['versions']['pandas']}, sklearn {config['versions']['sklearn']}. 30회 본/민감도 적합과 fold1 RF 두 모델의 결정적 재실행 검증 2회, 총 32회 적합했다. 재실행 예측 차이는 0이다. 원본/분할 파일 SHA-256은 실행 전후 동일했다.

## fold별 편차

{table(fold[['model','fold','train_n','n','rmse','mae','r2']])}

{table(ranges)}

표본 수와 AGBD 분포가 다른 공간 fold 간 편차를 함께 본다. `rmse_sample_sd`는 5개 fold RMSE의 기술적 표준편차(ddof=1)이며 독립 반복 실험의 표준오차나 신뢰구간이 아니다. 음의 R²는 해당 평가집합 평균 기준 제곱오차보다 더 큰 오차를 뜻한다. Dummy 자체는 학습 평균을 쓰므로 R²=0일 필요가 없다.

## EVI 이상점 민감도

RMSE 기준 기존 RF는 5개 중 **{wins['week04_rf_baseline']}개**, 개선 RF는 **{wins['week05_spatial_B_fixed']}개** fold에서 Dummy보다 좋았다. 개선 RF는 기존 RF보다 모든 fold의 RMSE가 낮지만, fold4에서는 Dummy보다 나쁘다. 개선 RF의 fold2·4·5 R²는 음수다. 개선 RF의 MAE도 pooled 기준 낮지만 fold2·3·4에서는 기존 RF보다 높다. 따라서 지역 및 오차 지표 전반에서 안정적인 개선을 확보했다고 판단할 수 없다.

관측 ID: `{outlier_id}`. EVI={data.loc[outlier_id,'evi']:.6f}, AGBD={data.loc[outlier_id,'agbd']:.6f} Mg/ha, EVI 분모={data.loc[outlier_id,'evi_denominator']:.6f}. SCL=4가 지수 안정성을 보장하지 않는다. 값이 극단적이라는 이유만으로 물리적 오류라고 단정하지 않았다.

이 점은 fold2에서는 평가점이고 fold1/3/4/5에서는 학습점이다. 기본 실험은 포함, 민감도만 해당 ID를 모든 역할에서 제외했다. 원래 fold 및 buffer_excluded 배정은 유지했다. 따라서 fold2는 학습 349점/평가 78점, 다른 fold는 학습 1점 감소/평가 수 유지다. 500m 버퍼를 다시 산정하거나 이전 제외점을 복원하지 않았다.

기본 실험에서 이상점의 예측:

{table(outlier)}

{table(sp)}

`base_common_rmse`는 **재학습 없이** 기본 예측에서 해당 평가점만 뺀 공통 440점의 지표다. `evaluation_removal_delta_rmse`는 평가집합 변화 효과, `refit_common_delta_rmse`는 공통 440점에서 학습점 제외 후 재학습 효과다. 두 차이의 합은 전체 RMSE 변화와 같다. RMSE는 비선형 지표이므로 이 차이를 개별 관측점의 인과 기여로 해석하지 않는다. fold2의 공통 78점에서는 세 모델 모두 예측 변화가 정확히 0이며, 지표 변화는 평가점 제외 때문이다. 다른 fold의 변화는 학습에서 이 점을 제외한 영향이며 RF 부트스트랩 표본의 구성 변화도 함께 포함한다. 따라서 변화 전부를 EVI 수치 자체의 효과로 단정하지 않는다. 전체 및 fold별 RMSE·MAE·R² 분해는 `sensitivity_results.csv`에 저장했다. 민감도 결과로 기본 실험을 대체하거나 설정을 다시 선택하지 않았다.

## AGBD 구간별 편향

이상점 제외 후 기존 RF RMSE는 **{sb.loc['week04_rf_baseline','excluded_rmse']:.4f} Mg/ha**로 {sb.loc['week04_rf_baseline','total_delta_rmse']:+.4f} 변하고, 개선 RF는 **{sb.loc['week05_spatial_B_fixed','excluded_rmse']:.4f} Mg/ha**로 {sb.loc['week05_spatial_B_fixed','total_delta_rmse']:+.4f} 변한다. 개선 RF의 평가점 제거 효과는 {sb.loc['week05_spatial_B_fixed','evaluation_removal_delta_rmse']:+.4f}, 공통점 재학습 효과는 {sb.loc['week05_spatial_B_fixed','refit_common_delta_rmse']:+.4f} Mg/ha로 일부 상쇄된다. 이상점 제거만으로 성능 문제가 해결되지 않으며 기본 모델 순위도 유지된다.

잔차는 **관측 GEDI AGBD − 예측**이다. 양수는 과소 예측, 음수는 과대 예측이다. 다음은 기본 OOF 결과이며, 작은 구간 표본 수를 함께 확인해야 한다.

{table(bins.loc[bins.case.eq('main_all441'),['model','agbd_bin','n','mean_agbd','mean_prediction','mean_residual','rmse']])}

고 AGBD 구간의 표본 부족과 평균 회귀 경향을 확인하는 진단표다. 구간 경계는 성능 선택에 사용하지 않았으며 `[하한, 상한)`이다. 모든 조건의 상세 지표는 `agbd_bin_metrics.csv`에 있다.

## 해석의 한계와 다음 단계

실제로 50 Mg/ha 미만 표본이 311/441점(70.5%)이다. 개선 RF는 이 구간에서 평균 약 17.03 Mg/ha 과대 예측하며, 100–150 구간에서는 약 66.45, 150–200에서는 약 120.12, 200 이상에서는 약 178.30 Mg/ha 과소 예측한다. 200 이상은 4점뿐이며 이 구간의 RF RMSE는 두 모델 모두 Dummy보다 크다. 높은 AGBD를 구분하는 능력이 부족하고 평균 쪽으로 예측이 수축되는 문제가 남는다.

기존 90점은 Browser TIFF의 /65535 반사율 가정이 미검증이며, 신규 441점은 SAFE XML의 offset=-1000, quantification=10000을 확인해 처리했다. 품질 마스크 적용 기준도 다르다. 따라서 두 집합을 합치지 않았고 Week04/05의 기존 90점 성능 수치와 직접 개선률을 계산하지 않았다.

신규 표본도 기존과 동일한 2025-04-25·O36071 GEDI 관측에 속한다. 이번 평가는 독립 궤도 검증이 아니다. AGBD는 GEDI L4A 모델 추정치이며 독립 현장 실측값이 아니다. 대상은 산림 확정 표본이 아닌 SCL 식생 후보이며 농경지 등 다른 식생이 포함될 수 있다. 5km 격자·500m 버퍼는 공간 독립성을 보증하는 거리나 자기상관 범위 추정치가 아니다. 이 공간 CV를 독립 최종 일반화 성능이나 통계적 유의성의 증거로 해석하지 않는다.

다음 단계에서는 (1) EVI 분모·원본 픽셀과 주변 품질을 점검해 지수 안정성 처리 규칙을 사전에 정하고, (2) 고 AGBD 표본 부족 및 지역별 잔차를 점검하고, (3) 허용 시 독립 관측일·궤도와 검증된 산림 자료를 추가하며, (4) 기존 90점의 반사율/품질 조건을 통일한 뒤 결합 가능성을 재검토한다. RF가 Dummy보다 좋더라도 이 한계와 fold 편차는 남는다. 이번 단계에서는 추가 자료 확보나 튜닝을 수행하지 않았다.

## 재현 및 파일

저장소 루트의 PowerShell에서 다음을 실행한다. 출력 폴더에 파일이 있으면 덮어쓰지 않고 중단하므로 재현 시 새로운 경로를 지정한다.

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week06\\step8\\evaluate_models.py --output-dir .\\week06\\step8\\reproduction_outputs
```

의존성은 위 버전 및 기존 Week04 환경을 사용한다. 설치·네트워크 연결은 코드에서 하지 않는다. `fold_metrics.csv`/`overall_metrics.csv`/`predictions.csv`의 case는 기본 `main_all441`, 민감도 `sensitivity_exclude_evi`를 구분한다. 지표 단위는 RMSE·MAE·잔차 Mg/ha, R² 무차원이다. `validation.json`은 입력 해시, 분할 검증, 재현성 및 지표 독립 재계산 결과를 담는다. `experiment_config.json`은 전체 모델 파라미터와 버전을, `fold_validation.csv`는 학습/평가 수와 최소 거리를 기록한다.
'''
    (out/'report.md').write_text(text,encoding='utf-8')


if __name__ == '__main__':
    main()
