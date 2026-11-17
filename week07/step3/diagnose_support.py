"""Append training-support diagnostics to the existing feature comparison.

No model fitting or evaluation changes. Neighbors are descriptive probes only.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def md(df):
    def fmt(v):
        return f'{v:.3f}' if isinstance(v, (float, np.floating)) else str(v)
    return '| ' + ' | '.join(df.columns) + ' |\n| ' + ' | '.join(['---']*len(df.columns)) + ' |\n' + '\n'.join(
        '| ' + ' | '.join(fmt(v) for v in row) + ' |' for row in df.itertuples(index=False, name=None))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).parent/'outputs')
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    names = ['support_fold_targets.csv', 'support_neighbors.csv', 'support_point_diagnostics.csv',
             'support_summary.csv', 'support_scaling.csv', 'support_config.json', 'support_validation.json',
             'support_diagnostics.png', 'support_report.md']
    if any((out/n).exists() for n in names):
        raise RuntimeError('Support outputs already exist; use a new --output-dir.')
    base = Path(__file__).parent/'outputs'
    paths = dict(data=base/'evaluation_input.csv', predictions=base/'predictions.csv',
                 config=base/'experiment_config.json', validation=base/'validation.json',
                 assignments=ROOT/'week06/step7/outputs/split_assignments.csv',
                 top20=ROOT/'week07/step1/outputs/top20_absolute_errors.csv')
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths.values()}
    assert json.loads(paths['validation'].read_text(encoding='utf-8'))['status']=='passed'
    config = json.loads(paths['config'].read_text(encoding='utf-8'))
    features = config['feature_sets']
    data = pd.read_csv(paths['data'], dtype={'shot_number':str}).set_index('point_id', drop=False)
    pred = pd.read_csv(paths['predictions'])
    assignments = pd.read_csv(paths['assignments'])
    top = pd.read_csv(paths['top20'])[['point_id','error_rank']]
    assert len(data)==441 and data.index.is_unique
    assert np.isfinite(data[sum(features.values(), [])].to_numpy()).all()
    folds, summaries, neighbors, scales, point_rows = [], [], [], [], []
    seen = {name:[] for name in features}
    for fold in range(1,6):
        a = assignments[assignments.split_id.eq(f'new441_spatial5_fold{fold}')]
        assert len(a)==441 and a.point_id.is_unique and set(a.point_id)==set(data.index)
        tr = data.loc[sorted(a.loc[a.role.eq('train'),'point_id'])]
        te = data.loc[sorted(a.loc[a.role.eq('test'),'point_id'])]
        assert set(tr.index).isdisjoint(te.index) and set(tr.spatial_group).isdisjoint(te.spatial_group)
        distance = np.sqrt(((tr[['easting_m','northing_m']].to_numpy()[:,None]-
                            te[['easting_m','northing_m']].to_numpy()[None])**2).sum(axis=2))
        assert distance.min() >= 500
        folds.append(dict(fold=fold, train_n=len(tr), train_mean=float(tr.agbd.mean()),
            train_ge100_n=int(tr.agbd.ge(100).sum()), train_ge100_fraction=float(tr.agbd.ge(100).mean()),
            train_ge150_n=int(tr.agbd.ge(150).sum()), train_ge200_n=int(tr.agbd.ge(200).sum()),
            train_max=float(tr.agbd.max()), test_n=len(te), test_ge100_n=int(te.agbd.ge(100).sum()),
            test_ge100_fraction=float(te.agbd.ge(100).mean()), test_max=float(te.agbd.max()),
            test_above_train_max_n=int(te.agbd.gt(tr.agbd.max()).sum())))
        for name, cols in features.items():
            x = tr[cols].to_numpy()
            z = te[cols].to_numpy()
            center = np.median(x,axis=0)
            iqr = np.quantile(x,.75,axis=0)-np.quantile(x,.25,axis=0)
            sd = x.std(axis=0)
            scale = np.where(iqr>0,iqr,np.where(sd>0,sd,1))
            for j,col in enumerate(cols):
                scales.append(dict(fold=fold, feature_set=name, feature=col, train_median=center[j],
                                   train_iqr=iqr[j], scale=scale[j], fallback=bool(iqr[j]<=0)))
            # Train-only robust scaling prevents the single extreme EVI from
            # determining the distance scale. No model inputs are modified.
            difference = ((z-center)/scale)[:,None]-((x-center)/scale)[None]
            d = np.sqrt((difference**2).sum(axis=2))
            order = np.argsort(d,axis=1,kind='stable')[:,:5]
            pp = pred[pred.fold.eq(fold)&pred.feature_set.eq(name)].set_index('point_id')
            assert set(pp.index)==set(te.index)
            for i,row in enumerate(te.itertuples(index=False)):
                selected = tr.iloc[order[i]]
                assert len(selected)==5 and row.point_id not in set(selected.point_id)
                for rank, nr in enumerate(selected.itertuples(index=False),1):
                    neighbors.append(dict(feature_set=name, fold=fold, test_id=row.point_id,
                        neighbor_rank=rank, train_id=nr.point_id, train_agbd=nr.agbd,
                        train_spatial_group=nr.spatial_group, feature_distance=float(d[i,order[i,rank-1]])))
                actual_pred = float(pp.loc[row.point_id,'prediction'])
                assert np.isclose(pp.loc[row.point_id,'agbd'],row.agbd)
                point_rows.append(dict(feature_set=name,fold=fold,point_id=row.point_id,agbd=row.agbd,
                    prediction=actual_pred,residual=row.agbd-actual_pred,spatial_group=row.spatial_group,
                    train_max=float(tr.agbd.max()), above_train_max=bool(row.agbd>tr.agbd.max()),
                    outside_feature_range_n=int(((z[i]<x.min(axis=0))|(z[i]>x.max(axis=0))).sum()),
                    nearest_distance=float(d[i,order[i,0]]), neighbor5_mean_agbd=float(selected.agbd.mean()),
                    neighbor5_min_agbd=float(selected.agbd.min()), neighbor5_max_agbd=float(selected.agbd.max()),
                    neighbor5_ge100_n=int(selected.agbd.ge(100).sum()),
                    neighbor5_agbd_gap=float(row.agbd-selected.agbd.mean())))
                seen[name].append(row.point_id)
    detail = pd.DataFrame(point_rows).merge(top,on='point_id',how='left',validate='many_to_one')
    for name in features:
        assert len(seen[name])==len(set(seen[name]))==441
        g = detail[detail.feature_set.eq(name)]
        for scope, h in [('all441',g),('agbd_ge100',g[g.agbd.ge(100)]),('previous_top20',g[g.error_rank.notna()])]:
            summaries.append(dict(feature_set=name,scope=scope,n=len(h),
                above_train_max_n=int(h.above_train_max.sum()),
                outside_any_feature_range_n=int(h.outside_feature_range_n.gt(0).sum()),
                no_high_agbd_in_neighbor5_n=int(h.neighbor5_ge100_n.eq(0).sum()),
                neighbor5_high_agbd_fraction=float(h.neighbor5_ge100_n.mean()/5),
                mean_neighbor5_agbd=float(h.neighbor5_mean_agbd.mean()),
                mean_agbd=float(h.agbd.mean()),mean_prediction=float(h.prediction.mean()),
                mean_neighbor5_agbd_gap=float(h.neighbor5_agbd_gap.mean())))
    summary = pd.DataFrame(summaries)
    fold_df = pd.DataFrame(folds)
    nn = pd.DataFrame(neighbors)
    assert len(nn)==441*3*5 and not nn.test_id.eq(nn.train_id).any()
    # Reconstruct the local target summaries from the exported neighbor table.
    rebuilt=nn.groupby(['feature_set','test_id']).train_agbd.mean()
    aligned=detail.set_index(['feature_set','point_id']).neighbor5_mean_agbd
    rebuilt.index=rebuilt.index.set_names(aligned.index.names)
    assert np.allclose(rebuilt.sort_index(),aligned.sort_index())
    for filename, frame in [('support_fold_targets',fold_df),('support_neighbors',nn),
        ('support_point_diagnostics',detail),('support_summary',summary),('support_scaling',pd.DataFrame(scales))]:
        frame.to_csv(out/f'{filename}.csv',index=False)
    save = lambda name,value: (out/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    save('support_config.json',dict(k=5,high_agbd_threshold=100,
        distance='Euclidean after train-only median/IQR scaling; zero IQR -> train SD -> 1',
        feature_sets=features,input_sha256=hashes,script_sha256=sha(Path(__file__)),
        model_fit=False,feature_selection=False,evaluation_changed=False))
    high = summary[summary.scope.eq('agbd_ge100')]
    examples = detail[(detail.feature_set=='indices3')&detail.error_rank.le(5)].sort_values('error_rank')
    above = int(high.iloc[0].above_train_max_n)
    text = f'''## 고 AGBD 과소 예측의 학습 표본·입력 구분력 진단

추가 학습이나 튜닝 없이 동일한 평가점과 그 fold의 학습점만 사용했다. 100 Mg/ha 이상 기준은 기존 진단 구간을 유지했다. 전체 고 AGBD는 56/441점(12.7%)이며, 150 이상은 21점, 200 이상은 4점이다.

### fold별 고 AGBD 학습 지원

{md(fold_df)}

고 AGBD 평가점 56점 중 해당 fold 학습 최대 AGBD를 넘는 점은 {above}점이다. 따라서 학습 최댓값을 넘어선 외삽 제약만으로 56점 전체의 과소 예측을 설명할 수 없다. 표본 수·구성, 입력의 구분력, 모델 제약이 함께 작용할 가능성이 있으며 이 분석만으로 각각의 기여를 확정하지 않는다.

### 입력이 비슷한 학습점 5개와 비교

각 평가점에서 **그 fold 학습점만** 대상으로 입력별 median/IQR 표준화 후 유클리드 거리상 가까운 5점을 찾았다. IQR=0이면 학습 표준편차, 그것도 0이면 1을 사용했다. 모든 중심·척도는 학습점에서만 계산했으며 RF 입력이나 기존 예측에는 적용하지 않았다. k=5는 이번 기술적 진단을 위한 고정값이며 성능으로 탐색하지 않았다. 이웃 평균은 별도 예측 모델의 성능으로 평가하지 않는다.

{md(high)}

이번 결과에서 고 AGBD 56점의 평균은 143.34 Mg/ha인데 가까운 학습점 5개의 평균은 입력군별 54.89·57.81·55.52 Mg/ha였다. 56점 중 22~24점은 가까운 5개 안에 고 AGBD 학습점이 하나도 없었다. 고 AGBD 이웃 비율은 약 16%로 전체 평가점 이웃의 약 11%보다는 높아, 입력에 구분 정보가 전혀 없다는 뜻은 아니다. 다만 고 AGBD 예측을 지지하는 학습점이 충분히 모이지 않는 패턴이다. fold2는 학습 고 AGBD 비율 10.0%, 평가 24.1%로 구성 차이도 컸다. 현재 정보의 구분력 부족과 표본 구성 불균형을 함께 점검할 근거이며, 어느 하나를 단독 원인으로 확정하지 않는다.

`no_high_agbd_in_neighbor5_n`은 고 AGBD 평가점의 이웃 5개가 모두 100 미만인 경우의 수다. `outside_any_feature_range_n`은 하나 이상의 입력이 해당 학습 입력의 최솟값~최댓값 밖인 점 수이며, 범위 안이라고 다변량 지원이 충분하다는 의미는 아니다. 거리 크기는 입력 차원이 달라 입력군 사이 직접 비교하지 않는다.

### 기존 오차 상위 5점의 사례: 지수 3개 기준

{md(examples[['error_rank','agbd','prediction','train_max','neighbor5_mean_agbd','neighbor5_min_agbd','neighbor5_max_agbd','neighbor5_ge100_n']])}

이웃 선택은 AGBD를 보지 않고 입력값으로만 수행했다. 가까운 입력에서 다른 AGBD가 나타나는 것은 현재 표현의 구분력 한계와 일치하는 관측이지만, 식생지수 포화나 오표기·공간 위치 오차를 직접 입증하지 않는다. 고 AGBD 표본이 적다는 구성 효과도 섞여 있다. 상대적으로 가까운 5개라는 뜻이며 절대적으로 유사한 반사율·동일 피복임을 보증하지 않는다. 이웃이 여러 평가점에서 재사용되므로 독립 표본 검정으로 해석하지 않는다.

![학습 지원 진단](support_diagnostics.png)

### 이번 결과의 활용

기존 밴드 추가만으로 성능 개선이 작고 고 AGBD 과소 예측이 유지된다는 결과와 함께 해석한다. 이웃 표본의 AGBD 구성은 후속 자료 확보·모델 설계의 점검 근거로 사용하며, 이웃 평균이나 평가점 AGBD를 모델 입력으로 추가하지 않는다. 표본 재가중·고 AGBD 복제·타깃 변환 등 추가 조정을 현재 CV에 반복 적용해 최종 성능으로 주장하지 않는다.

새 자료 없이 할 수 있는 원인 진단을 진행했지만 원인을 확정하지는 못했다. 연구의 다음 우선순위는 별도 토지피복 자료로 산림·비산림을 구분하고 고 AGBD 표본의 피복 및 공간 대응을 확인하는 것이다. 독립 날짜·궤도 관측은 이후 일반화 검증에 필요하다. 이번 작업에서는 다운로드·신청·모델 변경을 하지 않았다.

### 재현과 검증

`support_*.csv`에 fold별 표본 수, 이웃 ID·거리·AGBD, 개별점 진단, 학습 스케일을 저장했다. 이웃은 기존 분할의 학습점에만 속하고 500m 버퍼를 충족한다. 3개 입력군 모두 441점이 정확히 한 번 포함되며, 저장 이웃표로 이웃 평균을 독립 재계산했다. 분석 입력과 기존 성능표·예측은 변경하지 않았다.

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\step3\\diagnose_support.py --output-dir .\\week07\\step3\\support_reproduction
```
'''
    (out/'support_report.md').write_text('# 고 AGBD 학습 지원 진단\n\n'+text,encoding='utf-8')
    fig, axes = plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    axes[0].bar(fold_df.fold-.18,fold_df.train_ge100_fraction,width=.36,label='Train')
    axes[0].bar(fold_df.fold+.18,fold_df.test_ge100_fraction,width=.36,label='Test')
    axes[0].set(xlabel='Fold',ylabel='Fraction with AGBD >=100',title='High AGBD representation',xticks=range(1,6))
    axes[0].legend()
    axes[1].bar(high.feature_set,high.no_high_agbd_in_neighbor5_n,color=['#737373','#3182bd','#31a354'])
    axes[1].set(ylabel='Count of high-AGBD test points',title='No high-AGBD point among 5 neighbors',ylim=(0,56))
    axes[1].tick_params(axis='x',labelrotation=20,labelsize=8)
    for name in features:
        h=detail[detail.feature_set.eq(name)&detail.agbd.ge(100)]
        axes[2].scatter(h.agbd,h.neighbor5_mean_agbd,s=17,alpha=.55,label=name)
    axes[2].plot([0,250],[0,250],'--',color='gray')
    axes[2].set(xlabel='Test GEDI AGBD (Mg/ha)',ylabel='5-neighbor training mean AGBD',
                title='Feature-neighbor target mismatch',xlim=(0,250),ylim=(0,250))
    axes[2].legend(fontsize=8)
    fig.savefig(out/'support_diagnostics.png',dpi=170)
    plt.close(fig)
    assert all(sha(ROOT/p)==h for p,h in hashes.items())
    save('support_validation.json',dict(status='passed',model_fits=0,new_downloads=False,
        unique_test_points_per_feature_set=441,neighbor_rows=len(nn),train_only_neighbors=True,
        buffer_verified=True,train_only_scaling=True,neighbor_means_recomputed=True,input_sha256_unchanged=hashes))
    # Only the default run extends the existing comparison report; reproduction
    # outputs remain isolated. The marker makes accidental duplicate appends fail.
    if out==base.resolve():
        report=base/'report.md'
        prior=report.read_text(encoding='utf-8')
        marker='<!-- support-diagnostics -->'
        assert marker not in prior
        report.write_text(prior+'\n'+marker+'\n\n'+text,encoding='utf-8')
    print(high.to_string(index=False))
    print(fold_df.to_string(index=False))


if __name__=='__main__':
    main()
