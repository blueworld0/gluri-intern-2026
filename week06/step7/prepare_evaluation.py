"""Read-only comparison and evaluation design. No model import, fit or prediction."""
from pathlib import Path
import hashlib
import json
import math
import os
import re
from datetime import datetime,timezone

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/"outputs"
os.environ["MPLCONFIGDIR"]=str(HERE/".matplotlib")
import numpy as np
import pandas as pd

OLD=ROOT/"week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv"
NEW=ROOT/"week06/step5/outputs/vegetation_candidates.csv"
FEATURES=["ndvi","evi","ndmi"]


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def nearest(a,b):
    """Nearest distances in 64-row blocks; never constructs a full NxN matrix."""
    result=np.full(len(a),np.inf)
    for start in range(0,len(a),64):
        squared=((a[start:start+64,None,:]-b[None,:,:])**2).sum(axis=2)
        result[start:start+64]=np.sqrt(squared.min(axis=1))
    return result


def standardize(data,cohort):
    new=cohort=="new441"
    frame=pd.DataFrame({"point_id":data.point_id,"gedi_file":data.source_file if new else data.gedi_file,
        "beam":data.beam,"shot_number":data.shot_number.astype("string"),
        "observed_at_utc":data.observed_at_utc,"latitude":data.latitude,"longitude":data.longitude,
        "easting_m":data.easting_m,"northing_m":data.northing_m,
        "agbd":data.agbd_Mg_ha if new else data.agbd,"agbd_se":data.agbd_se_Mg_ha if new else data.agbd_se,
        **{name:data[name] for name in FEATURES},"sensitivity":data.sensitivity,
        "cohort":cohort,"forest_status":"not_assessed",
        "radiometry_status":"SAFE_XML_verified" if new else "Browser_65535_provisional",
        "s2_tile":"T52SDD" if new else "T52SDE_association_not_independently_proven",
        "scl_native_20m":data.scl_native_20m,"scl_aligned_10m":data.scl_aligned_10m})
    frame["orbit"]=frame.gedi_file.str.extract(r"_(O\d+)_",expand=False)
    frame["observation_date_utc"]=pd.to_datetime(frame.observed_at_utc,utc=True).dt.strftime("%Y-%m-%d")
    frame["grid_5km_x"]=np.floor(frame.easting_m/5000).astype(int)
    frame["grid_5km_y"]=np.floor(frame.northing_m/5000).astype(int)
    frame["spatial_group"]=[f"X{x}_Y{y}" for x,y in zip(frame.grid_5km_x,frame.grid_5km_y)]
    frame["utm_10m_x"]=np.floor(frame.easting_m/10).astype(int)
    frame["utm_10m_y"]=np.floor(frame.northing_m/10).astype(int)
    frame["evi_denominator"]=data.reflectance_B08+6*data.reflectance_B04-7.5*data.reflectance_B02+1 if new else np.nan
    frame["evi_review_required"]=frame.evi.abs().gt(2)|(new & frame.evi_denominator.abs().lt(.05))
    return frame


def md(table):
    return "\n".join(["| "+" | ".join(table.columns)+" |","| "+" | ".join(["---"]*len(table.columns))+" |"]
                      +["| "+" | ".join(str(x) for x in row)+" |" for row in table.itertuples(index=False,name=None)])


def run():
    if OUT.exists():raise FileExistsError("Existing Step7 output; refusing overwrite")
    OUT.mkdir(parents=True)
    hashes={str(path.relative_to(ROOT)):digest(path) for path in [OLD,NEW]}
    a=pd.read_csv(OLD,dtype={"shot_number":"string"})
    b=pd.read_csv(NEW,dtype={"shot_number":"string"})
    old,new=standardize(a,"old90"),standardize(b,"new441")
    allrows=pd.concat([old,new],ignore_index=True)
    assert len(old)==90 and len(new)==441 and len(allrows)==531
    required=["agbd","agbd_se",*FEATURES,"latitude","longitude","easting_m","northing_m"]
    missing=[]
    for cohort,p in [("old90",old),("new441",new)]:
        for col in required:
            v=pd.to_numeric(p[col],errors="coerce")
            missing.append({"cohort":cohort,"column":col,"missing":int(v.isna().sum()),"nonfinite":int((~np.isfinite(v)).sum())})
    pd.DataFrame(missing).to_csv(OUT/"required_quality.csv",index=False)
    assert sum(x["nonfinite"] for x in missing)==0
    overlap={"point_id_overlap":len(set(old.point_id)&set(new.point_id)),
        "shot_identity_overlap":len(set(map(tuple,old[["gedi_file","beam","shot_number"]].to_numpy())) & set(map(tuple,new[["gedi_file","beam","shot_number"]].to_numpy()))),
        "exact_coordinate_overlap":len(set(map(tuple,old[["latitude","longitude"]].to_numpy())) & set(map(tuple,new[["latitude","longitude"]].to_numpy()))),
        "common_utm_10m_cell_overlap":len(set(map(tuple,old[["utm_10m_x","utm_10m_y"]].to_numpy())) & set(map(tuple,new[["utm_10m_x","utm_10m_y"]].to_numpy()))),
        "within_cohort_duplicate_point_ids":int(allrows.duplicated(["cohort","point_id"]).sum()),
        "within_cohort_duplicate_shot_keys":int(allrows.duplicated(["cohort","gedi_file","beam","shot_number"]).sum())}
    oldxy=old[["easting_m","northing_m"]].to_numpy();newxy=new[["easting_m","northing_m"]].to_numpy()
    gaps=nearest(newxy,oldxy)
    overlap.update(minimum_distance_between_cohorts_m=float(gaps.min()),new_to_old_nearest_distance_median_m=float(np.median(gaps)))
    (OUT/"duplicate_and_distance_audit.json").write_text(json.dumps(overlap,indent=2),encoding="utf-8")
    assert not any(overlap[x] for x in ["point_id_overlap","shot_identity_overlap","exact_coordinate_overlap","common_utm_10m_cell_overlap","within_cohort_duplicate_point_ids","within_cohort_duplicate_shot_keys"])
    allrows.to_csv(OUT/"comparison_table_531_not_harmonized.csv",index=False)
    new.to_csv(OUT/"new441_evaluation_input.csv",index=False)
    allrows.loc[allrows.evi_review_required].to_csv(OUT/"evi_review_points.csv",index=False)
    mapping=[{"logical_field":x,"old90":x if x not in ["gedi_file","agbd","agbd_se"] else x,
              "new441":{"gedi_file":"source_file","agbd":"agbd_Mg_ha","agbd_se":"agbd_se_Mg_ha"}.get(x,x)}
             for x in ["point_id","gedi_file","shot_number","agbd","agbd_se",*FEATURES,"latitude","longitude","easting_m","northing_m"]]
    pd.DataFrame(mapping).to_csv(OUT/"schema_mapping.csv",index=False)
    descriptions=[]
    for cohort,p in [("old90",old),("new441",new)]:
        for col in ["agbd","agbd_se",*FEATURES,"sensitivity"]:
            v=p[col]
            descriptions.append({"cohort":cohort,"variable":col,"n":int(v.notna().sum()),"min":float(v.min()),
                "p05":float(v.quantile(.05)),"p25":float(v.quantile(.25)),"median":float(v.median()),
                "mean":float(v.mean()),"std":float(v.std()),"p75":float(v.quantile(.75)),"p95":float(v.quantile(.95)),"max":float(v.max())})
    stats=pd.DataFrame(descriptions);stats.to_csv(OUT/"distribution_comparison.csv",index=False)
    allrows.groupby(["cohort","gedi_file","observation_date_utc","orbit","beam"]).size().rename("count").reset_index().to_csv(OUT/"counts_cohort_date_orbit_beam.csv",index=False)
    allrows.groupby(["cohort","spatial_group"]).size().rename("count").reset_index().to_csv(OUT/"counts_spatial_groups.csv",index=False)
    bins=pd.cut(allrows.agbd,[-np.inf,50,100,150,200,np.inf],right=True)
    allrows.groupby(["cohort",bins],observed=False).size().rename("count").reset_index().to_csv(OUT/"agbd_bin_counts.csv",index=False)
    # Primary: coordinate-only 5km blocks, count-balanced 5 folds, 500m buffer.
    groups=new.groupby("spatial_group").size().reset_index(name="count").sort_values(["count","spatial_group"],ascending=[False,True])
    loads=[0]*5;assignment={}
    for r in groups.itertuples():
        fold=min(range(5),key=lambda k:(loads[k],k));assignment[r.spatial_group]=fold+1;loads[fold]+=r.count
    assignments=[];splitstats=[]
    def add_split(name,data,test_mask,buffer,primary):
        xy=data[["easting_m","northing_m"]].to_numpy()
        test=np.flatnonzero(test_mask);candidates=np.flatnonzero(~test_mask)
        distances=nearest(xy[candidates],xy[test])
        train=candidates[distances>=buffer];excluded=candidates[distances<buffer]
        assert len(train)>=30 and len(test)>=10
        minimum=float(nearest(xy[train],xy[test]).min())
        assert minimum>=buffer
        role=np.full(len(data),"train",dtype=object);role[test]="test";role[excluded]="buffer_excluded"
        assert not(set(data.iloc[train].point_id)&set(data.iloc[test].point_id))
        assignments.extend({"split_id":name,"point_id":row.point_id,"cohort":row.cohort,"spatial_group":row.spatial_group,"role":r}
                           for row,r in zip(data.itertuples(),role))
        splitstats.append({"split_id":name,"train_n":len(train),"test_n":len(test),"buffer_excluded_n":len(excluded),
                           "buffer_m":buffer,"minimum_train_test_distance_m":minimum,"primary":primary,
                           "train_evi_review_n":int(data.iloc[train].evi_review_required.sum()),"test_evi_review_n":int(data.iloc[test].evi_review_required.sum())})
    for fold in range(1,6):
        test=new.spatial_group.map(assignment).eq(fold).to_numpy()
        add_split(f"new441_spatial5_fold{fold}",new,test,500,True)
    permutation=np.random.RandomState(42).permutation(len(new))
    testmask=np.zeros(len(new),bool);testmask[permutation[:math.ceil(.2*len(new))]]=True
    add_split("new441_random80_20_seed42",new,testmask,0,False)
    # Regional transfer is a radiometrically confounded secondary experiment.
    add_split("legacy90_to_new441_exploratory",allrows,allrows.cohort.eq("new441").to_numpy(),500,False)
    splitdf=pd.DataFrame(assignments);splitdf.to_csv(OUT/"split_assignments.csv",index=False)
    split_stats=pd.DataFrame(splitstats);split_stats.to_csv(OUT/"split_summary.csv",index=False)
    primary=splitdf.loc[splitdf.split_id.str.startswith("new441_spatial5")]
    testcounts=primary.loc[primary.role.eq("test")].groupby("point_id").size()
    assert len(testcounts)==441 and testcounts.eq(1).all()
    for split,frame in primary.groupby("split_id"):
        assert not (set(frame.loc[frame.role.eq("test"),"spatial_group"]) & set(frame.loc[frame.role.eq("train"),"spatial_group"]))
    plan={"status":"prepared_not_trained","research_population":"Sentinel-2 SCL vegetation-screened GEDI L4A model-estimated AGBD samples; forest unverified",
        "forest_validation":"deferred_by_user; no EGIS or WorldCover acquisition in this step",
        "primary_dataset":"new441_evaluation_input.csv","primary_split":"new441_spatial5_fold1..5",
        "primary_features":FEATURES,"target":"agbd","target_unit":"Mg/ha",
        "models":{
            "dummy_mean":{"class":"DummyRegressor","strategy":"mean"},
            "week04_rf_baseline":{"class":"RandomForestRegressor","n_estimators":100,"random_state":42,"max_depth":None,"min_samples_leaf":1,"max_features":1.0,"n_jobs":1},
            "week05_spatial_B_fixed":{"class":"RandomForestRegressor","n_estimators":100,"random_state":42,"max_depth":6,"min_samples_leaf":3,"max_features":.67,"n_jobs":1},
            "week05_random_B_fixed_secondary":{"class":"RandomForestRegressor","n_estimators":100,"random_state":42,"max_depth":3,"min_samples_leaf":4,"max_features":1.0,"n_jobs":1}},
        "model_parameter_source":"week05/improvement/experiment_design.json and model_comparison.md; fixed historical settings, not newly selected",
        "metrics":["RMSE","MAE","R2"],"pooled_cv":"one held-out prediction per point; report pooled metrics and each fold; not mean R2 as pooled R2",
        "tuning":"none for primary frozen-parameter comparison; any later tuning must be entirely within buffered outer training data",
        "evi_policy":{"main":"retain all441; no silent deletion or clipping","sensitivity":"same frozen split IDs, remove only preflagged evi_review_required point from each role; report changed test counts and compare all models on same subset",
                      "review_rule":"abs(EVI)>2 or abs(SAFE EVI denominator)<0.05; heuristic diagnostic, not a validated error label"},
        "secondary_splits":["new441_random80_20_seed42","legacy90_to_new441_exploratory"],
        "pooled_531_training":"deferred until legacy radiometry and quality policy harmonized; 531 table is comparison-only",
        "orbit_holdout":"not_possible: all 531 from same O36071 and 2025-04-25",
        "preprocessing":"No full-data imputation/scaling/selection. Fit any learned transformation on training only.",
        "legacy_transfer_limits":"legacy EVI calibration provisional; quality policy and scene/tile differ; not an isolated geography effect or field validation",
        "existing_model_objects":"No trained model objects loaded; future comparisons require refitting on designated train IDs",
        "training_or_prediction_performed":False}
    (OUT/"evaluation_plan.json").write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding="utf-8")
    summary={"old_n":90,"new_n":441,"comparison_n":531,"input_columns_old":len(a.columns),"input_columns_new":len(b.columns),
        "missing_or_nonfinite_required":0,"duplicates":overlap,"unique_gedi_orbits":allrows.orbit.nunique(),
        "unique_gedi_dates":allrows.observation_date_utc.nunique(),"new_evi_review_points":int(new.evi_review_required.sum()),
        "old_evi_review_points":int(old.evi_review_required.sum()),"forest_confirmed_sample_count":None,
        "training_performed":False,"prediction_performed":False,"inputs_sha256":hashes,
        "input_hashes_unchanged":all(digest(ROOT/path)==value for path,value in hashes.items())}
    assert summary["input_hashes_unchanged"]
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    plot(old,new,splitdf)
    report(summary,stats,split_stats,new)
    (OUT/"validation.json").write_text(json.dumps({"status":"passed","source_counts":True,"duplicates_zero":True,
        "required_values_finite":True,"all441_tested_once_in_spatial_cv":True,"spatial_groups_disjoint":True,
        "buffers_at_least_500m":bool(split_stats.loc[split_stats.primary,"minimum_train_test_distance_m"].ge(500).all()),
        "split_ids_cover_each_dataset":True,"original_input_hashes_unchanged":True,"training_performed":False},indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True);print(split_stats.to_string(index=False),flush=True)


def plot(old,new,splits):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,4,figsize=(17,5),constrained_layout=True)
    for cohort,color in [(old,"#275dad"),(new,"#e37a21")]:
        label=f"{cohort.cohort.iloc[0]} (n={len(cohort)})"
        for ax,col in zip(axes[:3],["agbd","ndvi","evi"]):
            # ECDF preserves the EVI outlier without pretending it is typical.
            values=np.sort(cohort[col].to_numpy());ax.plot(values,np.arange(1,len(values)+1)/len(values),label=label,color=color)
            ax.set(xlabel=col,ylabel="Empirical cumulative proportion");ax.grid(alpha=.2)
        axes[3].scatter(cohort.easting_m,cohort.northing_m,s=9,color=color,label=label)
    axes[2].set_xscale("symlog",linthresh=1)
    axes[2].set_title("EVI: symlog scale includes -41.47")
    axes[0].legend(fontsize=7);axes[3].set(xlabel="UTM 52N easting (m)",ylabel="Northing (m)",title="Same orbit/date, distinct spatial cohorts")
    fig.savefig(OUT/"cohort_comparison.png",dpi=160);plt.close(fig)
    primary=splits.loc[splits.split_id.str.startswith("new441_spatial5")&splits.role.eq("test")]
    merged=new.merge(primary[["point_id","split_id"]],on="point_id",validate="one_to_one")
    fig,ax=plt.subplots(figsize=(7,7),constrained_layout=True)
    for split,g in merged.groupby("split_id"):ax.scatter(g.easting_m,g.northing_m,s=12,label=f"{split}: {len(g)}")
    ax.set(xlabel="UTM 52N easting (m)",ylabel="Northing (m)",title="Coordinate-only 5 km test groups; 500 m train buffer per fold")
    ax.legend(fontsize=7);ax.grid(alpha=.2);fig.savefig(OUT/"spatial_fold_map.png",dpi=160);plt.close(fig)


def report(s,stats,splits,new):
    selected=stats.loc[stats.variable.isin(["agbd","ndvi","evi","ndmi"]),["cohort","variable","mean","median","min","max"]].copy()
    for c in ["mean","median","min","max"]:selected[c]=selected[c].round(6)
    outlier=new.loc[new.evi_review_required]
    text=f'''# Week06 Step 7 — 확대 데이터 모델 평가 준비

## 연구 범위와 진행 상태

사용자 결정에 따라 EGIS 확보와 산림 검증은 후속 연구로 남긴다. 연구 대상을 **Sentinel-2 SCL 식생 분류를 통과한 GEDI L4A 추정 AGBD 표본**으로 정의한다. 산림 표본 또는 현장 실측 AGBD 정답이라고 표현하지 않는다. WorldCover도 이번에는 조회·다운로드·적용하지 않았다.

이번에는 구조·결측·중복·분포 비교, 고정된 분할 ID와 RF 평가 계획만 만들었다. **학습·예측·성능 점수 계산은 하지 않았다.** 원본 Week03~05 및 Step5 CSV는 수정하지 않았다.

## 표본 및 중복 검사

- 기존 90점: 원본 {s['input_columns_old']}열. agbd·agbd_se·gedi_file 사용.
- 신규 441점: 원본 {s['input_columns_new']}열. agbd_Mg_ha·agbd_se_Mg_ha·source_file 사용.
- 공통 컬럼으로 명시 매핑한 비교표는 **531행**이다. shot_number는 정밀도 손실 방지를 위해 문자열로 읽었다.
- 필수 지수·AGBD·표준오차·좌표의 결측·비유한값은 **0개**다.
- point_id, granule+BEAM+shot, 정확한 좌표, 공통 UTM 10m 셀의 두 집합 간 중복은 **모두 0개**다.
- 두 집합의 최소 중심점 거리는 **{s['duplicates']['minimum_distance_between_cohorts_m']:.2f}m**다. 원본 영상 격자가 서로 다르므로 공통 10m 셀 검사는 같은 원본 픽셀이라는 뜻과 구분한다.

## 지역·날짜·궤도 분포

**531점 모두 같은 GEDI granule·2025-04-25·궤도 O36071이다.** 신규 지역을 추가했지만 독립 관측일 또는 궤도를 추가한 것은 아니다. 따라서 궤도 홀드아웃 및 시간 일반화 평가는 현재 불가능하다. cohort·날짜·BEAM 및 5km 격자 수는 별도 CSV로 남겼다.

## 데이터 분포와 주의점

{md(selected)}

신규 AGBD 평균은 약 42.75 Mg/ha로 기존 약 103.94 Mg/ha보다 낮다. 이는 실제 표본 구성이 다르다는 기술통계이며 원인을 산림 유형 하나로 확정하지 않는다. 표준편차·분위수·AGBD 구간별 수도 보존했다. 이번 비교로 통계적 유의성이나 인과관계를 주장하지 않는다.

신규 표본 **1점의 EVI가 약 -41.47**이다. SAFE 반사율로 계산한 분모가 약 -0.0133으로 작아 발생한다. SCL=4라도 지수의 안정성이 보장되는 것은 아니다. `evi_review_points.csv`에 식별자를 남겼으며 441점 원본을 삭제하거나 EVI를 클리핑하지 않았다. 진단 규칙은 abs(EVI)>2 또는 abs(분모)<0.05로 명시했으며 물리적 오류를 확정하는 보편 기준이 아니다.

주평가에는 441점을 유지하고, 고정 분할을 그대로 둔 채 이 1점을 제외한 민감도 분석을 추가하는 계획을 준비했다. 각 역할에서 해당 점만 제외하고 바뀐 학습·평가 수를 보고한다. 모든 모델은 조건별 동일 평가점을 사용하며 성능이 좋은 조건만 골라 보고하지 않는다.

## 두 자료를 곧바로 합쳐 학습하지 않는 이유

기존 90점은 Browser TIFF의 /65535 가정으로 계산한 지수이며 2025년 변환식·장면 연결이 미검증이다. 신규 441점은 원본 SAFE XML로 offset=-1000, quantification=10000을 확인했다. 특히 EVI는 공통 배율에도 변하고 offset은 NDVI·NDMI에도 영향을 준다.

기존 90점은 SCL 양 격자 식생을 통과했지만 신규는 MSK_QUALIT·검출기·CLASSI까지 추가 적용했다. 엄격한 품질 기준도 다르다. 따라서 `comparison_table_531_not_harmonized.csv`는 **비교용**이며 방사·품질 기준을 맞춘 학습표가 아니다. 531점 혼합 학습은 기존 자료의 재처리 또는 검증 후로 보류한다.

## 학습·평가 분할 설계

**주평가: 신규 441점의 공간 5-fold CV.** EPSG:32652 원점 0 기준 5km×5km 격자로 나누고, 격자별 표본 수 내림차순·격자 ID 순으로 5개 fold에 표본 수가 균형을 이루도록 배정했다. 좌표와 표본 수만 사용하며 AGBD·지수·성능으로 분할을 선택하지 않았다. 각 fold 평가점에서 500m 미만인 다른 격자 점은 학습에서 제외한다. 모든 441점은 정확히 한 번 평가점이 된다.

{md(splits.round({'minimum_train_test_distance_m':3}))}

5km·500m는 탐색 연구를 위한 사전 지정값이며 공간 독립성을 보장하는 추정 자기상관 거리로 해석하지 않는다. 5개 fold도 독립된 신규 궤도나 최종 외부 검증이 아니다. 실제 분할은 `split_assignments.csv`로 고정했다.

무작위 80/20은 seed42의 비교용 보조 평가다. 가까운 점이 양쪽에 들어갈 수 있으므로 주 성능으로 삼지 않는다. 기존90→신규441 지역 전이 분할도 만들었으나 **방사 변환·품질·공간 분포 차이가 함께 섞이는 탐색 평가**다. 여기서 나온 차이를 순수 지역 일반화 효과로 해석하지 않는다. 기존 평가점이 Week04~05에서 반복 사용된 이력도 명시한다.

## 모델 재평가 준비

공통 입력은 ndvi·evi·ndmi, 목표는 agbd(Mg/ha)다. 주 비교는 학습 평균 Dummy, 기존 RF(100 trees, seed42), Week05 Spatial B의 고정 설정(max_depth=6, min_samples_leaf=3, max_features=0.67)이다. Week05 Random B 설정(depth3·leaf4·features1.0)은 보조 비교용으로 기록했다.

Week05의 설정은 **과거 자료에서 선택된 파라미터**를 고정 재사용하는 것이며 신규 자료에서 튜닝된 개선 모델이라고 부르지 않는다. 기존 학습 코드가 90행을 강제 검사하므로 그대로 신규441점에 실행하지 않는다. 저장된 모델 객체를 불러오거나 예측한 작업도 없으며, 다음 단계는 지정된 train ID로 모델을 새로 적합하는 재평가다.

RMSE·MAE·R²를 동일 test ID에서 계산하고, 주 공간 CV는 각 점의 단 한 번의 홀드아웃 예측을 모아 pooled 지표와 fold별 지표를 함께 보고하도록 계획했다. fold별 R² 평균을 pooled R²로 대체하지 않는다. 전체 자료의 평균 대체·변수 선택·튜닝은 하지 않는다. 추후 튜닝이 필요하면 outer train 내부에서 버퍼를 둔 공간 분할을 별도로 만들어야 한다.

## 주요 산출물

| 파일 | 용도 |
| --- | --- |
| comparison_table_531_not_harmonized.csv | 출처·방사 상태를 보존한 531점 비교표 |
| new441_evaluation_input.csv | 신규441점 공통 모델 컬럼과 EVI 검토 플래그 |
| schema_mapping.csv / required_quality.csv | 컬럼 대응 및 결측 검사 |
| duplicate_and_distance_audit.json | 중복과 두 집합의 거리 |
| distribution_comparison.csv / agbd_bin_counts.csv | 변수·AGBD 분포 |
| counts_cohort_date_orbit_beam.csv / counts_spatial_groups.csv | 지역·날짜·궤도·BEAM·격자 집계 |
| evi_review_points.csv | 원본 유지한 지수 검토 대상 |
| split_assignments.csv / split_summary.csv | 고정 학습·평가·버퍼 제외 ID와 수 |
| evaluation_plan.json | 모델 파라미터·지표·민감도·금지 해석 |
| summary.json / validation.json | 입력 해시·집계 및 분할 검증 |

![표본 비교](cohort_comparison.png)

![공간 분할](spatial_fold_map.png)

다음 단계에서는 이 고정 계획으로 신규441점의 RF·Dummy 평가와 EVI 민감도를 실행할 수 있다. 산림 검증을 하지 않았다는 한계와 단일 궤도·관측일 제약은 보고서와 모델 결과에 계속 명시한다.
'''
    (OUT/"report.md").write_text(text,encoding="utf-8")
    assert "???" not in text and "\ufffd" not in text


if __name__=="__main__":run()
