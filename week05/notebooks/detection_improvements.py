"""Colab runtime extension. Insert after baseline evaluation and before sensitivity.
Do not execute this file on the local Windows machine.
"""
from scipy.ndimage import gaussian_filter, maximum_filter, label as ndi_label
from skimage.segmentation import watershed
from rasterio.features import shapes as raster_shapes
from shapely.geometry import shape
import copy


def finalize_detection_frame(records, method):
    cols=DET_COLUMNS if 'DET_COLUMNS' in globals() else ['det_id','plotID','siteID','field_year','score','H_m','CD_m','geometry']
    frame=gpd.GeoDataFrame(records,geometry='geometry',crs=WORK_CRS) if len(records) else gpd.GeoDataFrame(columns=cols,geometry='geometry',crs=WORK_CRS)
    frame['method']=method
    frame['DBH_jucker_cm']=estimate_dbh(frame.H_m,frame.CD_m)
    frame['AGB_default_kg']=estimate_agb(frame.DBH_jucker_cm)
    frame['C_default_kg']=carbon_kg(frame.AGB_default_kg,CFG['carbon_fraction'])
    return frame


def deepforest_patch_candidates(detector,rgb,chm,plot,patch_size,cache_dir):
    # Cache one inference per patch, then filter the same raw predictions at all five thresholds.
    saved=[]
    for obj in [detector.model,getattr(detector.model,'model',None)]:
        if obj is not None and hasattr(obj,'score_thresh'):
            saved.append((obj,'score_thresh',obj.score_thresh)); obj.score_thresh=0.05
    configs=[]
    for obj in [detector.config,getattr(detector.config,'retinanet',None)]:
        if obj is None: continue
        try:
            old=obj['score_thresh']; obj['score_thresh']=0.05; configs.append((obj,old))
        except (KeyError,TypeError,AttributeError): pass
    if not saved and not configs: raise RuntimeError('DeepForest 내부 score_thresh를 0.05로 설정하지 못했습니다. 낮은 점수 탐색을 숨겨 진행하지 않습니다.')
    oldcfg=CFG.copy(); oldout=globals()['OUT']
    directory=Path(cache_dir)/f'patch_{patch_size}'
    directory.mkdir(parents=True,exist_ok=True)
    fingerprint=dict(rgb=str(rgb),rgb_bytes=Path(rgb).stat().st_size,patch_size=int(patch_size),
        patch_overlap=CFG['patch_overlap'],nms_iou=CFG['nms_iou'],min_score=0.1,
        model_version=importlib.metadata.version('deepforest'),model=CFG['model_name'],revision=CFG['model_revision'])
    stamp=directory/f'{plot.plotID}_input.json'
    rawpath=directory/f'{plot.plotID}_raw_predictions.csv'
    if rawpath.exists() and (not stamp.exists() or json.loads(stamp.read_text())!=fingerprint): rawpath.unlink()
    try:
        CFG['patch_size']=int(patch_size); CFG['score_threshold']=0.1
        globals()['OUT']=directory
        rows=predict_crowns(detector,rgb,chm,plot)
        stamp.write_text(json.dumps(fingerprint,indent=2))
    finally:
        CFG.clear(); CFG.update(oldcfg); globals()['OUT']=oldout
        for obj,key,value in saved: setattr(obj,key,value)
        for obj,value in configs: obj['score_thresh']=value
    for row in rows: row['det_id']=f'{plot.plotID}_p{patch_size}_{row["det_id"].split("_")[-1]}'
    return finalize_detection_frame(rows,f'DF_patch_{patch_size}')


def variable_window_peaks(height,valid,resolution_m,diameter_a=2.51,diameter_b=0.00901,min_height=3.0):
    # Popescu & Wynne (2004) eq.3: combined crown WIDTH 2.51503 + .00901 H^2.
    # User-requested rounded coefficient 2.51 is used, diameter in metres (not radius).
    mask=valid & (height>=min_height)
    cross=np.array([[0,1,0],[1,1,1],[0,1,0]],dtype=bool)
    # At 1m, the smallest requested radius is ~1.30m: axial neighbours are inside
    # every variable window, but diagonal pixels are not necessarily inside it.
    if resolution_m>(diameter_a+diameter_b*min_height**2)/2:
        raise ValueError('해상도가 최소 가변 창 반지름보다 큽니다. 더 세밀한 CHM 또는 현지 창 계수를 사용하세요.')
    maxima=mask & (height>=maximum_filter(np.where(valid,height,-np.inf),footprint=cross,mode='constant',cval=-np.inf))
    components,count=ndi_label(maxima,structure=cross.astype(int))
    candidates=[]
    for k in range(1,count+1):
        coords=np.argwhere(components==k)
        order=np.lexsort((coords[:,1],coords[:,0],-height[coords[:,0],coords[:,1]]))
        r,c=coords[order[0]]; H=float(height[r,c])
        radius=(diameter_a+diameter_b*H*H)/2
        pixels=int(np.ceil(radius/resolution_m))
        r0=max(0,r-pixels); r1=min(height.shape[0],r+pixels+1)
        c0=max(0,c-pixels); c1=min(height.shape[1],c+pixels+1)
        rr,cc=np.ogrid[r0:r1,c0:c1]
        circle=((rr-r)*resolution_m)**2+((cc-c)*resolution_m)**2<=radius**2
        neighborhood=height[r0:r1,c0:c1][circle & valid[r0:r1,c0:c1]]
        if len(neighborhood) and H>=float(neighborhood.max())-1e-7:
            candidates.append((int(r),int(c),H,float(radius*2)))
    return sorted(candidates,key=lambda x:(-x[2],x[0],x[1]))


def chm_watershed_crowns(chm_path,plot,sigma_m=1.0):
    with rasterio.open(chm_path) as ds:
        raw=ds.read(1,masked=True); transform=ds.transform; source_crs=ds.crs
        resolution=float(abs(transform.a))
        if not np.isclose(abs(transform.a),abs(transform.e)): raise ValueError('CHM 정사각형 픽셀 필요')
        rawdata=np.asarray(raw.astype(float).filled(np.nan),dtype=float)
        valid=np.isfinite(rawdata)&(rawdata>=0)
        # Normalized Gaussian prevents nodata from biasing the smoothed canopy edge.
        weight=gaussian_filter(valid.astype(float),sigma=sigma_m/resolution)
        smooth=gaussian_filter(np.where(valid,rawdata,0),sigma=sigma_m/resolution)/np.maximum(weight,1e-12)
        peaks=variable_window_peaks(smooth,valid,resolution,CFG['chm_window_a'],CFG['chm_window_b'],CFG['min_height_m'])
        markers=np.zeros(rawdata.shape,dtype=np.int32)
        accepted=[]
        for r,c,H,diameter in peaks:
            # H is the unsmoothed CHM height at the detected (smoothed) treetop pixel.
            if rawdata[r,c]<CFG['min_height_m']: continue
            k=len(accepted)+1; markers[r,c]=k
            accepted.append((k,r,c,float(rawdata[r,c]),diameter))
        if not accepted: return finalize_detection_frame([],'B_CHM')
        segments=watershed(-smooth,markers,mask=valid&(rawdata>=CFG['min_height_m']))
        crown_parts={}
        for geom,value in raster_shapes(segments.astype('int32'),mask=segments>0,transform=transform):
            crown_parts.setdefault(int(value),[]).append(shape(geom))
    rows=[]
    for k,r,c,H,diameter in accepted:
        xy=transform*(c+0.5,r+0.5)
        point=gpd.GeoSeries([Point(*xy)],crs=source_crs).to_crs(WORK_CRS).iloc[0]
        if not plot.geometry.covers(point): continue
        polygon=gpd.GeoSeries([unary_union(crown_parts[k])],crs=source_crs).to_crs(WORK_CRS).iloc[0]
        area=float(polygon.area)
        if area<=0: continue
        # Keep the existing point-in-BOX Hungarian matching rule for every method.
        rows.append(dict(det_id=f'{plot.plotID}_CHM_{k}',plotID=plot.plotID,siteID=plot.siteID,
            field_year=plot.field_year,imagery_year=plot.imagery_year,score=1.0,H_m=H,
            CD_m=2*np.sqrt(area/np.pi),crown_area_m2=area,crown_geometry_wkt=polygon.wkt,
            peak_x=point.x,peak_y=point.y,window_diameter_m=diameter,geometry=polygon.envelope))
    return finalize_detection_frame(rows,'B_CHM')


def combine_boxes_and_peaks(boxes,chm,dedup_m=1.0):
    rows=[]; used=set()
    for _,boxrow in boxes.sort_values(['score','det_id'],ascending=[False,True]).iterrows():
        candidates=chm[chm.apply(lambda p:boxrow.geometry.covers(Point(p.peak_x,p.peak_y)),axis=1)] if len(chm) else chm
        had_seed=len(candidates)>0
        candidates=candidates[~candidates.det_id.isin(used)]
        if had_seed and candidates.empty: continue # the same CHM seed was already assigned to a higher-score box
        row=boxrow.to_dict(); row['det_id']='C_'+str(row['det_id'])
        if len(candidates):
            peak=candidates.sort_values(['H_m','det_id'],ascending=[False,True]).iloc[0]
            used.add(peak.det_id); row.update(H_m=peak.H_m,peak_x=peak.peak_x,peak_y=peak.peak_y)
        else:
            # Preserve DF boxes without a CHM seed; do not invent a peak.
            row.update(peak_x=boxrow.geometry.centroid.x,peak_y=boxrow.geometry.centroid.y)
        point=Point(row['peak_x'],row['peak_y'])
        if any(point.distance(Point(x['peak_x'],x['peak_y']))<dedup_m for x in rows): continue
        row['fusion_source']='DeepForest_box'; rows.append(row)
    for _,p in chm.iterrows():
        point=Point(p.peak_x,p.peak_y)
        if any(g.covers(point) for g in boxes.geometry): continue
        if any(point.distance(Point(x['peak_x'],x['peak_y']))<dedup_m for x in rows): continue
        row=p.to_dict(); row['det_id']='C_'+row['det_id']; row['fusion_source']='outside_CHM_peak'; rows.append(row)
    return finalize_detection_frame(rows,'C_fusion')


def detection_counts(detections,truth,plot_ids):
    tp=fp=fn=0
    for pid in plot_ids:
        d=detections[detections.plotID==pid]; f=truth[truth.plotID==pid]
        hits=len(match_one_to_one(d,f)); tp+=hits; fp+=len(d)-hits; fn+=len(f)-hits
    return dict(TP=tp,FP=fp,FN=fn,precision=tp/(tp+fp) if tp+fp else np.nan,
                recall=tp/(tp+fn) if tp+fn else np.nan,F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan)


def nested_group_parameter_selection(candidates,plot_table,truth,n_splits=5):
    ids=np.array(sorted(plot_table.plotID.unique()))
    outer=GroupKFold(n_splits=min(n_splits,len(ids)))
    predictions=[]; selections=[]; inner_rows=[]; split_rows=[]
    keys=sorted(candidates)
    for fold,(train,test) in enumerate(outer.split(ids,groups=ids)):
        trainids=ids[train]; testids=ids[test]
        assert set(trainids).isdisjoint(testids)
        for pid in trainids: split_rows.append(dict(outer_fold=fold,plotID=pid,role='train'))
        for pid in testids: split_rows.append(dict(outer_fold=fold,plotID=pid,role='test'))
        inner=GroupKFold(n_splits=min(CFG['search_inner_folds'],len(trainids)))
        scores={key:dict(TP=0,FP=0,FN=0) for key in keys}
        for inner_fold,(itr,ival) in enumerate(inner.split(trainids,groups=trainids)):
            assert set(trainids[itr]).isdisjoint(trainids[ival])
            assert set(trainids[ival]).isdisjoint(testids)
            # Fixed pretrained inference has no fitting here; only held-out inner labels score parameters.
            for key in keys:
                result=detection_counts(candidates[key],truth,trainids[ival])
                inner_rows.append(dict(outer_fold=fold,inner_fold=inner_fold,patch_size=key[0],score_threshold=key[1],**result))
                for k in ['TP','FP','FN']: scores[key][k]+=result[k]
        def ranking(key):
            counts=scores[key]; den=2*counts['TP']+counts['FP']+counts['FN']
            f1=2*counts['TP']/den if den else -1
            return (-f1,key[0],key[1]) # deterministic ties
        best=min(keys,key=ranking)
        chosen=candidates[best][candidates[best].plotID.isin(testids)].copy()
        chosen['outer_fold']=fold; chosen['selected_patch_size']=best[0]; chosen['selected_score_threshold']=best[1]
        predictions.append(chosen)
        selections.append(dict(outer_fold=fold,patch_size=best[0],score_threshold=best[1],
            train_plots=';'.join(trainids),test_plots=';'.join(testids),**detection_counts(chosen,truth,testids)))
    oof=gpd.GeoDataFrame(pd.concat(predictions,ignore_index=True),geometry='geometry',crs=WORK_CRS)
    # Deployment selection uses all 23 plots; these scores are tuning scores, NOT unbiased test scores.
    deployment_rows=[]
    for key in keys: deployment_rows.append(dict(patch_size=key[0],score_threshold=key[1],**detection_counts(candidates[key],truth,ids)))
    deployment=pd.DataFrame(deployment_rows).sort_values(['F1','patch_size','score_threshold'],ascending=[False,True,True])
    return oof,pd.DataFrame(selections),pd.DataFrame(inner_rows),pd.DataFrame(split_rows),deployment


def evaluate_detection_method(det,plot_table,all_truth,coordinate_truth):
    metric_rows=[]; pair_tables={}; joined=None
    for population in ['visible','all']:
        pairs=[]
        for _,p in plot_table.iterrows():
            d=det[det.plotID==p.plotID]; f=coordinate_truth[coordinate_truth.plotID==p.plotID]
            if population=='visible': f=f[f.visible]
            match=match_one_to_one(d,f); pairs.append(match)
            metric_rows.append(dict(plotID=p.plotID,siteID=p.siteID,field_year=p.field_year,population=population,
                **detection_counts(d,f,[p.plotID])))
        pair_tables[population]=pd.concat(pairs,ignore_index=True)
    visible=pair_tables['visible']
    joined=visible.merge(det.drop(columns='geometry'),on=['det_id','plotID'],validate='one_to_one')
    joined=joined.merge(coordinate_truth.drop(columns=['geometry','plotID','siteID','field_year']),on='individualID',validate='one_to_one')
    joined['DBH_regression_oof_cm'],joined['cv_fold'],_=grouped_dbh_regression(joined,CFG['cv_folds'])
    joined['AGB_regression_oof_kg']=estimate_agb(joined.DBH_regression_oof_cm)
    joined['AGB_known_species_kg']=estimate_agb(joined.DBH_jucker_cm,joined.species_beta0,joined.species_beta1)
    joined['species_effect_kg']=joined.AGB_known_species_kg-joined.AGB_default_kg
    joined['species_effect_pct']=100*(joined.AGB_known_species_kg/joined.AGB_default_kg-1)
    plotrows=[]
    for _,p in plot_table.iterrows():
        ds=det[det.plotID==p.plotID]; fa=all_truth[all_truth.plotID==p.plotID]; fc=coordinate_truth[coordinate_truth.plotID==p.plotID]
        pairs=pair_tables['all']; pairs=pairs[pairs.plotID==p.plotID]
        plotrows.append(plot_biomass_totals(ds.to_dict('records'),fa.to_dict('records'),fc.to_dict('records'),pairs.to_dict('records'),
            p.to_dict(),CFG['carbon_fraction'],CFG['carbon_alternative']))
    return dict(detections=det,matched=joined,metrics=pd.DataFrame(metric_rows),pairs=pair_tables,plots=pd.DataFrame(plotrows))


def method_comparison_rows(method_results,requested_sites):
    rows=[]
    for method,result in method_results.items():
        for site in ['ALL']+list(requested_sites):
            pc=result['plots']; mt=result['matched']; dm=result['metrics']
            if site!='ALL': pc=pc[pc.siteID==site]; mt=mt[mt.siteID==site]; dm=dm[dm.siteID==site]
            stats=comparison_statistics(pc.C_ref_tC_ha,pc.C_est_tC_ha)
            row=dict(method=method,siteID=site,n_plots=len(pc),n_field=int(pc.n_truth.sum()),n_detected=int(pc.n_detected.sum()),
                n_carbon_plots=stats['n'],C_est_mean_tC_ha=stats['mean_estimate'],C_ref_mean_tC_ha=stats['mean_reference'],
                estimate_reference_ratio=stats['estimate_reference_ratio'],carbon_r=stats['r'],carbon_RMSE_tC_ha=stats['RMSE'],carbon_MAE_tC_ha=stats['MAE'])
            for pop in ['visible','all']:
                q=dm[dm.population==pop]; tp=q.TP.sum(); fp=q.FP.sum(); fn=q.FN.sum()
                row['F1_'+pop]=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan
                row['precision_'+pop]=tp/(tp+fp) if tp+fp else np.nan; row['recall_'+pop]=tp/(tp+fn) if tp+fn else np.nan
            for name,a,b in [('height_Bias_m','field_H_m','H_m'),('crown_Bias_m','field_CD_m','CD_m'),
                ('DBH_Bias_cm','field_DBH_cm','DBH_jucker_cm'),('individual_AGB_Bias_kg','ref_AGB_kg','AGB_default_kg')]:
                row[name]=errors(mt[a],mt[b])['Bias']
            rows.append(row)
    return pd.DataFrame(rows)


def detection_algorithm_smoke_checks():
    # Colab-only: actual scipy/skimage/rasterio calls on two synthetic trees, no external data/GPU.
    from types import SimpleNamespace
    yy,xx=np.mgrid[:64,:64]
    canopy=20*np.exp(-((xx-18)**2+(yy-32)**2)/8)+30*np.exp(-((xx-46)**2+(yy-32)**2)/8)
    path=OUT/'synthetic_CHM_smoke.tif'
    with rasterio.open(path,'w',driver='GTiff',width=64,height=64,count=1,dtype='float32',
        crs=WORK_CRS,transform=rasterio.transform.from_origin(0,64,1,1),nodata=-9999) as ds:
        ds.write(canopy.astype('float32'),1)
    p=SimpleNamespace(plotID='SMOKE',siteID='SMOKE',field_year=2021,imagery_year=2021,geometry=box(0,0,64,64))
    result=chm_watershed_crowns(path,p,1.0)
    assert len(result)==2,'합성 두 나무의 꼭대기 2개를 찾아야 합니다'
    assert np.allclose(sorted(result.H_m),[20,30],atol=0.01)
    assert np.allclose(result.CD_m,2*np.sqrt(result.crown_area_m2/np.pi))
    assert np.isclose(result.crown_area_m2.sum(),np.count_nonzero(canopy>=CFG['min_height_m']))
    first_peak=result.iloc[0].to_dict()
    duplicate=dict(first_peak); duplicate['det_id']='duplicate_box'; duplicate['score']=0.8
    first_peak['score']=0.9
    boxes=gpd.GeoDataFrame([first_peak,duplicate],geometry='geometry',crs=WORK_CRS)
    fused=combine_boxes_and_peaks(boxes,result,1.0)
    assert len(fused)==2,'동일 꼭대기 중복 상자는 제거하고 상자 밖 두 번째 꼭대기는 추가해야 합니다'
    ids=['P0','P1','P2','P3']; drows=[]; frows=[]
    for i,pid in enumerate(ids):
        x=i*20
        drows.append(dict(det_id='d'+pid,plotID=pid,score=0.8,H_m=20,CD_m=5,geometry=box(x-3,-3,x+3,3)))
        frows.append(dict(individualID='f'+pid,plotID=pid,stemDiameter=30,geometry=Point(x,0)))
    ds=gpd.GeoDataFrame(drows,geometry='geometry',crs=WORK_CRS)
    fs=gpd.GeoDataFrame(frows,geometry='geometry',crs=WORK_CRS)
    plotdf=pd.DataFrame({'plotID':ids})
    grid={(300,0.1):ds,(400,0.2):ds.iloc[:2].copy()}
    _,first,_,_,_=nested_group_parameter_selection(grid,plotdf,fs,n_splits=2)
    changed=fs.copy(); changed.loc[changed.plotID=='P0','geometry']=Point(9999,9999)
    _,second,_,_,_=nested_group_parameter_selection(grid,plotdf,changed,n_splits=2)
    fold=first[first.test_plots.str.split(';').apply(lambda z:'P0' in z)].index[0]
    assert first.loc[fold,['patch_size','score_threshold']].equals(second.loc[fold,['patch_size','score_threshold']]),'바깥 평가 정답이 선택에 새어 들어갔습니다'
    for _,r in first.iterrows(): assert set(r.train_plots.split(';')).isdisjoint(r.test_plots.split(';'))
    path.unlink(missing_ok=True)
    print('합성 CHM 꼭대기·watershed 면적·CD, GroupKFold 정답 격리 검증 통과')


# Same plot cohort for all methods; never silently remove an inconvenient plot.
detection_algorithm_smoke_checks()
if CFG['expected_plot_count'] is not None and len(plots)!=CFG['expected_plot_count']:
    raise RuntimeError(f'이번 비교는 동일 23개 조사구가 필요합니다. 현재 {len(plots)}개: 영상 커버·실측 선택 설정 확인')
plots[['plotID','siteID','field_year','imagery_year']].sort_values('plotID').to_csv(OUT/'comparison_cohort.csv',index=False)
BASELINE_RESULT=dict(detections=detections.copy(),matched=matched.copy(),metrics=detmetrics.copy(),
    pairs={k:v.copy() for k,v in MATCHES.items()},plots=plot_comparison.copy())
all_candidates={}; patch_frames={}; CHM_FRAMES=[]
for patch in CFG['search_patch_sizes']:
    frames=[]
    for _,p in plots.iterrows():
        print('A 후보 추론:',p.plotID,'patch',patch)
        frames.append(deepforest_patch_candidates(model,*RASTERS[p.plotID],p,patch,OUT/'detection_search_cache'))
    patch_frames[patch]=gpd.GeoDataFrame(pd.concat(frames,ignore_index=True),geometry='geometry',crs=WORK_CRS)
    for threshold in CFG['search_score_thresholds']:
        q=patch_frames[patch]; all_candidates[(patch,threshold)]=q[q.score>=threshold].copy()
A_oof,A_outer,A_inner,A_splits,A_deployment=nested_group_parameter_selection(all_candidates,plots,field_geo[field_geo.visible],CFG['cv_folds'])
A_outer.to_csv(OUT/'A_outer_groupkfold.csv',index=False); A_inner.to_csv(OUT/'A_inner_grid_scores.csv',index=False)
A_splits.to_csv(OUT/'A_fold_plot_assignments.csv',index=False); A_deployment.to_csv(OUT/'A_deployment_grid.csv',index=False)
fig,ax=plt.subplots(figsize=(7,5))
for patch,q in A_deployment.groupby('patch_size'):
    q=q.sort_values('score_threshold'); ax.plot(q.recall,q.precision,'o-',label=f'patch {patch}')
    for _,r in q.iterrows(): ax.annotate(str(r.score_threshold),(r.recall,r.precision),fontsize=8)
oof_counts=detection_counts(A_oof,field_geo[field_geo.visible],plots.plotID)
ax.scatter(oof_counts['recall'],oof_counts['precision'],marker='*',s=180,c='black',label='Nested outer OOF')
ax.set(xlabel='Recall',ylabel='Precision',title='A: tuning curves (selection scores) and outer OOF')
ax.legend(); fig.savefig(OUT/'A_precision_recall.png',dpi=170,bbox_inches='tight'); plt.close(fig)
for _,p in plots.iterrows():
    print('B CHM 꼭대기/watershed:',p.plotID)
    CHM_FRAMES.append(chm_watershed_crowns(RASTERS[p.plotID][1],p,CFG['chm_sigma_m']))
B_chm=gpd.GeoDataFrame(pd.concat(CHM_FRAMES,ignore_index=True),geometry='geometry',crs=WORK_CRS)
fusion=[]
for _,p in plots.iterrows():
    q=combine_boxes_and_peaks(A_oof[A_oof.plotID==p.plotID],B_chm[B_chm.plotID==p.plotID],CFG['fusion_dedup_m'])
    fusion.append(q)
C_fusion=gpd.GeoDataFrame(pd.concat(fusion,ignore_index=True),geometry='geometry',crs=WORK_CRS)
METHOD_RESULTS={'baseline_015':BASELINE_RESULT}
for name,det in [('A_DF_nested_OOF',A_oof),('B_CHM',B_chm),('C_fusion_nested_OOF',C_fusion)]:
    METHOD_RESULTS[name]=evaluate_detection_method(det,plots,field_all,field_geo)
method_comparison=method_comparison_rows(METHOD_RESULTS,CFG['sites'])
method_comparison.to_csv(OUT/'method_comparison.csv',index=False)
for name,result in METHOD_RESULTS.items():
    folder=OUT/'methods'/name; folder.mkdir(parents=True,exist_ok=True)
    for key,filename in [('detections','detected_trees.csv'),('matched','matched_trees_visible.csv'),('metrics','detection_metrics.csv'),('plots','plot_comparison.csv')]:
        q=result[key].copy()
        if 'geometry' in q: q['geometry']=q.geometry.apply(lambda g:g.wkt)
        pd.DataFrame(q).to_csv(folder/filename,index=False)
    for pop,q in result['pairs'].items(): q.to_csv(folder/f'matches_{pop}.csv',index=False)

for pid in CFG['comparison_overlay_plots']:
    if pid not in set(plots.plotID): raise ValueError(f'필수 나란히 비교 조사구 없음: {pid}')
    fig,axes=plt.subplots(1,4,figsize=(24,6),sharex=True,sharey=True)
    for ax,(name,result) in zip(axes,METHOD_RESULTS.items()):
        with rasterio.open(RASTERS[pid][0]) as ds: show(ds,ax=ax); crs=ds.crs
        q=result['detections']; q=q[q.plotID==pid].to_crs(crs)
        if len(q):
            q.boundary.plot(ax=ax,color='cyan',linewidth=0.6)
            if 'peak_x' in q and 'peak_y' in q:
                peaks=q.dropna(subset=['peak_x','peak_y'])
                points=gpd.GeoSeries([Point(r.peak_x,r.peak_y) for _,r in peaks.iterrows()],crs=WORK_CRS).to_crs(crs)
                ax.scatter(points.x,points.y,c='yellow',s=10,marker='x')
        f=field_geo[field_geo.plotID==pid].to_crs(crs)
        for visible,color in [(True,'lime'),(False,'magenta')]:
            z=f[f.visible==visible]; ax.scatter(z.geometry.x,z.geometry.y,c=color,s=14,edgecolors='black',linewidths=0.3)
        p=plots[plots.plotID==pid].iloc[0]
        gpd.GeoSeries([p.geometry],crs=WORK_CRS).to_crs(crs).boundary.plot(ax=ax,color='yellow')
        ax.set_title(f'{name}\nn={len(q)}')
    fig.suptitle(pid+' | green: visible field stems; magenta: occluded')
    fig.savefig(OUT/f'{pid}_methods_side_by_side.png',dpi=160,bbox_inches='tight'); plt.close(fig)

# Choose on plot carbon error, never fit/change Jucker or Chojnacky coefficients.
ranking=method_comparison[method_comparison.siteID=='ALL'].copy()
ranking=ranking[ranking.n_carbon_plots==len(plots)].sort_values([CFG['default_method_metric'],'method'])
if ranking.empty: raise RuntimeError('23개 전체 조사구의 참조 AGB가 유효하지 않아 기본 방법을 선택할 수 없습니다.')
DEFAULT_DETECTION_METHOD=ranking.iloc[0]['method']
deployment=A_deployment.iloc[0]
selection=dict(default_method=DEFAULT_DETECTION_METHOD,selection_metric=CFG['default_method_metric'],
    selection_scope='same 23 plot comparison; selection error is not independent generalization error',
    deployment_DF_patch_size=int(deployment.patch_size),deployment_DF_score_threshold=float(deployment.score_threshold),
    current_evaluation='outer GroupKFold OOF for A and C; globally selected parameters are deployment-only',
    chosen_carbon_RMSE_tC_ha=float(ranking.iloc[0].carbon_RMSE_tC_ha))
(OUT/'selected_detection_method.json').write_text(json.dumps(selection,ensure_ascii=False,indent=2))
chosen=METHOD_RESULTS[DEFAULT_DETECTION_METHOD]
detections=chosen['detections']; matched=chosen['matched']; detmetrics=chosen['metrics']; MATCHES=chosen['pairs']; plot_comparison=chosen['plots']
_,_,selected_dbh_fit=grouped_dbh_regression(matched,CFG['cv_folds'])
if selected_dbh_fit:
    selected_dbh_fit['detection_method']=DEFAULT_DETECTION_METHOD
    (OUT/'dbh_regression_full_fit.json').write_text(json.dumps(selected_dbh_fit,indent=2))
# Rebuild the existing stage table for the chosen default before downstream sensitivity/report cells.
stage_rows=[]
for year,z in matched.groupby('field_year'):
    for name,unit,refcol,predcol in comparisons:
        stage_rows.append(dict(stage=name,unit=unit,field_year=year,population='visible_matched',n_plots=z.plotID.nunique(),**errors(z[refcol],z[predcol])))
for (year,image_year),z in plot_comparison.groupby(['field_year','imagery_year']):
    for name,unit,a,b in [('plot_AGB','t/ha','AGB_ref_t_ha','AGB_est_t_ha'),('plot_carbon','tC/ha','C_ref_tC_ha','C_est_tC_ha')]:
        for site,q in [('ALL',z)]+list(z.groupby('siteID')):
            stage_rows.append(dict(stage=name,unit=unit,siteID=site,field_year=year,imagery_year=image_year,population='plot',n_plots=len(q),**comparison_statistics(q[a],q[b])))
for (year,pop),z in detmetrics.groupby(['field_year','population']):
    tp=z.TP.sum(); fp=z.FP.sum(); fn=z.FN.sum()
    stage_rows.append(dict(stage='detection',unit='ratio',field_year=year,population=pop,n=int(tp+fn),n_plots=z.plotID.nunique(),
        precision=tp/(tp+fp) if tp+fp else np.nan,recall=tp/(tp+fn) if tp+fn else np.nan,F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan))
stage_errors=pd.DataFrame(stage_rows)
display(method_comparison.round(3)); print('이후 단계 기본 방법:',DEFAULT_DETECTION_METHOD)
