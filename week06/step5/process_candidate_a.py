"""Offline SAFE preprocessing and GEDI sampling; no training or downloads.

Refuses existing source/output directories. Only crops around actual GEDI
points are processed; source ZIP and previous weeks remain unchanged.
"""
from pathlib import Path, PurePosixPath
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
import shutil
import stat
import zipfile
import zlib
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["MPLCONFIGDIR"] = str(HERE / ".matplotlib")
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import rowcol, array_bounds
from rasterio.windows import Window, from_bounds
from rasterio.warp import reproject, transform

BANDS = ("B02", "B03", "B04", "B08", "B11", "B12")
PRODUCT = "S2C_MSIL2A_20250428T020711_N0511_R103_T52SDD_20250428T054014.SAFE"
ARCHIVE = ROOT / "week06/step4/downloads" / (PRODUCT + ".zip")
POINTS = ROOT / "week06/step3/candidate_a_20250425_gedi_points.csv"
OUT = HERE / "outputs"
SOURCES = HERE / "source"
NODATA = -9999.0


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def extract_safe():
    resume_empty = SOURCES.exists() and OUT.exists() and not any(p.is_file() for p in OUT.rglob("*"))
    if (SOURCES.exists() or OUT.exists()) and not resume_empty:
        raise FileExistsError("Step5 source/output already exists; use a new run directory; no overwrite performed.")
    verification = json.loads((ROOT / "week06/step4/checks/20261113T091309Z_download_cddb9b9a/download_verification.json").read_text())
    assert ARCHIVE.stat().st_size == verification["archive_bytes"]
    assert digest(ARCHIVE, "md5") == verification["catalogue_expected_md5"]
    with zipfile.ZipFile(ARCHIVE) as z:
        members = z.infolist()
        assert len({i.filename for i in members}) == len(members)
        if shutil.disk_usage(HERE).free < sum(i.file_size for i in members) + 5*1024**3:
            raise OSError("Insufficient extraction and processing free space")
        # Validate ALL paths before the first write.
        for item in members:
            parts = PurePosixPath(item.filename)
            assert not parts.is_absolute() and ".." not in parts.parts
            assert parts.parts[0] == PRODUCT and "\\" not in item.filename and ":" not in item.filename
            assert not stat.S_ISLNK(item.external_attr >> 16)
            target = SOURCES.joinpath(*parts.parts).resolve()
            assert target.is_relative_to(SOURCES.resolve())
        if resume_empty:
            # A failed pre-output run may reuse the complete extracted SAFE only
            # after checking each existing file against ZIP sizes and CRCs.
            expected={i.filename for i in members if not i.is_dir()}
            assert {p.relative_to(SOURCES).as_posix() for p in SOURCES.rglob("*") if p.is_file()} == expected
            for item in members:
                if item.is_dir():continue
                target=SOURCES.joinpath(*PurePosixPath(item.filename).parts)
                assert target.stat().st_size==item.file_size
                crc=0
                with target.open("rb") as f:
                    for block in iter(lambda:f.read(8*1024*1024),b""):crc=zlib.crc32(block,crc)
                assert crc==item.CRC
            print("Reused fully CRC-verified extracted SAFE; existing output contains no files.",flush=True)
            return SOURCES/PRODUCT
        SOURCES.mkdir(exist_ok=False)
        for item in members:
            target = SOURCES.joinpath(*PurePosixPath(item.filename).parts)
            if item.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(item) as src, target.open("xb") as dst:
                    shutil.copyfileobj(src, dst, 1024*1024)
    print("SAFE extracted to new Step5 source directory.", flush=True)
    OUT.mkdir(exist_ok=False)
    (OUT / "rasters").mkdir()
    return SOURCES / PRODUCT


def crop(path, bounds, indexes=1):
    with rasterio.open(path) as ds:
        win = from_bounds(*bounds, ds.transform)
        left = max(0, int(np.floor(win.col_off))-2)
        top = max(0, int(np.floor(win.row_off))-2)
        right = min(ds.width, int(np.ceil(win.col_off+win.width))+2)
        bottom = min(ds.height, int(np.ceil(win.row_off+win.height))+2)
        win = Window(left, top, right-left, bottom-top)
        data = ds.read(indexes, window=win)
        profile = {"crs": ds.crs, "transform": ds.window_transform(win),
                   "height": int(win.height), "width": int(win.width)}
        return data, profile


def warp(data, src, target, method, nodata=None):
    # Zero is VALID for quality/cloud flags. Reserving destination NoData=0
    # makes GDAL alter valid zeros to 1; reserve 255 for categorical data.
    dst_nodata = np.nan if np.issubdtype(data.dtype,np.floating) else 255
    dest = np.full((target["height"], target["width"]), dst_nodata, dtype=data.dtype)
    reproject(data, dest, src_crs=src["crs"], src_transform=src["transform"],
              dst_crs=target["crs"], dst_transform=target["transform"], resampling=method,
              src_nodata=nodata, dst_nodata=dst_nodata)
    return dest


def read_grid(path, target, indexes=1):
    data, profile = crop(path, array_bounds(target["height"],target["width"],target["transform"]), indexes)
    assert profile["crs"] == target["crs"]
    if data.ndim == 2:
        return warp(data, profile, target, Resampling.nearest)
    return np.stack([warp(a, profile, target, Resampling.nearest) for a in data])


def save(name, data, profile, tags, categorical=False):
    path = OUT / "rasters" / (name + ".tif")
    dtype = "uint8" if categorical else "float32"
    with rasterio.open(path, "w", driver="GTiff", width=profile["width"], height=profile["height"],
                       crs=profile["crs"], transform=profile["transform"], count=1, dtype=dtype,
                       nodata=0 if categorical else NODATA, compress="deflate", tiled=True) as ds:
        ds.write(data.astype("uint8") if categorical else np.where(np.isfinite(data),data,NODATA).astype("float32"),1)
        ds.update_tags(**tags)
        if not categorical: ds.set_band_unit(1,"reflectance" if name.startswith("reflectance") else "unitless")


def md(table):
    return "\n".join(["| " + " | ".join(table.columns) + " |", "| " + " | ".join(["---"]*len(table.columns)) + " |"]
                      + ["| " + " | ".join(str(x) for x in row) + " |" for row in table.itertuples(index=False,name=None)])


def run():
    baseline_hash = digest(ROOT / "week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv")
    points_hash = digest(POINTS)
    safe = extract_safe()
    product_xml = safe / "MTD_MSIL2A.xml"
    tree = ET.parse(product_xml)
    local = lambda tag: tag.rsplit("}",1)[-1]
    def vals(tree, tag): return [e.text for e in tree.iter() if local(e.tag)==tag]
    assert vals(tree,"PRODUCT_URI") == [PRODUCT]
    q = float(vals(tree,"BOA_QUANTIFICATION_VALUE")[0])
    ids = {e.attrib["bandId"]: e.attrib["physicalBand"] for e in tree.iter() if local(e.tag)=="Spectral_Information"}
    offsets = {ids[e.attrib["band_id"]]:float(e.text) for e in tree.iter() if local(e.tag)=="BOA_ADD_OFFSET"}
    # XML physicalBand uses B2/B3/B4/B8 rather than padded filenames.
    offsets = {b: offsets["B"+str(int(b[1:]))] for b in BANDS}
    nodata_dn, saturated_dn = None,None
    for e in tree.iter():
        if local(e.tag)=="Special_Values":
            v = {local(c.tag):c.text for c in e}
            if v.get("SPECIAL_VALUE_TEXT")=="NODATA": nodata_dn=int(v["SPECIAL_VALUE_INDEX"])
            if v.get("SPECIAL_VALUE_TEXT")=="SATURATED": saturated_dn=int(v["SPECIAL_VALUE_INDEX"])
    assert nodata_dn == 0 and saturated_dn == 65535
    granules = list((safe/"GRANULE").iterdir())
    assert len(granules)==1
    granule=granules[0]
    tiletree=ET.parse(granule/"MTD_TL.xml")
    tile_time = vals(tiletree,"SENSING_TIME")[0]
    catalogue_time="2025-04-28T02:07:11.025Z"
    p=pd.read_csv(POINTS,dtype={"shot_number":"string"})
    assert len(p)==1055 and p.point_id.nunique()==1055 and p.quality_valid.all()
    paths={b:next((granule/"IMG_DATA"/("R20m" if b in {"B11","B12"} else "R10m")).glob(f"*_{b}_*.jp2")) for b in BANDS}
    with rasterio.open(paths["B04"]) as ds:
        assert ds.crs.to_epsg()==32652 and ds.res==(10,10)
        x,y=map(np.asarray,transform("EPSG:4326",ds.crs,p.longitude.tolist(),p.latitude.tolist()))
        requested=(x.min()-1000,y.min()-1000,x.max()+1000,y.max()+1000)
        win=from_bounds(*requested,ds.transform)
        left=max(0,int(np.floor(win.col_off)));top=max(0,int(np.floor(win.row_off)))
        right=min(ds.width,int(np.ceil(win.col_off+win.width)));bottom=min(ds.height,int(np.ceil(win.row_off+win.height)))
        w=Window(left,top,right-left,bottom-top)
        ref={"crs":ds.crs,"transform":ds.window_transform(w),"width":int(w.width),"height":int(w.height)}
    bounds=array_bounds(ref["height"],ref["width"],ref["transform"])
    rr,cc=map(np.asarray,rowcol(ref["transform"],x,y))
    assert ((rr>=0)&(rr<ref["height"])&(cc>=0)&(cc<ref["width"])).all()
    p["easting_m"],p["northing_m"]=x,y
    p["raster_row"],p["raster_col"]=rr,cc
    p["sentinel_product"]=PRODUCT
    p["sentinel_tile"]="T52SDD"
    p["sentinel_tile_sensing_time_utc"]=tile_time
    p["sentinel_catalogue_start_utc"]=catalogue_time
    observed=pd.to_datetime(p.observed_at_utc,utc=True)
    p["time_gap_hours_tile_metadata"]=(pd.Timestamp(tile_time)-observed).dt.total_seconds()/3600
    p["forest_status"]="not_assessed"
    tags={"product":PRODUCT,"tile_sensing_time_utc":tile_time,"processing_baseline":vals(tree,"PROCESSING_BASELINE")[0],
          "radiometry":"(DN + per-band BOA_ADD_OFFSET) / BOA_QUANTIFICATION_VALUE", "quantification":str(q),
          "source":"native SAFE JP2 and XML; no Browser export scaling"}
    arrays={}; all_valid=np.ones((ref["height"],ref["width"]),bool); bandstats=[]; qualitystats=[]
    for band in BANDS:
        raw,native=crop(paths[band],bounds)
        quality_path=granule/"QI_DATA"/f"MSK_QUALIT_{band}.jp2"
        detector_path=granule/"QI_DATA"/f"MSK_DETFOO_{band}.jp2"
        quality=read_grid(quality_path,native,list(range(1,9)))
        assert set(np.unique(quality)).issubset({0,1})
        detector=read_grid(detector_path,native)
        assert np.isin(detector,np.arange(13)).all()
        quality_bad=np.any(quality!=0,axis=0)
        valid=(raw!=nodata_dn)&(raw!=saturated_dn)&~quality_bad&(detector>0)
        rho=(raw.astype("float32")+offsets[band])/q
        rho[~valid]=np.nan
        aligned=warp(rho,native,ref,Resampling.bilinear,nodata=np.nan)
        # Require all bilinear contributors to be valid; do not interpolate across bad source pixels.
        support=warp(valid.astype("float32"),native,ref,Resampling.bilinear)
        aligned_valid=np.isfinite(aligned)&(support>=1-1e-6)
        aligned[~aligned_valid]=np.nan
        arrays[band]=aligned;all_valid &= aligned_valid
        p[f"reflectance_{band}"]=aligned[rr,cc]
        p[f"{band}_valid"]=aligned_valid[rr,cc]
        nr,nc=map(np.asarray,rowcol(native["transform"],x,y))
        p[f"{band}_native_dn"]=raw[nr,nc]
        p[f"{band}_native_detector_id"]=detector[nr,nc]
        for j,label in enumerate(["anc_lost","anc_degraded","msi_lost","msi_degraded","defective","nodata","partial_crosstalk","saturated"]):
            p[f"{band}_qi_{label}"]=quality[j,nr,nc].astype(int)
            qualitystats.append({"band":band,"mask_layer":j+1,"meaning":label,"gedi_flagged":int(quality[j,nr,nc].sum())})
        finite=aligned[np.isfinite(aligned)]
        bandstats.append({"band":band,"native_resolution_m":20 if band in {"B11","B12"} else 10,
                          "offset":offsets[band],"quantification":q,"aligned_valid_pixels":len(finite),
                          "minimum":float(finite.min()),"maximum":float(finite.max()),"median":float(np.median(finite)),
                          "negative_pixels":int((finite<0).sum()),"above_one_pixels":int((finite>1).sum()),
                          "gedi_valid":int(p[f"{band}_valid"].sum())})
        save("reflectance_"+band,aligned,ref,dict(tags,band_offset=str(offsets[band])))
        print(f"{band}: sampled valid={int(p[f'{band}_valid'].sum())}/1055",flush=True)
    p=p.copy()  # Consolidate quality columns before adding index/scene fields.
    # Load the existing pure index module directly, without importing Week03 main/plotting.
    spec=importlib.util.spec_from_file_location("existing_indices",ROOT/"week03/src/vegetation_indices.py")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    indices,missing=module.calculate(arrays);assert not missing
    for name,data in indices.items():
        data[~all_valid]=np.nan
        save(name,data,ref,tags)
        p[name]=data[rr,cc]
    scl_path=next((granule/"IMG_DATA/R20m").glob("*_SCL_20m.jp2"))
    scl=read_grid(scl_path,ref)
    assert np.isin(scl,np.arange(12)).all()
    p["scl_aligned_10m"]=scl[rr,cc]
    with rasterio.open(scl_path) as ds:
        p["scl_native_20m"]=[int(v[0]) for v in ds.sample(zip(x,y))]
    p["scl_boundary_disagreement"]=p.scl_native_20m.ne(p.scl_aligned_10m)
    save("scl",scl,ref,tags,categorical=True)
    classi_path=granule/"QI_DATA/MSK_CLASSI_B00.jp2"
    classi=read_grid(classi_path,ref,[1,2,3])
    classi_bad=np.any(classi!=0,axis=0)
    for j,label in enumerate(["opaque_cloud","cirrus","snow_ice"]):p[f"classi_{label}"]=classi[j,rr,cc].astype(int)
    for label in ["CLDPRB","SNWPRB"]:
        path=granule/"QI_DATA"/f"MSK_{label}_20m.jp2"
        prob=read_grid(path,ref)
        p[label.lower()+"_percent"]=prob[rr,cc]
    p["all_bands_valid"]=p[[f"{b}_valid" for b in BANDS]].all(axis=1)
    p["indices_finite"]=np.isfinite(p[["ndvi","evi","ndmi"]]).all(axis=1)
    p["ndvi_ndmi_in_range"]=p.ndvi.between(-1,1)&p.ndmi.between(-1,1)
    p["evi_outside_minus1_plus1"]=~p.evi.between(-1,1)&p.evi.notna()
    p["evi_small_denominator"]=abs(p.reflectance_B08+6*p.reflectance_B04-7.5*p.reflectance_B02+1)<.05
    p["scl_cloud_shadow_either_grid"]=p.scl_native_20m.isin([2,3,8,9,10])|p.scl_aligned_10m.isin([2,3,8,9,10])
    p["clear_surface_both_scl_grids"]=p.scl_native_20m.isin([4,5,6])&p.scl_aligned_10m.isin([4,5,6])
    p["classi_clear"]=~classi_bad[rr,cc]
    p["time_gap_within_168h"]=p.time_gap_hours_tile_metadata.abs()<=168
    p["clear_valid_index_candidate"]=(p.all_bands_valid&p.indices_finite&p.ndvi_ndmi_in_range
        &p.clear_surface_both_scl_grids&p.classi_clear&p.time_gap_within_168h)
    p["vegetation_candidate"]=(p.clear_valid_index_candidate&p.scl_native_20m.eq(4)&p.scl_aligned_10m.eq(4))
    p["pixel_observation_count"]=p.groupby(["raster_row","raster_col"]).point_id.transform("size")
    p["duplicate_pixel"]=p.pixel_observation_count.gt(1)
    reason=[]
    for row in p.itertuples():
        reasons=[]
        if not row.all_bands_valid:reasons.append("band_nodata_quality_or_invalid_interpolation_support")
        if not row.indices_finite:reasons.append("nonfinite_indices")
        if row.indices_finite and not row.ndvi_ndmi_in_range:reasons.append("ndvi_ndmi_out_of_range")
        if row.scl_cloud_shadow_either_grid:reasons.append("scl_cloud_or_shadow")
        if not row.clear_surface_both_scl_grids and not row.scl_cloud_shadow_either_grid:reasons.append("scl_invalid_unclassified_or_snow")
        if not row.classi_clear:reasons.append("classi_cloud_cirrus_or_snow")
        if row.clear_valid_index_candidate and not row.vegetation_candidate:reasons.append("nonvegetation_scl")
        if not row.time_gap_within_168h:reasons.append("time_gap_exceeds_168h")
        reason.append("accepted_vegetation_candidate" if row.vegetation_candidate else ";".join(reasons))
    p["exclusion_reasons"]=reason
    stages=[];mask=np.ones(len(p),bool)
    for name,check in [("gedi_quality_input",np.ones(len(p),bool)),("all_six_bands_valid",p.all_bands_valid),
            ("finite_indices",p.indices_finite),("ndvi_ndmi_range",p.ndvi_ndmi_in_range),
            ("time_gap_within_168h",p.time_gap_within_168h),("clear_surface_both_scl_grids",p.clear_surface_both_scl_grids),
            ("classi_clear",p.classi_clear),("vegetation_both_grids",p.scl_native_20m.eq(4)&p.scl_aligned_10m.eq(4))]:
        before=int(mask.sum());mask &= np.asarray(check)
        stages.append({"stage":name,"remaining":int(mask.sum()),"excluded_at_stage":before-int(mask.sum())})
    assert np.array_equal(mask,p.vegetation_candidate.to_numpy())
    stage_table=pd.DataFrame(stages)
    stage_table.to_csv(OUT/"filter_counts.csv",index=False)
    pd.DataFrame(bandstats).to_csv(OUT/"band_statistics.csv",index=False)
    pd.DataFrame(qualitystats).to_csv(OUT/"quality_mask_counts.csv",index=False)
    sclcounts=pd.DataFrame([{"scl_class":k,"native_20m_points":int(p.scl_native_20m.eq(k).sum()),
                            "aligned_10m_points":int(p.scl_aligned_10m.eq(k).sum())} for k in range(12)])
    sclcounts.to_csv(OUT/"scl_counts.csv",index=False)
    p.to_csv(OUT/"gedi_sentinel_all_points.csv",index=False)
    p.loc[p.clear_valid_index_candidate].to_csv(OUT/"clear_valid_index_candidates.csv",index=False)
    final=p.loc[p.vegetation_candidate].copy()
    final.to_csv(OUT/"vegetation_candidates.csv",index=False)
    summary={"created_at_utc":datetime.now(timezone.utc).isoformat(),"product":PRODUCT,
        "gedi_input":len(p),"all_six_bands_valid":int(p.all_bands_valid.sum()),
        "clear_valid_index_candidates":int(p.clear_valid_index_candidate.sum()),"vegetation_candidates":len(final),
        "unique_vegetation_pixels":len(final[["raster_row","raster_col"]].drop_duplicates()),
        "scl_native_aligned_disagreements":int(p.scl_boundary_disagreement.sum()),
        "cloud_shadow_points_either_scl_grid":int(p.scl_cloud_shadow_either_grid.sum()),
        "classi_cloud_cirrus_snow_points":int((~p.classi_clear).sum()),
        "forest_confirmed_new_samples":None,"forest_status":"not_assessed; SCL vegetation is not forest",
        "quantification":q,"band_offsets":offsets,"tile_sensing_time_utc_from_xml":tile_time,
        "catalogue_start_utc":catalogue_time,"time_gap_hours_tile_min":p.time_gap_hours_tile_metadata.min(),
        "time_gap_hours_tile_max":p.time_gap_hours_tile_metadata.max(),
        "crop_grid":{"crs":str(ref["crs"]),"width":ref["width"],"height":ref["height"],
                     "transform":list(ref["transform"])[:6],"bounds_utm":list(bounds),"resolution_m":10},
        "quality_policy":"all 8 MSK_QUALIT layers zero and detector ID 1..12; mask native DN before bilinear; full valid interpolation support",
        "scl_policy":"clear surface 4/5/6 on native and aligned grids; final vegetation 4/4; CLASSI cloud/cirrus/snow all zero",
        "cloud_probability_policy":"recorded only; no arbitrary probability threshold applied",
        "evi_policy":"finite required; no clipping; abs>1 and denominator<0.05 flagged only",
        "input_points_sha256":points_hash,"source_zip_md5":digest(ARCHIVE,"md5"),
        "original_baseline_csv_unchanged":baseline_hash==digest(ROOT/"week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv"),
        "input_points_unchanged":points_hash==digest(POINTS),"model_training_performed":False}
    (OUT/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    # Verify emitted rasters at the actual sample cells, not just in-memory arrays.
    for name in ["ndvi","evi","ndmi"]:
        with rasterio.open(OUT/"rasters"/(name+".tif")) as ds:
            sampled=np.asarray([v[0] for v in ds.sample(zip(x,y))],dtype=float)
            sampled[sampled==NODATA]=np.nan
            assert np.allclose(sampled,p[name],equal_nan=True,atol=1e-6)
    assert final[["ndvi","evi","ndmi","agbd_Mg_ha"]].notna().all().all()
    assert summary["original_baseline_csv_unchanged"] and summary["input_points_unchanged"]
    assert set(final.point_id).issubset(set(p.point_id))
    validation={"status":"passed","unique_input_ids":True,"stage_counts_match_final_selection":True,
                "written_index_rasters_match_all_1055_samples":True,"required_final_values_finite":True,
                "source_zip_preserved":True,"previous_data_unchanged":True,"training_performed":False}
    (OUT/"validation.json").write_text(json.dumps(validation,indent=2),encoding="utf-8")
    plot(p,arrays,indices,scl,ref)
    report(summary,stage_table,sclcounts,p)
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)


def plot(p,arrays,indices,scl,ref):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    b=array_bounds(ref["height"],ref["width"],ref["transform"])
    extent=(b[0],b[2],b[1],b[3])
    rgb=np.stack([arrays[x] for x in ["B04","B03","B02"]],axis=-1)
    rgb=np.clip(rgb/0.3,0,1);rgb=np.nan_to_num(rgb,nan=0)
    fig,axes=plt.subplots(1,3,figsize=(17,6),constrained_layout=True)
    axes[0].imshow(rgb,extent=extent,origin="upper");axes[0].set_title("SAFE RGB (display stretch only)")
    im=axes[1].imshow(indices["ndvi"],extent=extent,origin="upper",vmin=-1,vmax=1,cmap="RdYlGn")
    fig.colorbar(im,ax=axes[1],shrink=.7);axes[1].set_title("NDVI from XML-calibrated reflectance")
    axes[2].imshow(scl,extent=extent,origin="upper",cmap="tab20",vmin=0,vmax=11)
    for selected,color,label in [(p.vegetation_candidate,"lime","Vegetation candidate"),(~p.vegetation_candidate,"magenta","Excluded")]:
        axes[2].scatter(p.loc[selected,"easting_m"],p.loc[selected,"northing_m"],s=7,color=color,label=f"{label}: {selected.sum()}")
    axes[2].legend(fontsize=7);axes[2].set_title("SCL and actual GEDI filtering; forest unassessed")
    for ax in axes:ax.set(xlabel="Easting (EPSG:32652, m)",ylabel="Northing (m)")
    fig.savefig(OUT/"quality_overview.png",dpi=160);plt.close(fig)


def report(s,stages,sclcounts,p):
    text=f'''# Week06 Step 5 — 후보 A 영상 전처리 및 GEDI 결합

원본 SAFE ZIP을 보존하고 새 `source/`에 안전하게 압축 해제했다. 다운로드·API 인증·모델 학습은 하지 않았고 Week03~05 및 이전 단계 입력을 수정하지 않았다.

## 결과

- GEDI 품질 통과 입력: **{s['gedi_input']:,}점**
- 여섯 밴드 유효: **{s['all_six_bands_valid']:,}점**
- 구름·그림자·무효 SCL 등을 제외한 지수 후보: **{s['clear_valid_index_candidates']:,}점**
- 두 SCL 격자에서 모두 식생인 후보: **{s['vegetation_candidates']:,}점**, 고유 10m 픽셀 **{s['unique_vegetation_pixels']:,}개**
- **산림으로 확정한 추가 학습점 수는 미확정**이다. SCL=4는 식생이며 산림 경계가 아니다.

## 단계별 필터

{md(stages)}

단계별 수는 누적 조건이다. 개별 제외 사유는 서로 중복될 수 있으므로 합산하지 않는다. 전체 1,055점의 값·품질·탈락 사유는 `gedi_sentinel_all_points.csv`에 보존했다. 식생 후보는 `vegetation_candidates.csv`에 별도 저장했다. 같은 픽셀의 여러 GEDI 샷은 자동 평균하거나 삭제하지 않고 duplicate_pixel로 표시했다.

## 반사율과 격자

제품은 `{PRODUCT}`다. 원본 XML에서 quantification={s['quantification']}, 여섯 밴드 offset=-1000을 읽었으며 **(DN-1000)/10000**으로 보정했다. DN=0 NoData와 65535 포화를 먼저 제외했다. 음수와 1 초과 반사율을 임의로 자르지 않았다. 이는 새 SAFE 원제품의 변환이며 기존 Browser TIFF의 /65535 검증을 완료했다는 뜻은 아니다.

B04 원본 10m 격자를 기준으로 GEDI 좌표 범위에 사방 1km 여유를 둔 영역만 처리했다. CRS는 EPSG:32652, 크기는 {s['crop_grid']['width']}×{s['crop_grid']['height']} 픽셀이다. 타일 밖 요청 영역은 타일 경계에서 자른다. 전체 110km 타일을 메모리에 적재하지 않는다.

B11/B12는 원본 20m에서 반사율을 보정하고 품질 불량 픽셀을 제거한 뒤 bilinear로 10m에 정렬했다. 보간에 기여한 원본 픽셀이 모두 유효한 경우만 인정했다. SCL·품질·검출기 마스크는 최근린으로 정렬했다. 20m 밴드를 10m로 바꿔도 실제 관측 해상도가 높아지는 것은 아니다.

기존 Week03 `vegetation_indices.py`의 순수 계산 함수를 재사용하여 NDVI·EVI·NDMI를 계산했다. 모든 여섯 밴드가 유효한 영역을 분석 마스크로 사용했다. NDVI·NDMI 범위 이탈은 후보에서 제외하고, EVI는 유한값만 요구하며 범위 밖 값과 작은 분모는 진단 플래그로 남겼다. 원본 계산 결과를 클리핑하지 않았다.

## 품질 마스크와 SCL

여섯 밴드 각각 MSK_QUALIT의 8개 층을 검사했다: ancillary lost/degraded, MSI lost/degraded, defective, NoData, partially corrected crosstalk, saturated. 하나라도 1이면 제외하는 보수적 기준을 적용했다. MSK_DETFOO는 검출기 ID 1~12만 허용했다. 정의는 [공식 S2 Processing 마스크 설명](https://sentiwiki.copernicus.eu/web/s2-processing)에 근거한다.

SCL 원본 20m와 정렬된 10m를 모두 샘플링했다. 두 격자 모두 4·5·6인 점만 명확한 표면 후보로 보고 2(그늘), 3(구름 그림자), 8·9·10(구름), 0·1·7·11을 제외했다. 최종 식생 후보는 두 격자 모두 4인 점이다. 원본/정렬 클래스 불일치 {s['scl_native_aligned_disagreements']}점이며, 어느 격자에서든 구름·그림자인 점은 {s['cloud_shadow_points_either_scl_grid']}점이다.

범주형 마스크의 재투영 NoData는 255로 예약해 정상 플래그 0을 보존했다. 정상 0 플래그 및 bilinear 유효성 지원 영역 보존 회귀 검사 두 개를 통과했다.

MSK_CLASSI의 opaque·cirrus·snow 층도 확인하여 세 층 모두 0인 경우만 통과시켰다. 그 외 점은 {s['classi_cloud_cirrus_snow_points']}점이다. CLDPRB·SNWPRB는 수치만 기록하고 임의 확률 임계값을 추가하지 않았다. 모든 밴드·SCL·마스크는 동일한 제품에서 가져왔다.

{md(sclcounts)}

## 관측시각

카탈로그 시작 시각은 {s['catalogue_start_utc']}이고 **실제 MTD_TL.xml의 SENSING_TIME은 {s['tile_sensing_time_utc_from_xml']}**이다. 서로 다른 메타데이터 필드이므로 함께 보존했다. 타일 메타데이터 시각 기준 GEDI와의 경과시간은 {s['time_gap_hours_tile_min']:.4f}~{s['time_gap_hours_tile_max']:.4f}시간이다. 이는 개별 GEDI 위치의 픽셀별 정밀 촬영시각을 독립 검증했다는 의미는 아니다.

이번 데이터에 실제 경과시간 168시간 이내 조건을 적용했다. 관측점은 중심 위치가 속한 단일 지수 픽셀에서 추출했으며 GEDI footprint 면적 평균은 계산하지 않았다.

## 산출물 및 검증

| 파일 | 내용 |
| --- | --- |
| summary.json / validation.json | 품질 집계·처리 설정·입력 보존 및 검증 |
| gedi_sentinel_all_points.csv | 1,055점 전부의 반사율·지수·SCL·품질·제외 사유 |
| clear_valid_index_candidates.csv | 구름 등 제외 후 지수 후보 |
| vegetation_candidates.csv | SCL 식생 후보; 산림 학습점으로 확정하지 않음 |
| filter_counts.csv / scl_counts.csv | 누적 필터 및 SCL 분포 |
| band_statistics.csv / quality_mask_counts.csv | 반사율 범위와 밴드별 품질 플래그 수 |
| rasters/ | 여섯 반사율·세 지수·SCL GeoTIFF, 원본 좌표 격자 기반 |
| quality_overview.png | RGB·NDVI·SCL과 GEDI 판정 지도 |

저장된 지수 GeoTIFF를 다시 열어 1,055점 모두에서 CSV 값과 일치함을 확인했다. 식생 후보의 필수 지수·AGBD 결측이 없고, 누적 필터와 최종 선택 수가 일치하며 기존 90점 CSV·입력 CSV·원본 ZIP을 보존했다.

![품질 확인](quality_overview.png)

## 다음 단계

식생 후보 지도와 분포를 검토한 뒤 산림 경계·토지피복 자료로 산림 표본을 구분한다. 그 후 기존 90점과 지역·궤도·날짜·반사율 계보를 맞춰 통합하고 공간 버퍼와 미사용 지역 홀드아웃을 설계한다. 새 SAFE와 기존 Analytical TIFF는 반사율 근거가 다르므로 확인 없이 합쳐 학습하지 않는다. 모델 학습은 이번 범위에 포함하지 않았다.
'''
    (OUT/"report.md").write_text(text,encoding="utf-8")


if __name__=="__main__":
    run()
