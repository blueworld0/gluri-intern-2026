"""Inventory local project geodata without downloads or sample filtering."""
from pathlib import Path
from collections import Counter
import csv
import json
import os
import re
import zipfile
from datetime import datetime, timezone

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/"data_inventory"
SKIP={".git",".venv","venv","env","__pycache__",".matplotlib"}
GEO={".shp",".shx",".dbf",".prj",".gpkg",".geojson",".jsonl",".gdb",".kml",".kmz",".gml",".fgb",".tif",".tiff",".jp2",".img",".asc",".vrt",".nc",".h5",".hdf5",".zip",".7z",".rar"}


def run():
    rows=[];archives=[];counts=Counter();filegdb=[];unknown=[]
    for folder,dirs,files in os.walk(ROOT):
        dirs[:]=[d for d in dirs if d not in SKIP and (Path(folder)/d).resolve()!=OUT.resolve()]
        for d in dirs:
            if d.lower().endswith(".gdb"):filegdb.append((Path(folder)/d).relative_to(ROOT).as_posix())
        for name in files:
            path=Path(folder)/name;ext=path.suffix.lower()
            if ext not in GEO:continue
            rel=path.relative_to(ROOT).as_posix();counts[ext]+=1
            if ext==".geojson":
                data=json.loads(path.read_text(encoding="utf-8"))
                props=data["features"][0].get("properties",{}) if data.get("features") else {}
                role="arbitrary_candidate_grid" if "boundary_type" in props else "satellite_product_footprint" if "product_uuid" in props else "needs_review"
            elif ext in {".h5",".hdf5"}:role="GEDI_L4A_source"
            elif ext==".zip":
                role="Sentinel2_SAFE_archive" if ".SAFE.zip" in name else "needs_review"
                with zipfile.ZipFile(path) as z:
                    for item in z.infolist():
                        if not item.is_dir():archives.append({"archive":rel,"member":item.filename,"uncompressed_bytes":item.file_size})
            elif ext in {".tif",".tiff",".jp2"}:
                if re.search(r"MSK_|SCL",name,re.I):role="Sentinel2_scene_or_quality_mask"
                elif re.search(r"ndvi|ndmi|evi",name,re.I):role="derived_vegetation_index"
                elif re.search(r"B(?:0[1-9]|1[012]|8A)|AOT|WVP|TCI|PVI",name,re.I):role="Sentinel2_band_or_auxiliary"
                else:role="needs_review"
            else:role="needs_review"
            rows.append({"path":rel,"extension":ext,"bytes":path.stat().st_size,"role":role})
            if role=="needs_review":unknown.append(rel)
    OUT.mkdir(parents=True,exist_ok=True)
    for name,data,columns in [("geodata_inventory.csv",rows,["path","extension","bytes","role"]),
                              ("archive_members.csv",archives,["archive","member","uncompressed_bytes"])]:
        with (OUT/name).open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(data)
    with (ROOT/"week06/step5/outputs/vegetation_candidates.csv").open(encoding="utf-8",newline="") as f:p=list(csv.DictReader(f))
    bounds={"west":min(float(x["longitude"]) for x in p),"south":min(float(x["latitude"]) for x in p),
            "east":max(float(x["longitude"]) for x in p),"north":max(float(x["latitude"]) for x in p)}
    summary={"checked_at_utc":datetime.now(timezone.utc).isoformat(),"scope":"workspace only; hidden and Git-ignored source data included",
        "excluded_directories":sorted(SKIP),"geospatial_file_counts_by_extension":dict(sorted(counts.items())),
        "file_geodatabase_directories":filegdb,"unclassified_geodata_files":unknown,
        "forest_boundary_dataset_found":False,"landcover_dataset_found":False,"forest_inventory_map_found":False,
        "existing_vector_features":{"candidate_grids":5,"sentinel2_product_footprints":16},
        "vegetation_points":len(p),"bounds_wgs84":bounds,"forest_status_counts":dict(Counter(x["forest_status"] for x in p)),
        "new_downloads":False,"external_requests":False,"samples_filtered":False,
        "note":"No dataset is inferred from a document citation or SCL vegetation class."}
    assert len(p)==441 and summary["forest_status_counts"]=={"not_assessed":441}
    assert not unknown and not filegdb
    assert counts[".geojson"]==2
    (OUT/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    report=f'''# Week06 Step 6 — 로컬 산림·토지피복 자료 조사

## 결론

**현재 저장소에는 후보 A의 441점을 산림·비산림으로 판정할 산림 경계, 임상도 또는 독립적인 토지피복 원자료가 없다.** 파일명뿐 아니라 벡터 속성, 압축 파일 목록, 기존 연구 계획과 처리 코드를 대조했다. 이번에는 다운로드·외부 검색·표본 선별·모델 학습을 수행하지 않았다.

## 조사 범위와 실제 파일

프로젝트 전체를 조사했으며 Git에서 제외된 HDF5, TIFF, JP2와 다운로드 ZIP, Step5의 압축 해제 SAFE도 포함했다. .git·가상환경·캐시는 제외했다. 파일별 경로와 역할은 `geodata_inventory.csv`, ZIP 내부 목록은 `archive_members.csv`에 저장했다.

| 자료 | 확인 결과 | 산림 판정에 사용 가능 여부 |
| --- | --- | --- |
| SHP·GeoPackage·FileGDB·KML 등 산림 벡터 | 확인되지 않음 | 적용할 원자료 없음 |
| 토지피복·임상도 래스터 | 확인되지 않음 | 분류 코드·작성 연도·좌표계 확인 불가 |
| Step2 recommended_candidate_grids.geojson | 임의 격자 5개, boundary_type 속성 확인 | 행정·산림 경계가 아님 |
| Step3 product_footprints.geojson | 위성 제품 footprint 16개, 제품 UUID 속성 확인 | 촬영 범위이며 산림 경계가 아님 |
| Sentinel-2 SCL·품질 마스크 | 원본 SAFE 및 기존 영상에 존재 | 구름·유효성·식생 구분에 활용 가능하나 산림 확정 불가 |
| NDVI·EVI·NDMI | 기존 및 Step5 결과에 존재 | 식생활력 지수이며 산림 경계·토지피복 분류자료가 아님 |
| GEDI L4A와 식생 후보 CSV | AGBD·좌표·품질·식생 후보 보유 | AGBD가 있다고 해당 지점을 산림으로 확정할 수 없음 |

검토가 필요한 미분류 공간자료 파일은 {len(unknown)}개다. 현재 확보한 모든 공간자료는 GEDI, Sentinel-2, 식생지수, 임의 후보 격자, 위성 footprint로 분류됐다.

## 기존 문서의 근거

`week04/research/paper1_reproduction_plan.md`의 현재 자료 표는 DEM·토지피복도와 NFI 기반 임상도를 각각 **없음**으로 기록한다. EGIS 중분류 토지피복과 코드표, 임상도 속성·갱신 연도를 후속 확보 대상으로 적었다. 이는 참고 자료 이름이지 실제 파일이 확보됐다는 뜻이 아니다.

`week05/gedi_l2a_feasibility.md`도 토지피복·임상도 등 보조 입력이 필요하다고 기록한다. Week03 SCL 코드와 Week06 Step5 코드는 식생 후보를 선별하지만 독립 산림 레이어를 불러오거나 산림 코드로 공간 결합하는 단계는 없다.

## 다음 자료가 덮어야 할 범위

Step5 후보 **441점**의 실제 WGS84 좌표 범위:

- 서쪽: {bounds['west']:.9f}°E
- 동쪽: {bounds['east']:.9f}°E
- 남쪽: {bounds['south']:.9f}°N
- 북쪽: {bounds['north']:.9f}°N

441점 모두 forest_status=not_assessed다. 점 중심만 판정한다면 위 범위 전체의 피복이 필요하고, GEDI footprint 면적과 경계 혼합을 검사하려면 범위 밖으로 여유를 둔 자료도 필요하다. 기존 Sentinel-2 지수 격자는 EPSG:32652, 10m이며 새 자료는 원래 CRS를 확인한 뒤 변환해야 한다. 범주형 래스터는 최근린 정렬을 사용하되 원자료 해상도와 작성 시점을 기록해야 한다.

## 다음 작업 제안

1. 기존 프로젝트가 언급한 **EGIS 토지피복도 또는 산림 임상도**의 공식 제공 경로와 대상 지역 자료 확보 가능성을 별도로 조사한다. 이번에는 외부 사이트의 최신 접근·로그인·다운로드 조건을 확인하지 않았다.
2. 실제 자료를 받기 전에 작성 연도·공간 범위·CRS·원본 해상도 또는 축척·분류 코드표·이용 조건을 확인한다. 2025년 관측과의 시차도 기록한다.
3. 산림으로 인정할 코드와 농경지·초지·도시 수목·공원 등의 처리 방식을 먼저 정한다. 토지피복상의 산림과 법적 산림 경계, 수목 피복은 서로 같은 정의가 아니다.
4. 확보한 자료를 보존하고 441점을 산림 후보·비산림·경계 또는 판정 불가로 나눈다. 경계에 걸리는 GEDI footprint는 별도 표시한다.
5. 기존 90점에도 같은 기준을 적용한 뒤 비교한다. 식생 후보 441점을 바로 산림 학습점으로 이름만 바꾸거나 기존 90점과 자동 합산하지 않는다.

현재는 새 레이어가 없으므로 **산림 표본 수를 확정할 수 없다.** 다음으로 필요한 일은 공식 자료의 접근·확보 방법 조사다.
'''
    (OUT/"report.md").write_text(report,encoding="utf-8")
    assert "???" not in report and "\ufffd" not in report
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=="__main__":run()
