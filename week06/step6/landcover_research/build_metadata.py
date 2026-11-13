"""Write research metadata and local AOI only; no network or data filtering."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import pandas as pd
from rasterio.warp import transform, transform_bounds

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def dump(name,value):
    (HERE/name).write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding="utf-8")


def run():
    source=ROOT/"week06/step5/outputs/vegetation_candidates.csv"
    before=hashlib.sha256(source.read_bytes()).hexdigest()
    p=pd.read_csv(source,dtype={"shot_number":"string"})
    assert len(p)==p.point_id.nunique()==441
    assert p.forest_status.eq("not_assessed").all()
    bbox=[float(p.longitude.min()),float(p.latitude.min()),float(p.longitude.max()),float(p.latitude.max())]
    x,y=transform("EPSG:4326","EPSG:32652",p.longitude.tolist(),p.latitude.tolist())
    padded=transform_bounds("EPSG:32652","EPSG:4326",min(x)-100,min(y)-100,max(x)+100,max(y)+100,densify_pts=21)
    features=[]
    for label,b in [("actual_441_point_bbox",bbox),("request_helper_bbox_100m_margin_utm",padded)]:
        w,s,e,n=b
        features.append({"type":"Feature","properties":{"role":label,"source_points":441,
            "boundary_type":"arbitrary_request_helper_not_admin_or_forest"},
            "geometry":{"type":"Polygon","coordinates":[[[w,s],[e,s],[e,n],[w,n],[w,s]]]}})
    dump("target_aoi.geojson",{"type":"FeatureCollection","features":features})
    sources=[
        {"id":"egis_land","title":"EGIS 토지피복지도 소개·분류체계·제작현황","url":"https://egis.me.go.kr/intro/land.do",
         "supports":["2025년 중분류·세분류 전국 제공 안내","310/320/330 및 311/321/331 코드·정의","2025년 세분류의 2024년 기본영상","SHP 및 중분류 GeoTIFF 등 형식"],
         "access":"official indexed page text; direct open timed out"},
        {"id":"egis_request","title":"EGIS 자료신청 안내","url":"https://egis.me.go.kr/req/intro.do",
         "supports":["로그인·신청 후 내역에서 다운로드","메타/지도·사각형·다각형 도엽 검색","공개제한 자료의 별도 절차"],
         "access":"official indexed page text; direct open timed out"},
        {"id":"egis_wms","title":"EGIS 토지피복지도 맵 서비스","url":"https://egis.me.go.kr/api/land.do",
         "supports":["lv3_2025y 및 lv2_2025y","WMS 지도 예제의 EPSG:3857"],"access":"official indexed page text; no WMS request sent"},
        {"id":"egis_map","title":"EGIS 통합지도","url":"https://egis.me.go.kr/map/map.do",
         "supports":["현재 지도 안내의 2024년 영상 기준 현행화 설명"],"access":"official indexed page text; no UI interaction"},
        {"id":"egis_statistics","title":"EGIS 국가 토지피복 통계 안내","url":"https://egis.me.go.kr/intro/stats.do",
         "supports":["지도 제작연도와 통계 기준연도 구분","임상도 등 참조자료 활용"],"access":"official indexed page text"},
        {"id":"forest_service","title":"산림공간정보서비스","url":"https://map.forest.go.kr/forest/",
         "supports":["2025 대축척 임상도 신청 항목","2024 임상도(2020~2024) 구분","유통신청 본인인증","일반인 원본 TIFF 안내","임상·수종·갱신년도 속성 안내"],
         "access":"official indexed page text; direct root yielded no readable text; no identity verification"},
        {"id":"forest_open_data","title":"산림청 공공데이터 임상도 안내","url":"https://www.forest.go.kr/kfsweb/kfi/kfs/trail/fTMap.do?mn=NKFS_06_08_02&pblicDataId=PBD0000059&subTitle=%EC%9E%84%EC%83%81%EB%8F%84&tabs=1",
         "supports":["무료신청·일반인/공공연구 경로","과거 제5차 임상도 SHP 제공 사례"],
         "access":"official indexed page text; direct open returned HTTP 403"},
        {"id":"forest_legacy_schema","title":"산림청 과거 임상도 활용정보","url":"https://www.forest.go.kr/kfsweb/opda/dataMng/pblicDataPrusContView.do?pblicDataId=PBD0000059&mn=NKFS_01_01",
         "supports":["과거 imsang_5th Polygon SHP ZIP·GRS80 TM 중부 설명; 2025에 그대로 적용 불가"],
         "access":"official indexed page text observed on session-bearing legacy URL; canonical non-session reference stored"},
    ]
    dump("sources.json",{"checked_at_utc":datetime.now(timezone.utc).isoformat(),"sources":sources,
                          "datasets_downloaded":False,"login_or_application_performed":False})
    code_url=sources[0]["url"]
    codes=[{"class_name":name,"medium_code":m,"fine_code":f,"definition_summary":definition,
            "source_url":code_url,"status":"official_classification_verified_not_file_schema"}
        for name,m,f,definition in [("활엽수림",310,311,"활엽수림 비중이 전체 임분 면적의 75% 이상"),
                                    ("침엽수림",320,321,"침엽수림 비중이 전체 임분 면적의 75% 이상"),
                                    ("혼효림",330,331,"활엽수·침엽수 임분이 혼재하며 각각 75% 이하")]]
    with (HERE/"forest_codes.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(codes[0]));w.writeheader();w.writerows(codes)
    access=[
        {"dataset":"EGIS_2025_fine","official_route":"https://egis.me.go.kr/req/intro.do","access":"login_then_application_then_download","analysis_format":"SHP with code attributes","native_epsg":"not_confirmed","status":"official_documentation_only"},
        {"dataset":"EGIS_2025_medium","official_route":"https://egis.me.go.kr/req/intro.do","access":"login_then_application_then_download","analysis_format":"classification GeoTIFF or SHP","native_epsg":"not_confirmed","status":"official_documentation_only"},
        {"dataset":"EGIS_WMS_2025","official_route":"https://egis.me.go.kr/api/land.do","access":"map_service_docs_only","analysis_format":"rendered_map_not_source_classification","native_epsg":"3857 example is not source CRS","status":"no_service_request_sent"},
        {"dataset":"Forest_2025_large_scale","official_route":"https://map.forest.go.kr/forest/","access":"identity_verification_and_application","analysis_format":"TIFF listed for public access; vector/attribute availability needs confirmation","native_epsg":"not_confirmed","status":"no_login_or_application"},
    ]
    with (HERE/"access_comparison.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(access[0]));w.writeheader();w.writerows(access)
    result={"gedi_vegetation_points":441,"input_sha256":before,"input_path":str(source.relative_to(ROOT)),
        "actual_bbox_wgs84":bbox,"request_helper_bbox_100m_margin_utm_wgs84":list(padded),
        "gedi_observation_date_utc":"2025-04-25","recommended_dataset":"EGIS_2025_fine_SHP",
        "egis_fine_map_production_year":2025,"egis_fine_base_imagery_year":2024,
        "medium_forest_codes":[310,320,330],"fine_forest_codes":[311,321,331],
        "egis_coverage":"national coverage announced; actual local tile supply and point inclusion not tested",
        "forest_2025_tile_coverage":"not_tested","native_epsg_2025_egis":None,"native_epsg_2025_forest":None,
        "actual_tile_ids":None,"file_schema":None,"forest_classified_sample_count":None,
        "input_unchanged":before==hashlib.sha256(source.read_bytes()).hexdigest(),
        "datasets_downloaded":False,"login_performed":False,"application_submitted":False,"model_trained":False}
    dump("research_summary.json",result)
    report=(HERE/"report.md").read_text(encoding="utf-8")
    assert "\ufffd" not in report and "???" not in report
    assert result["input_unchanged"]
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":run()
