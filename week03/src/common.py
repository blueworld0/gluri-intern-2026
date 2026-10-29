"""여러 단계에서 같이 쓰는 함수들."""
import os
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds

# 클라우드 저장소의 GeoTIFF를 통째로 받지 않고 필요한 부분만 읽기 위한 설정
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif,.TIF")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "5")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")
# (회사/학교 프록시 환경이면 HTTPS_PROXY 를 GDAL 에도 넘겨 줌)
if os.environ.get("HTTPS_PROXY") and not os.environ.get("GDAL_HTTP_PROXY"):
    os.environ["GDAL_HTTP_PROXY"] = os.environ["HTTPS_PROXY"].split("://")[-1]
if os.path.exists("/root/.ccr/ca-bundle.crt"):
    os.environ.setdefault("CURL_CA_BUNDLE", "/root/.ccr/ca-bundle.crt")


def clip_remote(urls, out_path, bbox):
    """원격 GeoTIFF(들)에서 bbox(경위도) 부분만 잘라 out_path 에 저장.

    urls 에 여러 타일을 넣으면 이어 붙여서 자릅니다 (bbox 가 타일 경계에 걸칠 때).
    """
    from rasterio.merge import merge
    if isinstance(urls, str):
        urls = [urls]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        print("  이미 있음:", out_path.name)
        return out_path
    srcs = [rasterio.open("/vsicurl/" + u) for u in urls]
    b = transform_bounds("EPSG:4326", srcs[0].crs, *bbox, densify_pts=21)
    arr, tr = merge(srcs, bounds=b)   # 타일 밖은 nodata(없으면 0)로 채움
    prof = srcs[0].profile.copy()
    prof.update(driver="GTiff", height=arr.shape[1], width=arr.shape[2],
                transform=tr, compress="deflate")
    for k in ("blockxsize", "blockysize", "tiled", "photometric", "interleave"):
        prof.pop(k, None)
    with rasterio.open(out_path, "w", **prof) as dst:
        dst.write(arr)
    for s_ in srcs:
        s_.close()
    print("  저장:", out_path.name, arr.shape)
    return out_path
