"""Read-only public catalogue metadata search; never accesses imagery or tokens."""
from pathlib import Path
from datetime import datetime, timezone
import json
import urllib.request
import urllib.parse
import urllib.error
import hashlib
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    source = ROOT / "week06/step2/recommended_candidate_points.csv"
    data = pd.read_csv(source, dtype={"shot_number": "string"})
    points = data.loc[data.region_id.eq("G025_515_140") & data.observation_date_utc.eq("2025-04-25")].copy()
    assert len(points) == 1055 and points.point_id.nunique() == 1055
    assert points.longitude.between(128.75, 129, inclusive="left").all()
    assert points.latitude.between(35, 35.25, inclusive="left").all()
    points.to_csv(HERE / "candidate_a_20250425_gedi_points.csv", index=False)
    summary = {"count": len(points), "bounds_wgs84": {
        "west": points.longitude.min(), "south": points.latitude.min(),
        "east": points.longitude.max(), "north": points.latitude.max()},
        "observed_at_utc_min": points.observed_at_utc.min(),
        "observed_at_utc_max": points.observed_at_utc.max(),
        "beam_counts": points.beam.value_counts().sort_index().to_dict(),
        "source_file": str(source.relative_to(ROOT)), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "quality_status": "GEDI quality only; Sentinel-2 and forest suitability not evaluated"}
    (HERE / "gedi_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    polygon = "POLYGON((128.75 35.00,129.00 35.00,129.00 35.25,128.75 35.25,128.75 35.00))"
    params = {
        "$filter": "Collection/Name eq 'SENTINEL-2' and contains(Name,'MSIL2A') "
            "and ContentDate/Start ge 2025-04-18T00:00:00.000Z "
            "and ContentDate/Start lt 2025-05-03T00:00:00.000Z "
            f"and OData.CSC.Intersects(area=geography'SRID=4326;{polygon}')",
        "$expand": "Attributes", "$top": "1000", "$orderby": "ContentDate/Start asc", "$count": "true"}
    url = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?" + urllib.parse.urlencode(params)
    manifest = {"retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "authentication": "none; no credentials or API keys read or sent", "search_bbox": [128.75,35,129,35.25],
        "time_interval_utc": "[2025-04-18T00:00:00Z,2025-05-03T00:00:00Z)",
        "cloud_filter": None, "odata_parameters": params, "requests": []}
    pages = []
    for page in range(10):
        if not url:
            break
        assert urllib.parse.urlparse(url).hostname == "catalogue.dataspace.copernicus.eu"
        log = {"url": url, "method": "GET"}
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "GLURI-Week06-metadata-research/1.0"})
            with urllib.request.urlopen(request, timeout=45) as response:
                log["status"] = response.status
                payload = response.read(20_000_001)
                if len(payload) > 20_000_000:
                    raise ValueError("Metadata size limit exceeded")
                result = json.loads(payload)
            name = f"catalogue_odata_page_{page+1:02d}.json"
            (HERE / name).write_bytes(payload)
            log.update(response_file=name, records=len(result.get("value", [])))
            pages.extend(result.get("value", []))
            url = result.get("@odata.nextLink")
        except Exception as exc:
            log.update(error_type=type(exc).__name__, error=str(exc))
            if isinstance(exc, urllib.error.HTTPError):
                log["status"] = exc.code
            manifest["requests"].append(log)
            break
        manifest["requests"].append(log)
    manifest.update(records=len(pages), pagination_complete=not bool(url))
    (HERE / "catalogue_query.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"gedi": summary, "catalogue": manifest}, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
