"""List public SAFE directory metadata only. Never call /$value or read imagery."""
from pathlib import Path
from datetime import datetime, timezone
import json
import urllib.request
import urllib.parse
import urllib.error

HERE = Path(__file__).resolve().parent


def run():
    products = json.loads((HERE / "catalogue_odata_page_01.json").read_text())["value"]
    chosen = [v for v in products if ("20250428T020711" in v["Name"] and "T52SDD" in v["Name"])]
    logs, files = [], []
    def listing(url, product, relative):
        parsed = urllib.parse.urlparse(url)
        assert parsed.hostname == "download.dataspace.copernicus.eu"
        assert "$value" not in url and url.endswith("/Nodes")
        entry = {"url": url, "product_uuid": product["Id"], "relative_directory": relative}
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=35) as response:
                entry["status"] = response.status
                payload = response.read(5_000_001)
                assert len(payload) <= 5_000_000
                data = json.loads(payload)
            out = f"nodes_response_{len(logs)+1:02d}.json"
            (HERE / out).write_bytes(payload)
            entry["response_file"] = out
            nodes = data.get("result", data.get("value", []))
        except Exception as exc:
            entry.update(error_type=type(exc).__name__, error=str(exc))
            if isinstance(exc, urllib.error.HTTPError): entry["status"] = exc.code
            logs.append(entry)
            return
        logs.append(entry)
        for node in nodes:
            name = node["Name"]
            path = f"{relative}/{name}" if relative else name
            if node.get("ChildrenNumber", 0) == 0:
                files.append({"product_uuid": product["Id"], "product_name": product["Name"],
                              "relative_path": path, "bytes": node.get("ContentLength"),
                              "evidence": out})
                continue
            # Restrict directory traversal to band, quality, and provenance data.
            if (name.endswith(".SAFE") or name in {"GRANULE", "IMG_DATA", "R10m", "R20m", "QI_DATA", "DATASTRIP"}
                    or (relative.endswith("/GRANULE") or relative.endswith("/DATASTRIP"))):
                child = node.get("Nodes", {}).get("uri")
                if not child:
                    child = url + "(" + urllib.parse.quote(name, safe="._-") + ")/Nodes"
                listing(child, product, path)
    for product in chosen:
        listing(f"https://download.dataspace.copernicus.eu/odata/v1/Products({product['Id']})/Nodes", product, "")
    (HERE / "nodes_query.json").write_text(json.dumps({"retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "authentication": "none", "image_or_xml_content_fetched": False, "requests": logs,
        "listed_files": len(files)}, indent=2), encoding="utf-8")
    (HERE / "listed_product_files.json").write_text(json.dumps(files, indent=2), encoding="utf-8")
    print(json.dumps({"requests": len(logs), "files": len(files), "errors": [x for x in logs if x.get("status") != 200]}, indent=2), flush=True)


if __name__ == "__main__":
    run()
