"""Read-only verification of the user-downloaded SAFE ZIP; no extraction/network."""
from pathlib import Path, PurePosixPath
from datetime import datetime, timezone
import csv
import hashlib
import json
import uuid
import zipfile
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def run():
    config = json.loads((HERE / "product.json").read_text(encoding="utf-8"))
    name = config["product_name"]
    archive = HERE / "downloads" / (name + ".zip")
    metadata_path = HERE / "checks/20261113T045634Z_4eb39837/product_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected_md5 = next(v["Value"] for v in metadata["Checksum"] if v["Algorithm"] == "MD5")
    checksum = hashlib.md5()
    with archive.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    print("Archive MD5 scan completed.", flush=True)
    with (ROOT / "week06/step3/download_manifest_primary.csv").open(encoding="utf-8", newline="") as f:
        required = list(csv.DictReader(f))
    with zipfile.ZipFile(archive) as z:
        entries = z.infolist()
        names = [e.filename for e in entries]
        unsafe = [s for s in names if PurePosixPath(s).is_absolute() or ".." in PurePosixPath(s).parts or "\\" in s or ":" in s]
        missing = [r["relative_path"] for r in required if r["relative_path"] not in names]
        wrong_sizes = [r["relative_path"] for r in required if r["relative_path"] in names
                       and z.getinfo(r["relative_path"]).file_size != int(r["bytes"])]
        bad_crc = z.testzip()
        print("ZIP CRC scan completed.", flush=True)
        product_xml = z.read(name + "/MTD_MSIL2A.xml")
        tree = ET.fromstring(product_xml)
        def local(tag):
            return tag.rsplit("}", 1)[-1]
        def values(tag):
            return [e.text for e in tree.iter() if local(e.tag) == tag]
        offsets = {e.attrib.get("band_id"): e.text for e in tree.iter() if local(e.tag) == "BOA_ADD_OFFSET"}
        bands = {e.attrib.get("bandId"): e.attrib.get("physicalBand") for e in tree.iter() if local(e.tag) == "Spectral_Information"}
        spectral = {bands.get(k, k): float(v) for k, v in offsets.items()}
        result = {
            "checked_at_utc": datetime.now(timezone.utc).isoformat(),
            "archive": str(archive.resolve()), "product_name": name,
            "archive_bytes": archive.stat().st_size, "catalogue_expected_bytes": metadata["ContentLength"],
            "size_matches": archive.stat().st_size == metadata["ContentLength"],
            "md5": checksum.hexdigest(), "catalogue_expected_md5": expected_md5,
            "md5_matches": checksum.hexdigest() == expected_md5,
            "zip_crc_all_members_pass": bad_crc is None, "first_crc_failure": bad_crc,
            "zip_entries": len(entries), "uncompressed_bytes": sum(e.file_size for e in entries),
            "unsafe_paths": unsafe, "duplicate_member_names": len(names) != len(set(names)),
            "only_expected_safe_root": all(PurePosixPath(s).parts[0] == name for s in names),
            "required_file_count": len(required), "missing_required_files": missing,
            "required_file_size_mismatches": wrong_sizes,
            "product_uri_from_xml": values("PRODUCT_URI"),
            "processing_baseline_from_xml": values("PROCESSING_BASELINE"),
            "boa_quantification_value_from_xml": values("BOA_QUANTIFICATION_VALUE"),
            "boa_add_offsets_by_band_from_xml": spectral,
            "cloud_coverage_assessment_from_xml": values("Cloud_Coverage_Assessment"),
            "archive_extracted": False, "image_pixels_processed": False,
        }
    passed = (result["size_matches"] and result["md5_matches"] and result["zip_crc_all_members_pass"]
              and not unsafe and not result["duplicate_member_names"] and result["only_expected_safe_root"]
              and not missing and not wrong_sizes and name in result["product_uri_from_xml"])
    result["status"] = "passed" if passed else "needs_review"
    out = HERE / "checks" / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_download_" + uuid.uuid4().hex[:8])
    out.mkdir(parents=True, exist_ok=False)
    with (out / "download_verification.json").open("x", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    print("Verification report:", out / "download_verification.json", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(run())
