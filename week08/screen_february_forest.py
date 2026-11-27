"""Screen the February holdout candidates against the week07 EGIS tiles."""

import argparse
import json
from pathlib import Path
import struct
import sys

import numpy as np
import pandas as pd
from rasterio.crs import CRS
from rasterio.warp import transform
from shapely.geometry import Point, box
from shapely.strtree import STRtree

ROOT = Path(__file__).resolve().parents[1]
WEEK07 = ROOT / "week07" / "landcover"
sys.path.insert(0, str(WEEK07))
from classify_points import FOREST, LABELS, polygon_from_record, read_dbf, sha


def load_polygons():
    provenance = json.loads((WEEK07 / "provenance.json").read_text(encoding="utf-8"))
    polygons, attributes, invalid = [], [], []
    for tile in provenance["tiles"]:
        for file in tile["files"]:
            if sha(ROOT / file["path"]) != file["sha256"]:
                raise ValueError(f"Week07 tile checksum mismatch: {file['path']}")
        tile_id = tile["tile"]
        stem = WEEK07 / "source" / tile_id / f"lv2_{tile_id}"
        crs = CRS.from_wkt(stem.with_suffix(".prj").read_text(encoding="utf-8"))
        if crs.to_epsg() != 5186:
            raise ValueError(f"Unexpected CRS for tile {tile_id}: {crs}")
        shp = stem.with_suffix(".shp").read_bytes()
        offset = 100
        for index, row in enumerate(read_dbf(stem.with_suffix(".dbf")), start=1):
            record_id, words = struct.unpack_from(">II", shp, offset)
            end = offset + 8 + 2 * words
            record = shp[offset + 8:end]
            offset = end
            code = row["L2_CODE"]
            if record_id != index or LABELS.get(code) != row["L2_NAME"]:
                raise ValueError(f"Unexpected record or code in tile {tile_id}: {index}")
            geom, _ = polygon_from_record(record)
            if geom is None or geom.is_empty or not geom.is_valid:
                invalid.append(box(*struct.unpack_from("<dddd", record, 4)))
                continue
            polygons.append(geom)
            attributes.append((code, row["L2_NAME"], tile_id, record_id))
        if offset != len(shp):
            raise ValueError(f"Record length mismatch in tile {tile_id}")
    return STRtree(polygons), attributes, STRtree(invalid)


def classify_points(data, tree, attributes, invalid_tree):
    x, y = transform("EPSG:4326", "EPSG:5186", data.longitude.tolist(), data.latitude.tolist())
    results = []
    for east, north in zip(x, y):
        point = Point(east, north)
        hits = tree.query(point, predicate="intersects")
        distances = [point.distance(tree.geometries[i].boundary) for i in hits]
        selected = (None, None, None, None)
        status = "unclassified"
        if len(invalid_tree.query(point, predicate="intersects")):
            reason = "invalid_geometry_bbox_at_point"
        elif len(hits) > 1:
            reason = "multiple_polygon_matches"
        elif distances and distances[0] <= 0.001:
            reason = "within_1mm_polygon_boundary"
        elif len(hits) == 1:
            selected = attributes[hits[0]]
            status = "forest" if selected[0] in FOREST else "nonforest"
            reason = "single_valid_polygon_interior"
        else:
            reason = "no_polygon_match"
        results.append((status, reason, len(hits), *selected, east, north))
    return pd.DataFrame(
        results,
        index=data.index,
        columns=[
            "landcover_status", "classification_reason", "center_match_n",
            "L2_CODE", "L2_NAME", "landcover_tile", "record_id",
            "landcover_x_m", "landcover_y_m",
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "outputs")
    out = parser.parse_args().output_dir
    if out.exists() and any(out.iterdir()):
        parser.error("--output-dir must be empty")

    candidate_path = ROOT / "week06" / "step2" / "recommended_candidate_points.csv"
    candidates = pd.read_csv(candidate_path, dtype={"shot_number": str})
    feb = candidates.loc[
        candidates.region_id.eq("G025_515_140")
        & candidates.observation_date_utc.eq("2025-02-04")
        & candidates.orbit.eq("O34820")
    ].copy()
    if (
        len(feb) != 533
        or not feb.point_id.is_unique
        or not feb.quality_valid.eq(True).all()
        or not np.isfinite(feb[["longitude", "latitude", "agbd_Mg_ha"]]).all().all()
    ):
        raise ValueError("February candidates differ from the expected 533 quality-filtered points")

    spatial_index = load_polygons()
    feb = feb.join(classify_points(feb, *spatial_index))
    feb["high_agbd_ge100"] = feb.agbd_Mg_ha.ge(100)
    april = pd.read_csv(WEEK07 / "outputs" / "landcover_points.csv", dtype={"L2_CODE": str})
    if len(april) != 441 or not april.point_id.is_unique or feb.point_id.isin(april.point_id).any():
        raise ValueError("Unexpected April reference or overlap with February candidate IDs")
    columns = ["landcover_status", "classification_reason", "L2_CODE", "L2_NAME"]
    pd.testing.assert_frame_equal(
        classify_points(april, *spatial_index)[columns].fillna(""),
        april[columns].fillna(""),
    )

    def counts(group):
        return {"n": len(group), "high_agbd_n": int(group.high_agbd_ge100.sum())}

    summary = {
        "candidate_source": candidate_path.relative_to(ROOT).as_posix(),
        "region_id": "G025_515_140",
        "date_utc": "2025-02-04",
        "orbit": "O34820",
        "high_agbd_threshold_Mg_ha": 100,
        "total": counts(feb),
        "status_counts": {
            status: counts(feb.loc[feb.landcover_status.eq(status)])
            for status in ("forest", "nonforest", "unclassified")
        },
        "forest_code_counts": {
            code: counts(feb.loc[feb.L2_CODE.eq(code)]) for code in sorted(FOREST)
        },
        "classification_reasons": feb.classification_reason.value_counts().to_dict(),
        "week07_april_reference_matches": len(april),
        "sentinel2_screening": "pending",
    }
    out.mkdir(parents=True, exist_ok=True)
    feb.to_csv(out / "february_landcover_points.csv", index=False, encoding="utf-8-sig")
    summary_json = json.dumps(summary, ensure_ascii=False, indent=2)
    (out / "february_forest_summary.json").write_text(summary_json + "\n", encoding="utf-8")
    print(summary_json)


if __name__ == "__main__":
    main()
