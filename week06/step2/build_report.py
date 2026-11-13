"""Build offline candidate manifests, plots, and Korean report from the inventory."""
from pathlib import Path
import json
import os
import hashlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["MPLCONFIGDIR"] = str(HERE / ".matplotlib")
import numpy as np
import pandas as pd
import h5py


def md(table):
    cols = list(table.columns)
    return "\n".join(["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
                     + ["| " + " | ".join(str(v) for v in row) + " |" for row in table.itertuples(index=False, name=None)])


def run():
    s = json.loads((HERE / "summary.json").read_text(encoding="utf-8"))
    p = pd.read_csv(HERE / "quality_valid_points.csv", dtype={"shot_number": "string"})
    assignments = pd.read_csv(HERE / "point_region_assignment.csv")
    p = p.merge(assignments, on="point_id", validate="one_to_one")
    grids = pd.read_csv(HERE / "regions_korea_window_025deg.csv").set_index("region_id")
    # Explicit exploratory selection: near-area expansion, a farther NE area,
    # a distant SW area, and another orbit/date NW of the baseline.
    choices = [
        ("A", "G025_515_140", "금정산 서쪽 좌표권", "1595점으로 조사창 내 최다; 3궤도·3날짜. 낮은 AGBD 중앙값(13.34) 때문에 산림 통과율을 먼저 확인"),
        ("B", "G025_516_141", "금정산 북쪽 좌표권", "799점, 3궤도·3날짜; 기존 지역에 인접하여 처리 절차를 확장하기 용이"),
        ("C", "G025_517_142", "금정산 동북쪽 좌표권", "1109점, 2궤도·2날짜; A·B보다 떨어진 공간 검증 후보"),
        ("D", "G025_506_133", "금정산 남서쪽 원거리 좌표권", "841점, 2궤도·2날짜; 지역 간 전이 검증용으로 우선 조사. 산림·계절 차이는 별도 확인"),
        ("E", "G025_513_143", "금정산 북서쪽 좌표권", "297점, 기존 O36071과 다른 O34794; 2월 단일 시기·낮은 AGBD 중앙값(9.58)이므로 예비 후보"),
    ]
    candidates = []
    requests = []
    features = []
    for label, rid, desc, reason in choices:
        r = grids.loc[rid].to_dict()
        r.update(candidate=label, region_id=rid, coordinate_label=desc, selection_reason=reason,
                 administrative_or_forest_boundary=False, final_training_ready_count=None)
        candidates.append(r)
        part = p.loc[p.region_id.eq(rid)]
        # Request box = occupied point bounds + 0.01 degree on each side (~1 km).
        bbox = [round(part.longitude.min() - .01, 6), round(part.latitude.min() - .01, 6),
                round(part.longitude.max() + .01, 6), round(part.latitude.max() + .01, 6)]
        for day, group in part.groupby("observation_date_utc"):
            date = pd.Timestamp(day)
            requests.append({"candidate": label, "region_id": rid, "gedi_date_utc": day,
                "gedi_orbits": ";".join(sorted(group.orbit.unique())), "gedi_quality_count": len(group),
                "bbox_west": bbox[0], "bbox_south": bbox[1], "bbox_east": bbox[2], "bbox_north": bbox[3],
                "preferred_search_start_utc": (date - pd.Timedelta(days=7)).strftime("%Y-%m-%d"),
                "preferred_search_end_utc": (date + pd.Timedelta(days=7)).strftime("%Y-%m-%d"),
                "collection": "Sentinel-2 L2A", "bands": "B02;B03;B04;B08;B11;B12;SCL;dataMask",
                "status": "planning_only_no_search_or_download", "product_id": "not_yet_identified"})
        w, so, e, n = (r[k] for k in ["grid_west", "grid_south", "grid_east", "grid_north"])
        features.append({"type": "Feature", "properties": {"candidate": label, "region_id": rid,
            "quality_valid": int(r["quality_valid"]), "boundary_type": "arbitrary_0.25_degree_grid"},
            "geometry": {"type": "Polygon", "coordinates": [[[w, so], [e, so], [e, n], [w, n], [w, so]]]}})
    c = pd.DataFrame(candidates)
    c.to_csv(HERE / "recommended_candidates.csv", index=False)
    req = pd.DataFrame(requests)
    req.to_csv(HERE / "sentinel2_acquisition_plan.csv", index=False)
    (HERE / "recommended_candidate_grids.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": features}, indent=2), encoding="utf-8")
    selected = p.loc[p.region_id.isin(c.region_id)]
    selected.to_csv(HERE / "recommended_candidate_points.csv", index=False)
    regional = p.loc[p.in_korea_screening_window]
    regional.to_csv(HERE / "korea_window_quality_valid_points.csv", index=False)
    region_dates = regional.groupby(["observation_date_utc", "orbit"]).agg(
        quality_valid=("point_id", "size"), existing_90=("existing_week04_05_90", "sum"),
        within_original_rectangle=("within_original_sentinel_rectangle", "sum")).reset_index()
    region_dates["additional_vs_90"] = region_dates.quality_valid - region_dates.existing_90
    region_dates["outside_original_rectangle"] = region_dates.quality_valid - region_dates.within_original_rectangle
    region_dates.to_csv(HERE / "counts_korea_window_date_orbit.csv", index=False)

    # Independent second scan: a different chunk size and directly expressed
    # equivalent predicates, plus per-granule cumulative stages. No NxN matrix.
    checks = []
    for path in sorted((ROOT / "week03/GEDI04_A").glob("*.h5")):
        counts = np.zeros(8, dtype=np.int64)
        with h5py.File(path, "r") as f:
            for beam in sorted(k for k in f if k.startswith("BEAM")):
                g = f[beam]
                for start in range(0, len(g["shot_number"]), 131_072):
                    sl = slice(start, start + 131_072)
                    a, se = g["agbd"][sl], g["agbd_se"][sl]
                    masks = [np.isfinite(a) & (a >= 0), np.isfinite(se) & (se >= 0),
                             g["l2_algrunflag"][sl] == 1, g["l2a_quality_flag_rel3"][sl] == 1,
                             g["l4a_quality_flag_rel3"][sl] == 1, g["degrade_include_flag"][sl] == 1,
                             g["elev_highestreturn_outlier_flag"][sl] == 0]
                    ok = np.ones(len(a), dtype=bool)
                    counts[0] += len(a)
                    for j, mask in enumerate(masks, 1):
                        ok &= mask
                        counts[j] += int(ok.sum())
        checks.append(dict(source_file=path.name, **dict(zip(
            ["all", "agbd_nonnegative_finite", "se_nonnegative_finite", "algorithm_run", "l2a_quality", "l4a_quality", "degrade_include", "no_elevation_outlier"], map(int, counts)))))
    independent = pd.DataFrame(checks)
    independent.to_csv(HERE / "quality_filter_counts_by_granule.csv", index=False)
    expected = pd.read_csv(HERE / "counts_by_granule.csv").set_index("source_file")
    assert independent.set_index("source_file")["all"].sort_index().equals(expected["all"].sort_index())
    assert independent.set_index("source_file").no_elevation_outlier.sort_index().equals(expected.quality_valid.sort_index())
    assert int(c.quality_valid.sum()) == len(selected) == 4641
    assert selected.point_id.nunique() == 4641 and not selected.existing_week04_05_90.any()
    assert selected.within_original_sentinel_rectangle.sum() == 0
    assert req.gedi_quality_count.sum() == 4641
    assert region_dates.quality_valid.sum() == 12444
    assert s["quality_valid"] == s["existing_week04_05_90"] + s["additional_quality_valid_vs_90"]
    assert s["additional_quality_valid_vs_90"] == s["outside_original_rectangle"] + 7
    for filename, key in [(ROOT / "week03/analyze_gedi.py", "source_quality_script_sha256"),
                          (ROOT / "week03/outputs/sentinel_2025/scl_validation/gedi_scl_valid_points.csv", "baseline_90_sha256"),
                          (ROOT / "week03/outputs/gedi/gedi_valid_points.csv", "old_97_sha256")]:
        assert hashlib.sha256(filename.read_bytes()).hexdigest() == s[key]
    validation = {"status": "passed", "independent_hdf5_chunk_size": 131072,
        "per_granule_all_and_quality_counts_match": True, "existing_90_and_97_recovered_exactly": True,
        "candidate_points_unique": 4641, "candidate_points_outside_original_rectangle": 4641,
        "candidate_date_counts_reconcile": True, "source_and_existing_csv_hashes_unchanged_since_scan": True,
        "grid_count_sum": int(pd.read_csv(HERE / "regions_global_025deg.csv").quality_valid.sum()),
        "component_count_sum": int(pd.read_csv(HERE / "occupied_grid_components.csv").quality_valid.sum())}
    assert validation["grid_count_sum"] == validation["component_count_sum"] == 503053
    (HERE / "validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig, ax = plt.subplots(figsize=(11, 9), constrained_layout=True)
    for orbit, group in regional.groupby("orbit"):
        ax.scatter(group.longitude, group.latitude, s=2, alpha=.6, label=orbit)
    base = p.loc[p.existing_week04_05_90]
    ax.scatter(base.longitude, base.latitude, s=15, marker="x", c="black", label="Existing 90")
    for row in c.itertuples():
        ax.add_patch(Rectangle((row.grid_west, row.grid_south), .25, .25, fill=False, color="black", lw=1.2))
        ax.annotate(f"{row.candidate}: {row.quality_valid:,}", (row.grid_east, row.grid_north), xytext=(5, 5), textcoords="offset points", fontsize=10,
                    bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none"})
    ax.set(xlim=(126.2, 130), ylim=(33, 36.5), xlabel="Longitude (WGS84 degrees)", ylabel="Latitude (WGS84 degrees)",
           title="Expansion candidates: arbitrary 0.25-degree cells\nCounts = GEDI quality only; no administrative/forest boundaries or basemap")
    ax.grid(alpha=.3)
    ax.legend()
    fig.savefig(HERE / "recommended_candidates.png", dpi=160)
    plt.close(fig)

    dates = pd.read_csv(HERE / "counts_by_date.csv")
    dates["궤도"] = pd.read_csv(HERE / "counts_by_orbit.csv").orbit
    dates = dates.rename(columns={"observation_date_utc": "관측일 UTC", "all": "전체", "quality_valid": "품질 통과"})
    date_table = md(dates[["관측일 UTC", "궤도", "전체", "품질 통과"]])
    beams = pd.read_csv(HERE / "counts_by_beam.csv")
    beam_table = md(beams[["beam", "all", "quality_valid"]])
    ct = c[["candidate", "region_id", "grid_west", "grid_east", "grid_south", "grid_north", "quality_valid", "date_count", "orbit_count"]]
    plan_table = md(req[["candidate", "gedi_date_utc", "gedi_orbits", "gedi_quality_count", "preferred_search_start_utc", "preferred_search_end_utc"]])
    rt = md(region_dates[["observation_date_utc", "orbit", "quality_valid", "additional_vs_90", "outside_original_rectangle"]])
    report = f"""# Week06 Step 2 — 보유 GEDI L4A V003의 연구 범위 확대 조사

작성 기준: 2026-11-13. 로컬 HDF5 5개 전수 조사 결과이며 새 다운로드·외부 API 호출·모델 학습을 하지 않았다. Week03–Week05 코드·산출물을 수정하지 않았다.

## 1. 결론과 표본의 의미

- 전체 **6,356,011점** 중 기존 품질 조건 통과는 **503,053점(7.9146%)**이다.
- 기존 Week04·Week05의 90점과 정확히 일치하는 식별자를 제외하면 **502,963점**이다. 이 중 **502,956점은 기존 Sentinel-2 사각형 밖**, 나머지 **7점은 기존 사각형 안에서 이전 SCL 검증에 탈락한 점**이다. 이 7점은 새 학습점으로 재수용하지 않는다.
- 국내 주변 조사창 [124,132)°E × [33,39)°N 안에는 **12,444점**, 기존 90점을 제외하면 **12,354점**, 기존 사각형 밖에는 **12,347점**이다. 이 창은 국가·육지·행정 경계가 아니며 일본·해상·북쪽 지역도 포함하므로 대한민국 표본 수로 해석할 수 없다.
- 추천한 서로 겹치지 않는 5개 좌표 격자에는 **4,641점**이 있다. 모두 기존 90점 및 기존 사각형 밖의 추가 GEDI 후보다. 12,347점 및 502,956점의 부분집합이며 서로 더하지 않는다.
- **새 최종 학습 가능점 수는 미확정(null)**이다. 추가 후보에 Sentinel-2 유효 반사도·구름/그림자·산림 경계·시차 조건을 적용하지 않았다. GEDI 품질 통과 자체로 산림이나 독립 현장 AGBD 정답을 보증하지 않는다.

## 2. 입력·재현 방법과 품질 조건

입력은 `../../week03/GEDI04_A/`의 HDF5 5개다. 전체 파일명과 크기는 `summary.json`의 inputs, 파일별 집계는 `counts_by_granule.csv`에 기록했다. HDF5 `short_name=GEDI_L4A`, 파일명 V003, delta_time 단위를 검사했다.

`../../week03/analyze_gedi.py`의 **physical_values()와 quality_mask() 원본 함수 AST를 그대로 실행**했다. 옛 모듈을 import하여 Week03 캐시에 쓰는 부작용을 피했다. source SHA-256: `{s['source_quality_script_sha256']}`.

1. agbd 및 agbd_se의 비유한값·-9999·음수를 결측 처리하고 유한값만 통과
2. l2_algrunflag=1
3. l2a_quality_flag_rel3=1
4. l4a_quality_flag_rel3=1
5. degrade_include_flag=1
6. elev_highestreturn_outlier_flag=0

degrade_flag=0 또는 sensitivity≥0.95는 추가하지 않았다. sensitivity와 원시 degrade_flag는 보존했다. 좌표·시간 유효성은 품질 조건과 별도로 검사했으며 품질 통과점 중 무효 좌표·시간은 각각 **0점**이다.

BEAM 내 최대 **100,000점씩** 읽고 통과점만 CSV에 순차 저장했다. 공간 요약 단계에는 503,053개 통과점만 메모리에 적재한다. 전체 거리행렬·DBSCAN 거리행렬은 만들지 않았다. 입력량에 선형인 스캔·해시 그룹 집계를 사용했다. 글로벌 지도는 고정 시드 표본 160,000점을 표시하며 조사창 지도는 12,444점 전부를 표시한다. 모든 통계는 지도 표본이 아닌 전체 통과점 기준이다.

## 3. 전체 품질 필터 집계

{md(pd.read_csv(HERE / 'quality_filter_counts.csv'))}

cumulative_count는 위에서부터 순서대로 모두 만족한 점 수다. independent_pass_count는 다른 조건과 무관하게 해당 조건 하나를 만족한 수로, 최종 표본 수와 혼동하지 않는다. 6월 9일 파일은 L4A 품질 단계까지 103,371점이 남지만 **degrade_include_flag=1 조건에서 0점**이 된다(`quality_filter_counts_by_granule.csv`). 기존 기준을 완화하지 않았다.

## 4. 날짜·궤도·granule·BEAM 분포

{date_table}

이번 5개 파일에서는 날짜·궤도·granule이 각각 1:1로 대응한다. 파일명에서 Oxxxxx를 궤도 식별자로 읽고 날짜는 파일명 추정이 아니라 delta_time을 2018-01-01 기준 UTC로 변환했다. 품질 통과점의 약 {358068/503053*100:.2f}%가 2월 2일 O34794에 집중되어 있다.

{beam_table}

40개 granule×BEAM 조합의 전체·통과 수(통과 0인 조합 포함)는 `counts_granule_date_orbit_beam.csv`에 있다. `counts_by_granule.csv`, `counts_by_date.csv`, `counts_by_orbit.csv`, `counts_by_beam.csv`는 각각 주변 합계를 제공한다.

국내 주변 조사창 내 날짜·궤도 분포는 다음과 같다.

{rt}

O35415의 조사창 내 통과점은 7점뿐이다. 전 세계 통과점이 38,347점이라는 이유로 국내 봄철 표본이 풍부하다고 해석할 수 없다. 4월 25일 조사창 통과 4,852점 중 기존 90점을 빼면 4,762점, 기존 사각형 밖은 4,755점이다. 기존 연구와 계절을 가깝게 유지하려면 이 4월 집합부터 확인하는 것이 합리적이다.

## 5. 공간 범위와 분포

품질 통과 전체 범위는 경도 **-178.859543–152.843742°**, 위도 **-0.107886–51.833091°**다. 전 세계 궤도 구간을 포함하며 단일 연결 연구 지역이 아니다. 이 최소·최대 범위를 통째로 Sentinel-2 다운로드 영역으로 사용하면 안 된다.

조사창 통과점의 실제 범위는 경도 **124.888030–131.996823°**, 위도 **33.046986–38.323644°**다. 궤도별 가늘고 긴 선형 분포이며, 점 수가 많아도 독립 지역·관측일 수가 같은 비율로 증가하지 않는다.

경도·위도 원점 0° 기준 **0.25°×0.25° 격자**로 전 세계 통과점을 나눴다. ID `G025_x_y`의 서쪽/남쪽 경계는 x×0.25, y×0.25이며 반열린 구간을 사용한다. 국내 주변 위도에서 대략 동서 22–23km, 남북 28km지만 전 세계에서 같은 면적의 격자는 아니다.

- 통과점이 있는 전 세계 격자: **453개** (`regions_global_025deg.csv`)
- 조사창 격자: **{len(grids)}개** (`regions_korea_window_025deg.csv`)
- 점이 있는 격자의 8방향 인접 연결 성분: **26개** (`occupied_grid_components.csv`)

연결 성분은 **인접 격자의 위상적 연결**일 뿐, 실제 점 간 거리 기반 군집·독립 검증 그룹·행정구역·산림 경계가 아니다. 빈 공간과 긴 궤도가 같은 성분으로 이어질 수 있으며 격자 크기·원점을 바꾸면 결과도 달라진다. 먼 소수 점도 포함하므로 전체 범위는 대표적인 밀집 범위와 구분한다.

최대 성분은 312,478점, 93.147–121.028°E / 40.700–50.414°N으로 기존 지역과 크게 떨어져 있다. 국내 연구의 단순 연장으로 바로 합치기보다 기후·식생·품질 검토가 필요하다. 국가·행정·산림 지오메트리가 없으므로 성분에 국가명이나 행정구역명을 확정 부여하지 않았다.

![전체 및 조사창 분포](spatial_distribution.png)

## 6. 추천 조사 지역

추천은 표본 수, 기존 지역과의 상대 위치, 날짜·궤도 다양성, 후속 영상 확보 범위를 고려한 **탐색적 선정**이다. 산림 적합도 순위가 아니다. 표의 서·동·남·북은 WGS84 도 단위 격자 경계다. 라벨은 방향 설명이며 실제 행정구역명이 아니다.

{md(ct)}

"""
    for label, rid, desc, reason in choices:
        report += f"- **{label} — {desc}**: {reason}.\n"
    report += f"""

A·B로 기존 처리 절차를 넓힌 뒤 C의 공간 분리 평가, D의 원거리 전이 평가로 확장할 수 있다. E는 다른 궤도를 보강하는 예비 후보다. A·B는 인접하므로 별도 지역 라벨을 붙였다는 이유만으로 독립 train/test로 나누지 않는다. C·D도 실제 학습·검증 최소 거리와 버퍼를 확인해야 한다.

A–E 합계 **4,641점**은 모두 기존 사각형 밖이다. 후보별 날짜·궤도·BEAM 집계는 `counts_region_date_orbit_beam.csv`, 실제 좌표와 AGBD·품질 정보는 `recommended_candidate_points.csv`에서 확인한다. 숫자는 GEDI 관측점 수이며 동일 Sentinel-2 픽셀 수나 독립 표본 수가 아니다.

![추천 격자](recommended_candidates.png)

## 7. 다음 단계 Sentinel-2 범위

`sentinel2_acquisition_plan.csv`에는 **{len(req)}개 후보×관측일 조합**의 계획을 기록했다. 외부 카탈로그를 조회하지 않았으므로 실제 제품 ID·구름량·장면 확보 여부·타일 수는 미확정이다.

각 후보의 **통과점 실제 최소·최대 좌표에 사방 0.01° 여유**를 둔 bbox를 기록했다. 이는 약 1km 수준의 처리 여유이며 GEDI footprint 경계가 아니다. 여러 타일에 걸리면 같은 기준으로 정렬할 필요가 있다.

{md(req[['candidate','bbox_west','bbox_south','bbox_east','bbox_north']].drop_duplicates())}

우선 GEDI 관측일과 가장 가까운 장면을 찾도록 **±7일을 초기 검색 구간**으로 제안한다. 이는 확정 품질 기준이 아니며 구름 때문에 확보하지 못하면 허용 시차와 계절·생육 변화의 영향을 검토한 뒤 확장해야 한다. 특히 2월·3월·4월을 섞을 때 낙엽·생육 차이를 설명변수 또는 분석 설계에 반영해야 한다.

{plan_table}

필요 자료: **Sentinel-2 L2A B02·B03·B04·B08·B11·B12, 동일 제품 SCL, 유효 마스크와 제품 메타데이터**. NDVI·EVI·NDMI에는 B02/B04/B08/B11이 필요하고 B03/B12는 RGB·추가 변수·기존 실험 호환용이다. 반사도 밴드는 B04 기준 목표 격자에 정렬하고 SCL은 최근린을 사용한다. 제품 ID·UTC 촬영시각·단위·scale/offset 또는 evalscript·출력 형식·처리 기준·해시를 보존한다. 기존 `/65535`를 새 자료에 자동 적용하지 않는다.

국내 주변 조사창 전체나 전 세계 min/max를 일괄 수집하지 말고, **A–D의 4월 25일 후보 2,526점**에 대응하는 장면부터 조사하여 기존 봄철 자료와 비교한 뒤 다른 날짜로 확대하는 순서를 권한다. 이 2,526점은 4,641점의 부분집합이며 E에는 4월 후보가 없다. 산림 경계·토지피복 자료는 별도로 필요하며 SCL=4만으로 산림을 확정하지 않는다.

## 8. 검증과 산출물 사용법

`validation.json`의 검증은 모두 통과했다. 다른 크기(131,072점) 청크로 HDF5를 다시 읽고 품질 조건을 독립 표현하여 granule별 전체·통과 수를 대조했다. 기존 90점·97점은 파일명+BEAM+정수 shot_number 문자열로 정확히 재식별했다. 기존 사각형 내부 통과 97점을 재현했고 새 후보 4,641점의 중복·범위·날짜별 합계를 확인했다. 전 세계 격자 합계와 연결 성분 합계는 각각 503,053점이다. 원본 품질 코드와 기존 두 CSV의 SHA-256은 조사 전후 동일하다.

주요 파일:

| 파일 | 용도 |
| --- | --- |
| summary.json | 전체·추가·조사창 수, 좌표 범위, 입력 파일 목록 |
| quality_valid_points.csv | 전 세계 품질 통과 503,053점, 약 201MB; 기존 90/97 중복·기존 사각형 포함 플래그 |
| korea_window_quality_valid_points.csv | 국내 주변 **사각형** 통과 12,444점; 국내 확정 표본 아님 |
| counts_granule_date_orbit_beam.csv | 0점 조합을 포함한 세부 집계 |
| quality_filter_counts.csv / quality_filter_counts_by_granule.csv | 전체 및 파일별 단계별 품질 통과 수 |
| regions_global_025deg.csv / regions_korea_window_025deg.csv | 임의 좌표 격자별 수·범위·날짜·궤도 |
| point_region_assignment.csv / point_component_assignment.csv | point_id와 격자·연결 성분 연결 |
| occupied_grid_components.csv | 임의 인접 격자 연결 성분 요약 |
| recommended_candidates.csv / recommended_candidate_points.csv | A–E 후보 정의·선정 이유 및 4,641개 실제 점 |
| recommended_candidate_grids.geojson | GIS에서 열 수 있는 임의 후보 격자; 행정·산림 경계 아님 |
| sentinel2_acquisition_plan.csv | 후보×날짜별 bbox·초기 검색 기간; 조회·다운로드 미실행 |
| validation.json | 집계 대조·기존 표본 재현·변경 여부 검증 |

재실행(기존 산출물은 읽기만 하고 Step2 결과를 재작성):

```powershell
& week03/.venv/Scripts/python.exe -B week06/step2/analyze_expansion.py
& week03/.venv/Scripts/python.exe -B week06/step2/build_report.py
```

다음 실행 순서: 후보·계절·허용 시차 확정 → 계획 bbox와 날짜에 맞는 Sentinel-2 및 산림 자료 확보 → 장면별 반사도 검증·마스크·지수 계산 → GEDI 결합 → 산림·시간·중복 픽셀 필터 → 지역·궤도·날짜별 **최종 학습 가능점** 재집계 → 공간 버퍼·미사용 지역 홀드아웃 확정 → 이후 모델 평가.
"""
    (HERE / "report.md").write_text(report, encoding="utf-8")
    print(json.dumps(validation, indent=2), flush=True)
    print(f"Report built; candidates={len(selected):,}; plan rows={len(req)}", flush=True)


if __name__ == "__main__":
    run()
