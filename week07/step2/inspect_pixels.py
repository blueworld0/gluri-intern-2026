"""Offline radiometry and local-pixel audit. No fitting, filtering, or downloads."""
from pathlib import Path
from contextlib import ExitStack
import argparse
import hashlib
import json
import os
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window, from_bounds

BANDS = ['B02', 'B03', 'B04', 'B08', 'B11', 'B12']


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def md(df):
    def fmt(v):
        return f'{v:.3f}' if isinstance(v, (float, np.floating)) else str(v)
    return '| ' + ' | '.join(df.columns) + ' |\n| ' + ' | '.join(['---'] * len(df.columns)) + ' |\n' + '\n'.join(
        '| ' + ' | '.join(fmt(v) for v in row) + ' |' for row in df.itertuples(index=False, name=None))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).parent / 'outputs')
    out = parser.parse_args().output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Output directory must be empty; specify a new --output-dir.')
    base = ROOT / 'week06/step5'
    safe = next((base / 'source').glob('*.SAFE'))
    granule = next((safe / 'GRANULE').iterdir())
    raw_csv = base / 'outputs/gedi_sentinel_all_points.csv'
    diagnostic_csv = ROOT / 'week07/step1/outputs/all441_error_diagnostics.csv'
    paths = [raw_csv, diagnostic_csv, safe / 'MTD_MSIL2A.xml'] + sorted((base / 'outputs/rasters').glob('*.tif'))
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    xml = ET.parse(safe / 'MTD_MSIL2A.xml')
    local = lambda t: t.rsplit('}', 1)[-1]
    q = float(next(e.text for e in xml.iter() if local(e.tag) == 'BOA_QUANTIFICATION_VALUE'))
    ids = {e.attrib['bandId']: e.attrib['physicalBand'] for e in xml.iter() if local(e.tag) == 'Spectral_Information'}
    offsets = {ids[e.attrib['band_id']]: float(e.text) for e in xml.iter() if local(e.tag) == 'BOA_ADD_OFFSET'}
    offsets = {b: offsets['B' + str(int(b[1:]))] for b in BANDS}
    diag = pd.read_csv(diagnostic_csv)
    source = pd.read_csv(raw_csv, dtype={'shot_number': str})
    assert len(diag) == 441 and diag.point_id.is_unique and source.point_id.is_unique
    df = diag.merge(source, on='point_id', validate='one_to_one', suffixes=('', '_source'), indicator=True)
    assert df['_merge'].eq('both').all()
    df = df.drop(columns='_merge').copy()
    for name in ['ndvi', 'evi', 'ndmi']:
        assert np.allclose(df[name], df[name + '_source'], atol=1e-7, rtol=0)
    coords = list(zip(df.easting_m, df.northing_m))
    bounds = (df.easting_m.min()-150, df.northing_m.min()-150, df.easting_m.max()+150, df.northing_m.max()+150)
    raw_stats = {}

    def native_read(path, indexes=1):
        raw_stats[str(path.relative_to(ROOT))] = [path.stat().st_size, path.stat().st_mtime_ns]
        with rasterio.open(path) as ds:
            assert ds.crs.to_epsg() == 32652
            win = from_bounds(*bounds, ds.transform).round_offsets().round_lengths()
            data = ds.read(indexes, window=win)
            t = ds.window_transform(win)
            rr, cc = rasterio.transform.rowcol(t, df.easting_m.to_numpy(), df.northing_m.to_numpy())
            assert min(rr) >= 0 and min(cc) >= 0 and max(rr) < data.shape[-2] and max(cc) < data.shape[-1]
            return data, t, np.asarray(rr), np.asarray(cc)

    checks = []
    for b in BANDS:
        res = '20' if b in ['B11', 'B12'] else '10'
        raw, _, rr, cc = native_read(next((granule / f'IMG_DATA/R{res}m').glob(f'*_{b}_{res}m.jp2')))
        sampled = raw[rr, cc]
        assert np.array_equal(sampled, df[b + '_native_dn'])
        quality, _, qr, qc = native_read(granule / f'QI_DATA/MSK_QUALIT_{b}.jp2', list(range(1, 9)))
        detector, _, dr, dc = native_read(granule / f'QI_DATA/MSK_DETFOO_{b}.jp2')
        df[b + '_native_quality_bad_rechecked'] = (quality[:, qr, qc] != 0).any(axis=0)
        df[b + '_native_detector_rechecked'] = detector[dr, dc]
        assert np.array_equal(detector[dr, dc], df[b + '_native_detector_id'])
        assert not df[b + '_native_quality_bad_rechecked'].any()
        assert np.isin(sampled, [0, 65535]).sum() == 0
        assert np.isin(detector[dr, dc], np.arange(1, 13)).all()
        rho = (sampled.astype(float) + offsets[b]) / q
        df[b + '_native_reflectance_rechecked'] = rho
        if res == '10':
            assert np.allclose(rho, df['reflectance_' + b], atol=1e-7, rtol=0)
        checks.append(dict(band=b, native_resolution_m=int(res), dn_csv_match=True,
                           center_quality_bad_n=0, offset=offsets[b], quantification=q))
        print(f'Native DN and quality verified: {b}', flush=True)
    scl, st, sr, sc = native_read(next((granule / 'IMG_DATA/R20m').glob('*_SCL_20m.jp2')))
    assert np.array_equal(scl[sr, sc], df.scl_native_20m)
    classi, _, cr, cc = native_read(granule / 'QI_DATA/MSK_CLASSI_B00.jp2', [1, 2, 3])
    assert not (classi[:, cr, cc] != 0).any()
    df['evi_denominator_rechecked'] = df.reflectance_B08 + 6 * df.reflectance_B04 - 7.5 * df.reflectance_B02 + 1
    df['evi_formula_rechecked'] = 2.5 * (df.reflectance_B08 - df.reflectance_B04) / df.evi_denominator_rechecked
    assert np.allclose(df.evi, df.evi_formula_rechecked, atol=.001, rtol=1e-5)
    df['evi_review_rechecked'] = (df.evi.abs() > 2) | (df.evi_denominator_rechecked.abs() < .05)
    df['relative_agbd_se'] = df.agbd_se / df.agbd
    out.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        rasters = {b: stack.enter_context(rasterio.open(base / f'outputs/rasters/reflectance_{b}.tif')) for b in BANDS}
        ndvi_ds = stack.enter_context(rasterio.open(base / 'outputs/rasters/ndvi.tif'))
        for b, ds in rasters.items():
            vals = np.array([v[0] for v in ds.sample(coords)])
            assert np.allclose(vals, df['reflectance_' + b], atol=1e-7, rtol=0)
        rows = []
        selected = df[(df.error_rank <= 20) | df.evi_review_rechecked]
        review_panels = []
        for row in df.itertuples():
            result = {'point_id': row.point_id}
            # Native SCL fractions use pixel centers, not area-weighted footprint overlap.
            for radius in [12.5, 50.0]:
                rr, cc = np.meshgrid(np.arange(max(0, sr[row.Index]-4), min(scl.shape[0], sr[row.Index]+5)),
                                     np.arange(max(0, sc[row.Index]-4), min(scl.shape[1], sc[row.Index]+5)), indexing='ij')
                xx = st.c + (cc + .5) * st.a
                yy = st.f + (rr + .5) * st.e
                mask = (xx-row.easting_m)**2 + (yy-row.northing_m)**2 <= radius**2
                values = scl[rr[mask], cc[mask]]
                prefix = 'scl_r12p5' if radius == 12.5 else 'scl_r50'
                result[prefix + '_pixel_n'] = len(values)
                for label, classes in [('vegetation', [4]), ('nonvegetated', [5]), ('water', [6]), ('cloud_shadow', [2,3,8,9,10]), ('other', [0,1,7,11])]:
                    result[prefix + '_' + label + '_fraction'] = float(np.isin(values, classes).mean()) if len(values) else np.nan
            r, c = ndvi_ds.index(row.easting_m, row.northing_m)
            win = Window(c-5, r-5, 11, 11)
            arr = ndvi_ds.read(1, window=win, masked=True)
            yy, xx = np.meshgrid(np.arange(r-5, r+6), np.arange(c-5, c+6), indexing='ij')
            t = ndvi_ds.transform
            disk = ((t.c+(xx+.5)*t.a-row.easting_m)**2 + (t.f+(yy+.5)*t.e-row.northing_m)**2) <= 50**2
            vals = arr[disk].compressed()
            result.update(ndvi_r50_valid_n=len(vals), ndvi_r50_total_n=int(disk.sum()),
                          ndvi_r50_mean=float(vals.mean()), ndvi_r50_sd=float(vals.std()),
                          ndvi_r50_min=float(vals.min()), ndvi_r50_max=float(vals.max()))
            rows.append(result)
            if row.point_id in set(selected.point_id):
                win = Window(c-10, r-10, 21, 21)
                rgb = np.stack([rasters[b].read(1, window=win, masked=True).filled(np.nan) for b in ['B04','B03','B02']], axis=-1)
                rgb = np.nan_to_num(np.clip(rgb / .3, 0, 1)) ** (1/2.2)
                nd = ndvi_ds.read(1, window=win, masked=True)
                # SCL overview stays on its native 20m grid.
                ri, ci = rasterio.transform.rowcol(st, row.easting_m, row.northing_m)
                sw = scl[ri-5:ri+6, ci-5:ci+6]
                rgb_t = rasters['B04'].window_transform(win)
                extent = [rgb_t.c-row.easting_m, rgb_t.c+210-row.easting_m,
                          rgb_t.f-210-row.northing_m, rgb_t.f-row.northing_m]
                scl_extent = [st.c+(ci-5)*20-row.easting_m, st.c+(ci+6)*20-row.easting_m,
                              st.f-(ri+6)*20-row.northing_m, st.f-(ri-5)*20-row.northing_m]
                review_panels.append((row, rgb, nd, sw, extent, scl_extent))
    df = df.merge(pd.DataFrame(rows), on='point_id', validate='one_to_one')
    df.to_csv(out / 'all441_pixel_audit.csv', index=False)
    selected = df[(df.error_rank <= 20) | df.evi_review_rechecked]
    selected.to_csv(out / 'review21_pixel_audit.csv', index=False)
    cmap = ListedColormap(['black','#ff0000','#555555','#8b4513','#228b22','#e0c078','#1488cc','#888888','#bbbbbb','#eeeeee','#aaccff','#ffccff'])
    for start in range(0, len(review_panels), 7):
        panels = review_panels[start:start+7]
        fig, axes = plt.subplots(len(panels), 3, figsize=(12, 3.4*len(panels)), squeeze=False)
        fig.subplots_adjust(left=.07, right=.98, top=.965, bottom=.03, wspace=.35, hspace=.7)
        for axrow, (row, rgb, nd, sw, extent, se) in zip(axes, panels):
            axrow[0].imshow(rgb, extent=extent, interpolation='nearest')
            axrow[1].imshow(nd, extent=extent, vmin=-1, vmax=1, cmap='RdYlGn', interpolation='nearest')
            axrow[2].imshow(sw, extent=se, cmap=cmap, norm=BoundaryNorm(np.arange(-.5,12.5),12), interpolation='nearest')
            for ax in axrow:
                ax.plot(0, 0, '+', color='magenta', ms=9)
                ax.add_patch(Circle((0,0),12.5, fill=False, color='white', linewidth=.9))
                ax.add_patch(Circle((0,0),50, fill=False, color='cyan', linewidth=.6))
                ax.set(xlim=(-100,100), ylim=(-100,100), xlabel='East offset (m)', ylabel='North offset (m)')
            axrow[0].set_title(f'Rank {row.error_rank}: RGB\nAGBD {row.agbd:.1f}, predicted {row.prediction:.1f}', fontsize=10)
            axrow[1].set_title(f'NDVI {row.ndvi:.3f}; EVI {row.evi:.3f}', fontsize=10)
            axrow[2].set_title('Native SCL: green=4, tan=5, blue=6', fontsize=9)
        fig.suptitle('Local review: fixed RGB stretch 0..0.3, gamma 2.2; white r=12.5m, cyan r=50m', fontsize=10)
        fig.savefig(out / f'pixel_review_{start//7+1:02d}.png', dpi=140, bbox_inches='tight')
        plt.close(fig)
    groups = []
    for label, subset in [('top20', df[df.error_rank <= 20]), ('other421', df[df.error_rank > 20]),
                          ('high_agbd_all56', df[df.agbd >= 100]), ('evi_review', df[df.evi_review_rechecked])]:
        groups.append(dict(group=label, n=len(subset), mean_ndvi=float(subset.ndvi.mean()),
                           mean_evi=float(subset.evi.mean()), median_relative_agbd_se=float(subset.relative_agbd_se.median()),
                           mean_scl_vegetation_r50=float(subset.scl_r50_vegetation_fraction.mean()),
                           any_nonvegetation_r12p5_n=int((subset.scl_r12p5_vegetation_fraction < 1).sum()),
                           any_cloud_shadow_r50_n=int((subset.scl_r50_cloud_shadow_fraction > 0).sum()),
                           evi_review_n=int(subset.evi_review_rechecked.sum())))
    comparison = pd.DataFrame(groups)
    comparison.to_csv(out / 'group_comparison.csv', index=False)
    top = df[df.error_rank <= 20]
    review = df[df.evi_review_rechecked]
    summary = dict(n=441, reviewed_n=len(selected), top20_center_native_quality_bad_n=0,
                   top20_evi_review_n=int(top.evi_review_rechecked.sum()),
                   top20_any_nonvegetation_r12p5_n=int((top.scl_r12p5_vegetation_fraction < 1).sum()),
                   top20_any_cloud_shadow_r50_n=int((top.scl_r50_cloud_shadow_fraction > 0).sum()),
                   evi_review_error_ranks=review.error_rank.astype(int).tolist(), forest_status='not_assessed',
                   fitting_performed=False, samples_removed=0)
    report = f'''# Week07 Step 2 — 원본 반사율·주변 픽셀 점검

## 범위와 검증

전체 441점의 원본 SAFE DN·밴드별 8층 품질 플래그·검출기·SCL·CLASSI를 다시 읽었다. 기존 CSV와 DN이 일치하고 중심점 품질 플래그는 모두 정상이다. XML에서 quantification={q:g}, 밴드별 offset={offsets}을 확인했다. 10m 밴드의 (DN+offset)/quantification과 기존 반사율이 일치한다. B11/B12는 원본 20m 중심 DN과 CSV DN을 대조했으며, 10m 정렬 반사율은 저장된 처리 GeoTIFF와 대조했다. 이 단계에서 20→10m bilinear 재투영 전체를 독립 재구현한 것은 아니다. 지수 입력 반사율과 EVI 공식을 재계산했다.

오차 상위 20점과 EVI 점검 대상의 합집합 {len(selected)}점에 대해 RGB·NDVI·원본 SCL 상세도를 만들었다. 전체 441점에도 동일한 주변 통계를 계산해 비교했다. 새 데이터 다운로드·학습·표본 삭제는 하지 않았다.

## 주요 결과

- 상위 20점의 중심 원본 DN·품질 마스크에서 기록 불일치나 불량 플래그는 발견되지 않았다.
- 상위 20점의 EVI 점검 규칙(abs(EVI)>2 또는 abs(분모)<0.05) 해당은 {summary['top20_evi_review_n']}점이다.
- EVI 점검 대상의 절대오차 순위는 {summary['evi_review_error_ranks']}이다. 이 이상점을 고 AGBD 과소 예측의 직접 원인으로 볼 근거는 없다.
- 반경 12.5m 안의 SCL 픽셀 중심 중 식생(4) 이외가 있는 상위 오차점은 {summary['top20_any_nonvegetation_r12p5_n']}/20점이다.
- 반경 50m 안 SCL 구름·그늘(2/3/8/9/10)이 하나라도 있는 상위 오차점은 {summary['top20_any_cloud_shadow_r50_n']}/20점이다.
- 상위 20점 중 반경 50m SCL 픽셀 중심이 모두 식생인 점은 {int((top.scl_r50_vegetation_fraction == 1).sum())}점이다. 주변 혼합 피복만으로 모든 큰 오차를 설명할 수 없다.
- 오차 순위 1·7·9·17 주변 50m에는 SCL 수역 픽셀이 포함된다. 특히 1위는 식생 50%, 수역 35%, 비식생 15%이며, 경계 인접 표본으로 후속 점검할 가치가 있다. 이 비율은 면적 비율이 아니다.

[상세도 1: 오차 순위 1~7](pixel_review_01.png) · [상세도 2: 8~14](pixel_review_02.png) · [상세도 3: 15~20 및 EVI 점검점](pixel_review_03.png)

{md(comparison)}

각 집단은 AGBD 구성과 위치가 다르며, high_agbd_all56에는 top20이 포함된다. 평균 차이를 원인 효과로 해석하지 않는다. 반경 12.5m에 원본 20m 픽셀 중심이 없는 점은 전체 {int(df.scl_r12p5_pixel_n.eq(0).sum())}점이며 해당 비율을 결측으로 보존했다(상위 20점에는 없음). 결측은 비식생 발견 개수에 포함하지 않으며 정상 피복 판정으로 해석하지 않는다. `relative_agbd_se`는 GEDI 표준오차/AGBD 비율이며 현장 검증 오차나 품질 탈락 기준이 아니다.

## 개별 점검표

{md(selected[['error_rank','agbd','prediction','ndvi','evi','evi_denominator_rechecked','scl_r12p5_vegetation_fraction','scl_r50_vegetation_fraction','ndvi_r50_sd']])}

## 주변 피복 해석의 한계

반경 12.5m는 주변 민감도 점검을 위한 원형 근사이며 실제 GEDI footprint 모양·위치 오차·에너지 가중치를 복원한 것이 아니다. 원본 20m SCL의 픽셀 **중심**이 반경 안에 들어오는 값의 비율이며 면적 중첩 비율이 아니다. 50m도 진단용 거리이며 최적 거리로 선택하지 않았다. NDVI 주변 통계는 기존 품질 마스크를 통과한 저장 GeoTIFF의 10m 픽셀 중심 기준이며 유효 개수를 CSV에 기록했다.

SCL=4는 식생이고 산림·농경지 구분이 아니다. SCL=5는 비식생이며 세부 지목을 확정하지 않는다. 실제 산림·토지피복 분류 원자료가 로컬에 없어 `forest_status=not_assessed`를 유지한다. RGB는 0~0.3 고정 스트레치와 gamma 2.2를 사용한 육안 검토 보조 그림이며 표시에서 음수·상한을 자른 것은 원자료 수정이 아니다. 검은 픽셀만으로 수역·오류를 판정하지 않는다.

## 산출물과 다음 단계

- `all441_pixel_audit.csv`: 441점 전체 반사율·품질 재점검·주변 통계.
- `review21_pixel_audit.csv`: 상위 20점과 EVI 점검점의 합집합, 정확한 point_id·좌표 포함.
- `pixel_review_01.png`~`pixel_review_03.png`: 순위별 RGB·NDVI·SCL 상세도.
- `group_comparison.csv`, `summary.json`, `validation.json`: 집계와 검증 근거.

이번 점검에서 중심 픽셀의 명백한 처리 불일치는 확인되지 않았다. 주변 SCL 혼합 여부는 추가 피복 검토 대상이며 오류 확정 근거가 아니다. 다음은 기존 공간 분할·RF 설정을 고정하고 원본 밴드를 추가하는 비교 실험을 설계하는 것이다. 산림 여부는 별도 토지피복 자료 확보 후 판정한다.

## 재현

```powershell
& .\\week04\\.venv\\Scripts\\python.exe .\\week07\\step2\\inspect_pixels.py --output-dir .\\week07\\step2\\reproduction_outputs
```
'''
    (out / 'report.md').write_text(report, encoding='utf-8')
    for p, digest in hashes.items():
        assert sha(ROOT / p) == digest
    for p, stat in raw_stats.items():
        path = ROOT / p
        assert [path.stat().st_size, path.stat().st_mtime_ns] == stat
    validation = dict(status='passed', input_sha256_unchanged=hashes, raw_size_mtime_unchanged=raw_stats,
                      band_checks=checks, evi_formula_tolerance=dict(atol=.001, rtol=1e-5),
                      processed_reflectance_match=True, center_scl_match=True, center_classi_clear=True)
    for name, value in [('summary', summary), ('validation', validation)]:
        (out / f'{name}.json').write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
