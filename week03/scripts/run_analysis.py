"""Hongcheon GEDI L4A V2.1 analysis; run from any working directory."""
from __future__ import annotations

import argparse
import json
import netrc
import os
from pathlib import Path
import sys
import urllib.request
import urllib.parse

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')

ROOT = Path(__file__).resolve().parents[1]
BBOX = (127.99, 37.66, 128.21, 37.84)
SHORT_NAME = 'GEDI_L4A_AGB_Density_V2_1_2056'
TEMPORAL = ('2021-01-01T00:00:00Z', '2021-12-31T23:59:59Z')
BEAMS = ('BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011',
         'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011')
VARIABLES = ('shot_number', 'lat_lowestmode', 'lon_lowestmode', 'agbd',
             'agbd_se', 'l4_quality_flag', 'l2_quality_flag', 'degrade_flag',
             'sensitivity', 'solar_elevation', 'delta_time')
WC_NAME = 'ESA_WorldCover_10m_2021_v200_N36E126_Map.tif'
WC_URL = 'https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/' + WC_NAME
AUTH_HELP = '''Earthdata 인증 정보가 없어 중단했습니다. 비밀번호를 채팅이나 코드에 넣지 마세요.
방법 1: 본인 PowerShell에서 아래를 실행하고, 같은 창에서 분석을 실행하세요.
  $env:EARTHDATA_USERNAME = Read-Host 'Earthdata username'
  $secret = Read-Host 'Earthdata password' -AsSecureString
  $env:EARTHDATA_PASSWORD = [System.Net.NetworkCredential]::new('', $secret).Password
  .\\.venv\\Scripts\\python.exe scripts\\run_analysis.py
방법 2: 사용자 홈의 _netrc 파일에 다음 항목을 직접 입력하세요.
  machine urs.earthdata.nasa.gov
  login YOUR_USERNAME
  password YOUR_PASSWORD
파일은 평문이므로 본인만 접근 가능한 위치에 보관하세요. 저장 후 다시 실행하세요.
계정 생성: https://urs.earthdata.nasa.gov/users/new
'''


def credential_strategy():
    if os.getenv('EARTHDATA_USERNAME') and os.getenv('EARTHDATA_PASSWORD'):
        return 'environment'
    candidates = [Path(os.environ['NETRC'])] if os.getenv('NETRC') else []
    candidates += [Path.home() / '_netrc', Path.home() / '.netrc']
    for path in candidates:
        if not path.is_file():
            continue
        try:
            entry = netrc.netrc(str(path)).authenticators('urs.earthdata.nasa.gov')
            if entry and entry[0] and entry[2]:
                os.environ['NETRC'] = str(path)
                return 'netrc'
        except (OSError, netrc.NetrcParseError):
            continue
    return None


def read_bbox(paths):
    import h5py
    import numpy as np
    import pandas as pd
    frames = []
    beam_log = []
    west, south, east, north = BBOX
    for path in sorted(paths):
        print(f'Extract: {path.name}', flush=True)
        with h5py.File(path, 'r') as file:
            for beam in BEAMS:
                if beam not in file:
                    beam_log.append({'file': path.name, 'beam': beam, 'status': 'absent'})
                    continue
                group = file[beam]
                missing = [v for v in VARIABLES if v not in group]
                if missing:
                    raise ValueError(f'{path.name}/{beam}: missing {missing}')
                lat, lon = group['lat_lowestmode'][:], group['lon_lowestmode'][:]
                idx = np.flatnonzero((lat >= south) & (lat <= north) &
                                     (lon >= west) & (lon <= east))
                beam_log.append({'file': path.name, 'beam': beam,
                                 'status': 'read', 'bbox_count': int(len(idx))})
                if not len(idx):
                    continue
                columns = {v: group[v][idx] for v in VARIABLES}
                if not np.issubdtype(columns['shot_number'].dtype, np.integer):
                    raise ValueError('shot_number must be stored as an integer')
                columns['shot_number'] = columns['shot_number'].astype('uint64')
                agb = columns['agbd'].astype('float64')
                valid = np.isfinite(agb) & (agb >= 0)
                dataset = group['agbd']
                for key in ('_FillValue', 'missing_value'):
                    if key in dataset.attrs:
                        for fill in np.atleast_1d(dataset.attrs[key]):
                            valid &= agb != fill
                columns['agbd_valid'] = valid
                frame = pd.DataFrame(columns)
                frame['beam'] = beam
                frame['source_file'] = path.name
                frames.append(frame)
    if frames:
        result = pd.concat(frames, ignore_index=True)
    else:
        result = pd.DataFrame({v: pd.Series(dtype='uint64' if v == 'shot_number' else 'float64')
                               for v in VARIABLES})
        result['agbd_valid'] = pd.Series(dtype=bool)
        result['beam'] = pd.Series(dtype=str)
        result['source_file'] = pd.Series(dtype=str)
    result = result.sort_values(['source_file', 'beam', 'shot_number']).reset_index(drop=True)
    result['observation_time_utc'] = pd.Timestamp('2018-01-01', tz='UTC') + pd.to_timedelta(
        result['delta_time'], unit='s')
    return result, beam_log


def quality_filter(all_points):
    import pandas as pd
    a = pd.Series(True, index=all_points.index)
    b = a & all_points['agbd_valid']
    c = b & all_points['l4_quality_flag'].eq(1)
    d = c & all_points['degrade_flag'].eq(0)
    e = d & all_points['solar_elevation'].lt(0)
    rows, previous = [], len(all_points)
    for stage, mask in [('a. 전체 점', a), ('b. AGBD 유효값', b),
                        ('c. L4 품질 = 1', c), ('d. degrade = 0 (최종)', d),
                        ('e. 야간 (참고용)', e)]:
        count = int(mask.sum())
        rows.append({'stage': stage, 'count': count,
                     'percent_of_all': 100 * count / len(a) if len(a) else 0,
                     'percent_of_previous': 100 * count / previous if previous else 0})
        previous = count
    return all_points.loc[d].copy(), pd.DataFrame(rows)


def valid_h5(path):
    import h5py
    try:
        with h5py.File(path, 'r') as f:
            beams = [b for b in BEAMS if b in f]
            if not beams:
                return False
            for b in beams:
                for name in VARIABLES:
                    ds = f[b][name]
                    if ds.size:
                        ds[-1]  # Detect truncated downloads before accepting a cached file.
        return True
    except (OSError, KeyError, ValueError):
        return False


def download_gedi(auth_strategy, raw):
    import earthaccess
    auth = earthaccess.login(strategy=auth_strategy, persist=False)
    if not auth.authenticated:
        raise RuntimeError('Earthdata authentication failed')
    granules = earthaccess.search_data(short_name=SHORT_NAME,
                                      bounding_box=BBOX, temporal=TEMPORAL, count=-1)
    if not granules:
        raise RuntimeError('CMR returned no granules; no analysis outputs were written')
    print(f'CMR granules: {len(granules)} (user expected: 14)', flush=True)
    manifest, paths = [], []
    for granule in granules:
        links = granule.data_links()
        names = {Path(urllib.parse.unquote(urllib.parse.urlparse(link).path)).name
                 for link in links if urllib.parse.urlparse(link).path.lower().endswith('.h5')}
        if len(names) != 1:
            raise ValueError(f'Expected one HDF5 filename, found {sorted(names)}')
        name = next(iter(names))
        destination = raw / name
        if valid_h5(destination):
            print(f'Cached: {name}', flush=True)
        else:
            staging = raw / '_incomplete' / name.removesuffix('.h5')
            staging.mkdir(parents=True, exist_ok=True)
            partial = staging / name
            if partial.exists() and not valid_h5(partial):
                partial.unlink()  # Only this explicitly named incomplete download.
            downloaded = earthaccess.download([granule], local_path=str(staging), threads=1)
            if not downloaded or not valid_h5(partial):
                raise RuntimeError(f'Incomplete or invalid HDF5 download: {name}')
            partial.replace(destination)
        paths.append(destination)
        manifest.append({'file': name, 'bytes': destination.stat().st_size,
                         'cmr': dict(granule)})
    (raw / 'download_manifest.json').write_text(json.dumps({
        'short_name': SHORT_NAME, 'bbox': BBOX, 'temporal': TEMPORAL,
        'granule_count': len(granules), 'granules': manifest}, ensure_ascii=False, indent=2), encoding='utf-8')
    return sorted(set(paths)), len(granules)


def worldcover_values(points, folder):
    import numpy as np
    import rasterio
    from rasterio.warp import transform
    path = folder / WC_NAME

    def validate(candidate):
        try:
            with rasterio.open(candidate) as ds:
                ds.read(1, window=((0, 1), (0, 1)))
                return ds.crs is not None
        except (OSError, rasterio.errors.RasterioError):
            return False

    if not validate(path):
        partial = path.with_suffix('.tif.part')
        print('Downloading ESA WorldCover 2021 tile...', flush=True)
        with urllib.request.urlopen(WC_URL, timeout=180) as response, partial.open('wb') as out:
            expected = response.headers.get('Content-Length')
            total = 0
            while chunk := response.read(8 * 1024 * 1024):
                out.write(chunk)
                total += len(chunk)
        if expected and total != int(expected):
            raise RuntimeError('Incomplete WorldCover download')
        if not validate(partial):
            raise RuntimeError('Invalid WorldCover raster')
        partial.replace(path)
    values = []
    with rasterio.open(path) as ds:
        xs, ys = transform('EPSG:4326', ds.crs,
                           points['lon_lowestmode'].tolist(), points['lat_lowestmode'].tolist()) if len(points) else ([], [])
        for x, y in zip(xs, ys):
            if not (ds.bounds.left <= x < ds.bounds.right and ds.bounds.bottom < y <= ds.bounds.top):
                raise ValueError('GEDI point outside WorldCover tile')
        for sample in ds.sample(zip(xs, ys), indexes=1, masked=True):
            values.append(None if np.ma.is_masked(sample[0]) else int(sample[0]))
    import pandas as pd
    points['worldcover_class'] = pd.array(values, dtype='Int64')
    points['tree_cover'] = points['worldcover_class'].eq(10).fillna(False)
    return points


def save_results(all_points, filtered, summary, processed, figures, granule_count):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import numpy as np
    import pandas as pd
    all_points.to_csv(processed / 'gedi_l4a_2021_all.csv', index=False, encoding='utf-8-sig')
    filtered.to_csv(processed / 'gedi_l4a_2021_filtered.csv', index=False, encoding='utf-8-sig')
    summary.to_csv(processed / 'filter_summary.csv', index=False, encoding='utf-8-sig')
    west, south, east, north = BBOX
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(all_points.lon_lowestmode, all_points.lat_lowestmode, s=6, c='0.7', label='All points')
    ax.scatter(filtered.lon_lowestmode, filtered.lat_lowestmode, s=8, c='forestgreen', label='Quality filtered')
    ax.add_patch(Rectangle((west, south), east-west, north-south, fill=False, color='black'))
    ax.set(xlabel='Longitude (WGS84)', ylabel='Latitude (WGS84)', title='Hongcheon GEDI L4A V2.1 | 2021',
           xlim=(west-.003, east+.003), ylim=(south-.003, north+.003))
    ax.set_aspect(1 / np.cos(np.deg2rad((south+north)/2)))
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / 'map_points.png', dpi=200)
    plt.close(fig)
    tree = filtered.loc[filtered.tree_cover]
    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.histogram_bin_edges(filtered.agbd, bins=30) if len(filtered) else np.linspace(0, 1, 31)
    ax.hist(filtered.agbd, bins=bins, alpha=.45, color='steelblue', label='All quality-filtered points')
    ax.hist(tree.agbd, bins=bins, histtype='step', linewidth=2, color='forestgreen', label='Tree cover (WorldCover = 10)')
    ax.set(xlabel='AGBD (Mg/ha)', ylabel='Point count', title='AGBD distribution | 2021')
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / 'agbd_hist.png', dpi=200)
    plt.close(fig)

    def table(headers, rows):
        return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                          '| ' + ' | '.join(['---']*len(headers)) + ' |'] +
                         ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])

    def number(value):
        return f'{value:.3f}' if pd.notna(value) else '해당 없음'

    counts = table(['단계', '남은 점', '전체 대비(%)', '이전 단계 대비(%)'],
                   [[r.stage, r.count, f'{r.percent_of_all:.2f}', f'{r.percent_of_previous:.2f}']
                    for r in summary.itertuples()])
    stat_rows = []
    for label, frame in [('품질 필터 통과', filtered), ('그중 나무 피복', tree)]:
        series = frame.agbd
        se = frame.agbd_se.where(np.isfinite(frame.agbd_se) & (frame.agbd_se >= 0))
        stat_rows.append([label, len(frame)] + [number(v) for v in [series.mean(), series.median(),
                         series.min(), series.max(), series.std(ddof=1), se.mean()]])
    stats = table(['구분', '개수', '평균', '중앙값', '최소', '최대', '표준편차', 'agbd_se 평균'], stat_rows)
    months = pd.period_range('2021-01', '2021-12', freq='M').astype(str)
    monthly = pd.DataFrame(index=months)
    for label, frame in [('전체', all_points), ('필터 통과', filtered), ('나무 피복', tree)]:
        monthly[label] = frame.observation_time_utc.dt.strftime('%Y-%m').value_counts().reindex(months, fill_value=0)
    monthly.index.name = 'month_utc'
    monthly.to_csv(processed / 'monthly_counts.csv', encoding='utf-8-sig')
    monthly_table = table(['월(UTC)', *monthly.columns], monthly.reset_index().values.tolist())
    duplicate_count = int(all_points.shot_number.duplicated().sum())
    warning = '\n**점 부족 — 기간을 2021~2022로 넓히는 것 검토**\n' if len(filtered) < 300 else ''
    report = f'''# 홍천 GEDI L4A 2021 분석 결과

- 데이터: GEDI L4A Version 2.1 (`{SHORT_NAME}`); V3 사용 안 함.
- 구역: 경도 {west}~{east}, 위도 {south}~{north} (WGS84); 산림만으로 확정된 경계가 아닌 사각 구역.
- 조회 기간: 2021-01-01~2021-12-31 (UTC).
- CMR 검색 궤도 파일: {granule_count}개 (사전 예상 14개와 다를 수 있음).
- WorldCover: 2021 v200, N36E126; 점 중심 위치의 최근접 픽셀 분류값 사용.
- 나무 피복 점: {len(tree)}개 / 최종 {len(filtered)}개 ({100*len(tree)/len(filtered) if len(filtered) else 0:.2f}%).
- 중복 shot_number: {duplicate_count}개. 임의 제거 없이 기록 단위로 집계함.
- 고유 shot_number: {all_points.shot_number.nunique()}개; CSV에는 정수 원문 저장.
{warning}
## 단계별 필터 결과

{counts}

a~d는 누적 필터이며 d가 최종 결과다. e는 d 중 야간 관측만 집계한 참고값이다.
AGBD는 유한한 0 이상 값이며 원본 _FillValue/missing_value를 제외했다.
l2_quality_flag와 sensitivity는 보존하며 별도 추가 필터로 사용하지 않았다.
전체 대비 비율의 분모는 bbox 안 전체 점이다. 전체 점이 0개면 비율을 0으로 표시한다.

## AGBD 통계

AGBD와 agbd_se 단위: Mg/ha. 표준편차는 표본 표준편차(ddof=1).
agbd_se 평균에서는 음수·비유한 값을 제외했다. 단일 관측점의 표준편차는 해당 없음이다.

{stats}

## 월별 관측 개수

delta_time을 2018-01-01 00:00:00 UTC 기준으로 환산했다.

{monthly_table}

## 해석과 한계

- GEDI는 궤도를 따라 점으로 관측한다. 지도에서 궤도 사이 빈 구역과 편중 여부를 확인해야 하며, 점 개수만으로 구역 전체의 대표성을 보장할 수 없다.
- 산지 경사는 파형의 지면 판별과 바이오매스 추정에 영향을 줄 가능성이 있다. 이번에는 경사 자료를 결합하지 않아 그 영향을 직접 검증하지 않았다.
- GEDI AGBD는 현장에서 직접 측정한 정답이 아니라 모델 추정치다. 현장자료 없이 실제 정확도를 검증한 것으로 해석할 수 없다.
- agbd_se는 제품이 제공하는 추정오차 정보로, 모든 오차 원인을 포괄하는 보증값은 아니다.
- WorldCover 10은 나무 피복 분류다. 법적 산림 경계나 footprint 전체가 나무라는 뜻은 아니며 분류오차와 경계의 혼합 픽셀 문제가 있다.
- 야간 비교는 품질 통과 점의 부분집합을 집계한 것이며 정확도 향상을 입증하는 실험은 아니다.
- 궤도 파일의 공간적 교차가 bbox 내부 관측점 존재를 보장하지 않는다. 실제 추출 점을 기준으로 판단한다.
- CSV를 Excel에서 열면 긴 shot_number를 반올림할 수 있으므로 해당 열은 텍스트로 가져오거나 pandas의 uint64 dtype을 명시해야 한다.

## 출처

- https://doi.org/10.3334/ORNLDAAC/2056
- https://ornldaac.github.io/gedi_tutorials/notebooks/gedi-l4a-subsets/
- https://esa-worldcover.org/en/data-access
'''
    (ROOT / 'summary.md').write_text(report, encoding='utf-8')
    print(counts)
    print(f'Tree cover: {len(tree)} / {len(filtered)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true', help='Use existing HDF5 files only; WorldCover still needs a cached tile')
    args = parser.parse_args()
    raw = ROOT / 'data/raw/gedi_l4a_2021'
    processed = ROOT / 'data/processed'
    wc = ROOT / 'data/raw/worldcover'
    figures = ROOT / 'figures'
    for folder in (raw, processed, wc, figures):
        folder.mkdir(parents=True, exist_ok=True)
    if not args.offline:
        strategy = credential_strategy()
        if strategy is None:
            print(AUTH_HELP)
            return 2
        paths, count = download_gedi(strategy, raw)
    else:
        if not (wc / WC_NAME).is_file():
            raise RuntimeError('Offline mode requires the cached WorldCover tile')
        paths = sorted(raw.glob('GEDI04_A*.h5'))
        if not paths:
            raise RuntimeError('No cached HDF5 files')
        count = len(paths)
    all_points, beam_log = read_bbox(paths)
    # Granules overlap dates: enforce shot-level year as well as CMR query dates.
    all_points = all_points.loc[all_points.observation_time_utc.between(
        '2021-01-01T00:00:00Z', '2022-01-01T00:00:00Z', inclusive='left')].copy()
    filtered, summary = quality_filter(all_points)
    # Preserve successful extraction even if the independent raster download fails.
    all_points.to_csv(processed / 'gedi_l4a_2021_all.csv', index=False, encoding='utf-8-sig')
    filtered.to_csv(processed / 'gedi_l4a_2021_filtered.csv', index=False, encoding='utf-8-sig')
    summary.to_csv(processed / 'filter_summary.csv', index=False, encoding='utf-8-sig')
    (processed / 'beam_read_log.json').write_text(json.dumps(beam_log, ensure_ascii=False, indent=2), encoding='utf-8')
    filtered = worldcover_values(filtered, wc)
    save_results(all_points, filtered, summary, processed, figures, count)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'분석 중단 ({type(exc).__name__}). 인증 정보 노출 방지를 위해 서버 응답은 출력하지 않습니다.', file=sys.stderr)
        print('다운로드·네트워크·원본 파일 구조 및 실행 환경을 확인한 뒤 다시 실행하세요.', file=sys.stderr)
        sys.exit(1)
