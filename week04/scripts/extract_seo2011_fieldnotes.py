"""Local-only, loss-checked extraction of Seo (2011) Appendix 2 field notes.

No network, OCR guessing, species-name substitution, or biomass calculation.
Tree records retain duplicates and source spellings. BA uses pi*(DBH/200)^2.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import re
import sys
import xml.etree.ElementTree as ET

import pypdf

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / 'data/literature/서정범_2011_금정산_식물군집구조_동아대박사.pdf'
FIELD = ROOT / 'data/field'
OUT = ROOT / 'outputs'
TMP = ROOT / 'tmp/seo2011_pdf'
ZONES = ('중앙', '동측', '서측', '남측', '북측')
COLORS = {'중앙': '#0072B2', '동측': '#E69F00', '서측': '#009E73', '남측': '#CC79A7', '북측': '#6A3D9A'}
CAMPUS = (129.0738, 35.2299, 129.0847, 35.2386)
NUMBER = r'\d+(?:\.\d+)?'
TRIPLE = re.compile(rf'([^\s\d]+)\s+({NUMBER})\s+({NUMBER})(?=\s|$)')
SUSPECT_NAMES = {'개어서나무', '대팻집', '고뢰쇠나무'}
PLOT_COLUMNS = ['plot_id', 'zone', '위도', '경도', '조사지명', '구', '표고', '사면방위', '경사', '나무 수']
TREE_COLUMNS = ['plot_id', '수종', '흉고직경_cm', '수고_m']


def save_csv(path, records, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(path.suffix + '.tmp')
    with pending.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    pending.replace(path)


def search(pattern, text, label):
    match = re.search(pattern, text, re.S)
    if match is None:
        raise ValueError(f'Missing or unreadable {label}')
    return match


def numeric(pattern, text, label):
    return float(search(pattern, text, label).group(1))


def dms(text, hemisphere):
    pattern = rf'{hemisphere}\s*(\d+)\s*[°º]\s*(\d+)\s*[′\']\s*({NUMBER})\s*[″"]'
    match = search(pattern, text, f'{hemisphere} GPS')
    degrees, minutes, seconds = map(float, match.groups())
    if not 0 <= minutes < 60 or not 0 <= seconds < 60:
        raise ValueError(f'Invalid DMS: {match.group(0)}')
    return degrees + minutes / 60 + seconds / 3600, match.group(0)


def tree_block(text):
    return search(r'\n\s*매목\s*\n(.*?)\n\s*종조성\s*\n', text, 'tree table').group(1)


def parse_tree_table(text):
    records = []
    table_row = 0
    for line in tree_block(text).splitlines():
        if not line.strip() or '흉고직경' in line:
            continue
        matches = list(TRIPLE.finditer(line))
        if not 1 <= len(matches) <= 4 or TRIPLE.sub('', line).strip():
            raise ValueError(f'Unparsed tree row: {line!r}')
        table_row += 1
        for table_column, match in enumerate(matches, 1):
            records.append({'수종': match.group(1), '흉고직경_cm': float(match.group(2)),
                            '수고_m': float(match.group(3)), 'table_row': table_row,
                            'table_group': table_column, 'raw_cell': match.group(0)})
    if not records:
        raise ValueError('Empty tree table')
    return records


def parse_layers(text, plot_id, pdf_page):
    block = search(r'층별\s*식피율\s*및\s*수고범위(.*?)비\s*고\s*\(', text, 'layers').group(1)
    records = []
    for layer in ('교목층', '아교목층', '관목층', '초본층'):
        pattern = rf'(?m)^\s*{layer}\s*\(([^)]+)\)\s*({NUMBER})\s*%\s*({NUMBER})?\s*~\s*({NUMBER})?\s*m'
        match = search(pattern, block, layer)
        records.append({'plot_id': plot_id, '층': layer, '기호': match.group(1).replace(' ', ''),
                        '식피율_pct': float(match.group(2)),
                        '최소수고_m': float(match.group(3)) if match.group(3) else '',
                        '최대수고_m': float(match.group(4)) if match.group(4) else '', 'pdf_page': pdf_page})
    return records


def extract():
    reader = pypdf.PdfReader(PDF)
    plots, trees, provenance, sources, layers, issues, broken = [], [], [], [], [], [], []
    for pdf_page in range(147, len(reader.pages) + 1):
        page = reader.pages[pdf_page - 1]
        plain = page.extract_text()
        if not re.search(r'GPS\s*좌표', plain) or '매목' not in plain:
            continue
        layout = page.extract_text(extraction_mode='layout')
        (TMP / f'p{pdf_page:03}_plain.txt').write_text(plain, encoding='utf-8')
        (TMP / f'p{pdf_page:03}_layout.txt').write_text(layout, encoding='utf-8')
        plot_id = ''
        try:
            identity = search(r'(중앙|동측|서측|남측|북측)\s*[-－−]\s*(\d{1,2})', plain, 'plot ID')
            zone = identity.group(1)
            plot_id = f'{zone}-{int(identity.group(2)):02}'
            lat, lat_raw = dms(plain, 'N')
            lon, lon_raw = dms(plain, 'E')
            location = ' '.join(search(r'조사지명\s*:\s*(.*?)\s*도엽명', plain, 'location').group(1).split())
            district = search(r'(금정구|북구|양산시)', location, 'district').group(1)
            elevation = numeric(rf'표고\s*({NUMBER})\s*m', plain, 'elevation')
            aspect = numeric(rf'사면방위\s*({NUMBER})\s*°', plain, 'aspect')
            slope = numeric(rf'평균경사\s*({NUMBER})\s*°', plain, 'slope')
            area = numeric(rf'면\s*적\s*({NUMBER})\s*㎡', plain, 'area')
            table = parse_tree_table(plain)
            layout_table = parse_tree_table(layout)
            key = lambda r: (r['수종'], r['흉고직경_cm'], r['수고_m'])
            if [key(r) for r in table] != [key(r) for r in layout_table]:
                raise ValueError('Plain and layout tree extraction disagree')
            # Independently read location/coordinates from the layout version.
            if dms(layout, 'N')[0] != lat or dms(layout, 'E')[0] != lon:
                raise ValueError('Plain and layout GPS extraction disagree')
            if not 35 <= lat <= 36 or not 128 <= lon <= 130:
                raise ValueError('GPS outside the expected region')
            if area != 400 or not 0 <= aspect <= 360 or not 0 <= slope <= 90:
                raise ValueError('Unexpected area, aspect or slope')
            page_layers = parse_layers(plain, plot_id, pdf_page)
            plots.append({'plot_id': plot_id, 'zone': zone, '위도': lat, '경도': lon,
                          '조사지명': location, '구': district, '표고': elevation,
                          '사면방위': aspect, '경사': slope, '나무 수': len(table)})
            printed_page = int(search(r'^\s*(\d+)\s*\n', plain, 'printed page').group(1))
            sources.append({'plot_id': plot_id, 'pdf_page': pdf_page, 'printed_page': printed_page,
                            '위도_원문': lat_raw, '경도_원문': lon_raw, '면적_m2': area,
                            '매목_본수_plain': len(table), '매목_본수_layout': len(layout_table),
                            '추출상태': '일치', '측지계': '야장 좌표란에 미기재'})
            for tree in table:
                tree['plot_id'] = plot_id
                trees.append(tree)
                source = {'tree_csv_row': len(trees) + 1, 'plot_id': plot_id, 'pdf_page': pdf_page,
                          'printed_page': printed_page, 'table_row': tree['table_row'],
                          'table_group': tree['table_group'], 'raw_cell': tree['raw_cell']}
                provenance.append(source)
                if tree['수종'] in SUSPECT_NAMES:
                    issues.append(dict(source, issue='원문 수종 표기 확인',
                                       detail=f'{tree["수종"]}: 원문 그대로 보존, 정규화·통합하지 않음'))
            layers.extend(page_layers)
            print(f'PDF {pdf_page}, {plot_id}: {len(table)} records, GPS OK, plain/layout agree', flush=True)
        except (ValueError, IndexError) as exc:
            broken.append({'pdf_page': pdf_page, 'plot_id': plot_id, 'issue': str(exc)})
    return plots, trees, provenance, sources, layers, issues, broken, len(reader.pages)


def summarize(plots, trees):
    by_plot = defaultdict(list)
    for r in trees:
        by_plot[r['plot_id']].append(r)
    summaries = []
    for plot in plots:
        rows = by_plot[plot['plot_id']]
        counts, basal = Counter(), defaultdict(float)
        for tree in rows:
            counts[tree['수종']] += 1
            basal[tree['수종']] += math.pi * (tree['흉고직경_cm'] / 200) ** 2
        max_ba = max(basal.values())
        dominant = ' / '.join(sorted(name for name, value in basal.items()
                                    if math.isclose(value, max_ba, rel_tol=1e-12)))
        count_dominant = ' / '.join(sorted(name for name, count in counts.items() if count == max(counts.values())))
        summaries.append({'plot_id': plot['plot_id'], 'zone': plot['zone'], '면적_m2': 400,
                          '나무 수': len(rows), '1ha당_본수': len(rows) / .04,
                          '평균_흉고직경_cm': sum(t['흉고직경_cm'] for t in rows) / len(rows),
                          '평균_수고_m': sum(t['수고_m'] for t in rows) / len(rows),
                          '흉고단면적_합계_m2_ha': sum(basal.values()) / .04,
                          '우점_수종': dominant, '우점_기준': '흉고단면적 합계 최대(원문 수종 표기별)',
                          '우점종_흉고단면적_비율_pct': 100 * max_ba / sum(basal.values()),
                          '본수_기준_우점_수종': count_dominant})
    return summaries


def write_kml(plots, summaries, sources):
    ns = 'http://www.opengis.net/kml/2.2'
    ET.register_namespace('', ns)
    tag = lambda name: f'{{{ns}}}{name}'
    root = ET.Element(tag('kml'))
    doc = ET.SubElement(root, tag('Document'))
    ET.SubElement(doc, tag('name')).text = '서정범 2011 금정산 조사구 50개'
    ET.SubElement(doc, tag('description')).text = ('2009년 5월–2010년 9월 조사. 각 조사구 20×20m. '
        '좌표는 야장의 도분초를 십진수로 변환한 위치점이다. 야장 좌표란에는 측지계가 미기재되어 '
        '추가 측지계 변환 없이 KML에 표시했다. 외부 재배포용 자료가 아님.')
    for zone, color in COLORS.items():
        style = ET.SubElement(doc, tag('Style'), id=f'z{ZONES.index(zone)}')
        icon = ET.SubElement(style, tag('IconStyle'))
        color = color.lstrip('#')
        ET.SubElement(icon, tag('color')).text = 'ff' + color[4:6] + color[2:4] + color[0:2]
        ET.SubElement(icon, tag('scale')).text = '0.8'
    lookup = {r['plot_id']: r for r in summaries}
    pages = {r['plot_id']: r['pdf_page'] for r in sources}
    for plot in plots:
        r = lookup[plot['plot_id']]
        mark = ET.SubElement(doc, tag('Placemark'))
        ET.SubElement(mark, tag('name')).text = plot['plot_id']
        ET.SubElement(mark, tag('styleUrl')).text = f'#z{ZONES.index(plot["zone"])}'
        ET.SubElement(mark, tag('description')).text = (
            f'조사번호: {plot["plot_id"]}\n{plot["조사지명"]}\n'
            f'우점 수종(흉고단면적 기준): {r["우점_수종"]}\n'
            f'평균 흉고직경: {r["평균_흉고직경_cm"]:.2f} cm\n'
            f'평균 수고: {r["평균_수고_m"]:.2f} m\n'
            f'매목 기록: {plot["나무 수"]}건 / 400m²\n'
            f'흉고단면적: {r["흉고단면적_합계_m2_ha"]:.3f} m²/ha\n'
            f'원문 PDF {pages[plot["plot_id"]]}쪽; 2009–2010년 조사\n'
            '원문 수종 표기 보존. 위치점의 초 단위 GPS 정밀도 및 미기재 측지계에 유의.')
        point = ET.SubElement(mark, tag('Point'))
        ET.SubElement(point, tag('coordinates')).text = f'{plot["경도"]:.9f},{plot["위도"]:.9f},0'
    ET.ElementTree(root).write(OUT / 'seo2011_plots.kml', encoding='utf-8', xml_declaration=True)


def write_map(plots):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.transforms import Bbox
    from pyproj import Geod
    plt.rcParams['font.family'] = 'Malgun Gothic'
    plt.rcParams['axes.unicode_minus'] = False
    fig, ax = plt.subplots(figsize=(9, 11))
    for zone in ZONES:
        rows = [r for r in plots if r['zone'] == zone]
        ax.scatter([r['경도'] for r in rows], [r['위도'] for r in rows], s=42, c=COLORS[zone],
                   edgecolors='white', linewidths=.5, label=f'{zone} ({len(rows)}개)', zorder=3)
    w, s, e, n = CAMPUS
    ax.add_patch(Rectangle((w, s), e-w, n-s, facecolor='0.95', edgecolor='0.35', ls='--', lw=1))
    ax.plot(129.0784, 35.2343, marker='*', color='black', ms=12, label='부산대 부산캠퍼스(참고)', zorder=4)
    ax.annotate('부산대\n부산캠퍼스', (129.0784, 35.2343), xytext=(8, 5), textcoords='offset points', fontsize=9)
    xs, ys = [r['경도'] for r in plots] + [w, e], [r['위도'] for r in plots] + [s, n]
    ax.set(xlim=(min(xs)-.006, max(xs)+.014), ylim=(min(ys)-.006, max(ys)+.006),
           xlabel='경도 (십진수)', ylabel='위도 (십진수)',
           title='금정산 현장 조사구 50개\n서정범(2011) 부록 야장 · 2009–2010년 조사')
    ax.set_aspect(1 / math.cos(math.radians(35.26)))
    ax.ticklabel_format(axis='both', useOffset=False)
    ax.grid(alpha=.18)
    ax.legend(loc='upper right', fontsize=9, framealpha=.96)
    fig.canvas.draw()
    occupied = []
    factor = fig.dpi / 72
    for i, row in enumerate(plots):
        px, py = ax.transData.transform((row['경도'], row['위도']))
        offsets = [(6, 4), (-6, 4), (6, -6), (-6, -6), (12, 10), (-12, 10), (12, -12), (-12, -12),
                   (20, 0), (-20, 0), (0, 16), (0, -18), (25, 15), (-25, -15)]
        if i % 2:
            offsets = [(-a, b) for a, b in offsets]
        scored = []
        for dx, dy in offsets:
            left = px + dx*factor - (15 if dx < 0 else 0)
            bottom = py + dy*factor - (11 if dy < 0 else 0)
            box = Bbox.from_bounds(left-2, bottom-2, 19, 15)
            overlap = sum(box.overlaps(other) for other in occupied)
            scored.append((overlap*1000 + math.hypot(dx, dy), dx, dy, box))
        _, dx, dy, box = min(scored, key=lambda x: x[0])
        occupied.append(box)
        ax.annotate(row['plot_id'].split('-')[1], (row['경도'], row['위도']), xytext=(dx, dy),
                    textcoords='offset points', fontsize=8, color=COLORS[row['zone']],
                    ha='right' if dx < 0 else 'left', va='top' if dy < 0 else 'bottom',
                    bbox=dict(boxstyle='round,pad=.08', fc='white', ec='none', alpha=.85),
                    arrowprops=dict(arrowstyle='-', color=COLORS[row['zone']], lw=.45) if abs(dx)+abs(dy)>20 else None)
    geod = Geod(ellps='WGS84')
    start_lon, start_lat = min(xs)-.002, min(ys)-.0035
    end_lon, end_lat, _ = geod.fwd(start_lon, start_lat, 90, 1000)
    ax.plot([start_lon, end_lon], [start_lat, end_lat], color='black', lw=2)
    ax.text((start_lon+end_lon)/2, start_lat+.0007, '1 km', ha='center', fontsize=8)
    fig.text(.5, .035, '색은 구역, 숫자는 조사번호 끝자리(01–10) · 좌표는 원문 GPS의 도분초 변환값',
             ha='center', fontsize=9)
    fig.savefig(OUT / 'seo2011_plots_map.png', dpi=190, bbox_inches='tight')
    plt.close(fig)


def write_report(validation, plots, trees, issues):
    lines = ['# 서정범(2011) 금정산 부록 야장 데이터화', '',
             '외부 업로드 없이 로컬 파일만 생성했다. 바이오매스·탄소 및 상대생장식 계산은 수행하지 않았다.', '',
             '## 검증 결과', '',
             f'- 조사구: {len(plots)}개, 구역별 10개씩', f'- GPS: {validation["gps_count"]}/50개',
             f'- 매목 기록: {len(trees):,}건 (원문 중복 기록 및 등나무 포함)',
             f'- 이상 수치 행: {validation["abnormal_tree_count"]}건 '
             '(DBH ≤0 또는 >100cm, 수고 ≤0 또는 >40m, 비유한 값)',
             f'- 추출 실패 페이지: {validation["broken_pages"] or "없음"}',
             f'- plain/layout 매목값 일치: {validation["plain_layout_agree_plots"]}/50개',
             '- 출처: PDF 147–196쪽, 논문 인쇄 쪽수 135–184쪽. 조사구당 1쪽.',
             '- 모든 조사구 면적은 원문에서 400m²로 확인. 조사 기간은 사용자 제공 2009년 5월–2010년 9월.', '',
             '| 구역 | 조사구 | 매목 기록 |', '|---|---:|---:|']
    for zone in ZONES:
        ids = {p['plot_id'] for p in plots if p['zone'] == zone}
        lines.append(f'| {zone} | {len(ids)} | {sum(t["plot_id"] in ids for t in trees)} |')
    lines += ['', '## 계산 정의', '',
              '- 1ha당 본수 = 매목 기록 수 / 0.04ha.',
              '- 평균 흉고직경과 수고 = 매목 표의 산술평균. 동일한 수종·측정값의 반복행도 별개 개체 기록으로 유지.',
              '- 흉고단면적(m²/ha) = Σ[π × (DBH(cm)/200)²] / 0.04.',
              '- 우점 수종 = 원문 수종 표기별 흉고단면적 합계 최대인 수종. 동률이면 모두 기록.',
              '- 본수 기준 우점 수종도 별도 열로 제공. 원문의 종조성 D.S를 이용한 우점 판정과는 구분.',
              '- 수종 표기를 자동 통합하지 않았으므로 수종별 분석 전 확인이 필요하다.', '',
              '## 원문 표기 확인 목록', '',
              '다음 표기는 PDF 이미지에서도 확인했으며, 추출 오류로 간주하거나 임의 교정하지 않았다.', '',
              '| 조사구 | PDF 쪽 | 원문 표기 | 매목 행 위치(행, 묶음) |', '|---|---:|---|---|']
    for issue in issues:
        lines.append(f'| {issue["plot_id"]} | {issue["pdf_page"]} | {issue["raw_cell"]} | '
                     f'{issue["table_row"]}, {issue["table_group"]} |')
    lines += ['', '## 좌표와 파일 해석', '',
              '- 위도·경도는 D + M/60 + S/3600으로 변환. 초 단위 정밀도를 넘어선 추가 위치 정확도를 뜻하지 않는다.',
              '- 야장 GPS 좌표란에는 측지계가 명시되지 않았다. KML에는 원문 숫자를 추가 측지계 변환 없이 넣었다.',
              '- 구 필드는 조사지명의 금정구/북구/양산시를 그대로 추출했으며, 행정경계 공간분석 결과가 아니다.',
              '- 지도에서 부산대 캠퍼스는 이전 작업의 참고 사각 범위이며 논문 조사구 경계가 아니다.',
              '- 층별 식피율·수고 범위도 별도 CSV로 보존했다. 빈 수고 범위는 빈칸으로 유지했다.',
              '- CSV는 UTF-8 BOM 형식. main trees CSV의 행 번호는 provenance 파일의 tree_csv_row로 원문에 연결된다.', '',
              '## 생성 파일', '']
    lines += [f'- data/field/{p.name}' for p in sorted(FIELD.glob('seo2011_*')) if p.is_file()]
    lines += ['- outputs/seo2011_plots_map.png', '- outputs/seo2011_plots.kml', '- outputs/seo2011_report.md']
    (OUT / 'seo2011_report.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    FIELD.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    plots, trees, provenance, sources, layers, issues, broken, pdf_pages = extract()
    expected = {f'{zone}-{i:02}' for zone in ZONES for i in range(1, 11)}
    actual = [r['plot_id'] for r in plots]
    if len(actual) != 50 or set(actual) != expected or broken:
        (OUT / 'seo2011_failed_extraction.json').write_text(json.dumps(
            {'broken': broken, 'missing_ids': sorted(expected-set(actual)), 'parsed_plot_count': len(actual)},
            ensure_ascii=False, indent=2), encoding='utf-8')
        raise RuntimeError('Incomplete extraction; no final data tables generated. See failure report.')
    abnormal = []
    for row, source in zip(trees, provenance):
        reasons = []
        for key, maximum in [('흉고직경_cm', 100), ('수고_m', 40)]:
            value = row[key]
            if not math.isfinite(value) or value <= 0 or value > maximum:
                reasons.append(key)
        if reasons:
            abnormal.append(dict(row, **source, reason='; '.join(reasons)))
    summaries = summarize(plots, trees)
    validation = {'source_pdf': str(PDF.relative_to(ROOT)), 'source_sha256': hashlib.sha256(PDF.read_bytes()).hexdigest(),
                  'pdf_pages': pdf_pages, 'source_pdf_pages': [r['pdf_page'] for r in sources],
                  'plot_count': len(plots), 'gps_count': sum(math.isfinite(r['위도']) and math.isfinite(r['경도'])
                                                           for r in plots),
                  'zone_plot_counts': dict(Counter(p['zone'] for p in plots)),
                  'zone_tree_counts': {z: sum(r['plot_id'].startswith(z+'-') for r in trees) for z in ZONES},
                  'district_plot_counts': dict(Counter(p['구'] for p in plots)), 'tree_count': len(trees),
                  'abnormal_tree_count': len(abnormal), 'broken_pages': broken,
                  'plain_layout_agree_plots': len(sources), 'layer_count': len(layers),
                  'species_spelling_note_rows': len(issues), 'raw_species_label_count': len({r['수종'] for r in trees}),
                  'dbh_cm_range': [min(r['흉고직경_cm'] for r in trees), max(r['흉고직경_cm'] for r in trees)],
                  'height_m_range': [min(r['수고_m'] for r in trees), max(r['수고_m'] for r in trees)],
                  'liana_label_records': sum(r['수종'] == '등나무' for r in trees),
                  'duplicate_coordinate_count': len(plots)-len({(r['위도'],r['경도']) for r in plots}),
                  'extraction_library': f'pypdf {pypdf.__version__}',
                  'external_uploads': False, 'biomass_or_carbon_calculated': False}
    assert len(trees) == 1750, 'Unexpected total: review source before assuming the expected total.'
    save_csv(FIELD / 'seo2011_plots.csv', plots, PLOT_COLUMNS)
    save_csv(FIELD / 'seo2011_trees.csv', trees, TREE_COLUMNS)
    save_csv(FIELD / 'seo2011_plot_summary.csv', summaries, list(summaries[0]))
    save_csv(FIELD / 'seo2011_layers.csv', layers, list(layers[0]))
    save_csv(FIELD / 'seo2011_source_index.csv', sources, list(sources[0]))
    save_csv(FIELD / 'seo2011_tree_provenance.csv', provenance, list(provenance[0]))
    save_csv(FIELD / 'seo2011_abnormal_trees.csv', abnormal,
             ['tree_csv_row', 'plot_id', '수종', '흉고직경_cm', '수고_m', 'pdf_page', 'table_row', 'table_group', 'reason'])
    save_csv(FIELD / 'seo2011_extraction_notes.csv', issues, list(issues[0]))
    (FIELD / 'seo2011_validation.json').write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding='utf-8')
    write_kml(plots, summaries, sources)
    write_map(plots)
    write_report(validation, plots, trees, issues)
    assert not {'pandas', 'geopandas', 'rasterio', 'pyogrio'} & set(sys.modules)
    print(json.dumps(validation, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
