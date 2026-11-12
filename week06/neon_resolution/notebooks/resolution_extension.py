"""Standard-library-only notebook authoring; never executes scientific cells."""
import ast
import hashlib
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _source(cell):
    return ''.join(cell['source'])


def _set(cell, source):
    cell['source'] = source.strip().splitlines(True)


def _guard_legacy(source):
    """Keep original functions verbatim; skip all old work in resume-only mode."""
    lines = source.splitlines(True)
    definitions = []; kept_lines = set()
    for node in ast.parse(source).body:
        segment = ''.join(lines[node.lineno - 1:node.end_lineno])
        keep = isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Import, ast.ImportFrom))
        if isinstance(node, ast.Assign):
            keep |= any(isinstance(t, ast.Name) and t.id in ('DET_COLUMNS', 'default') for t in node.targets)
        if keep:
            definitions.append(segment.rstrip())
            kept_lines.update(range(node.lineno - 1, node.end_lineno))
    head = '\n\n'.join(definitions)
    # Keep execution text verbatim, including semicolon-separated statements.
    # Extracting each AST node separately would duplicate those source lines.
    tail = ''.join(line for i, line in enumerate(lines) if i not in kept_lines)
    return head + ('\n\nif not RESOLUTION_ONLY:\n' + textwrap.indent(tail, '    ') if tail else '') + '\n'


def apply_resolution_extension(cells):
    assert len(cells) == 46, 'Apply extension to the unmodified v3 builder cells exactly once.'
    cells[0]['source'] = ('''# NEON GRSM·ORNL 2021: 항공사진 해상도 실험 추가

**Google Colab / T4 GPU. 기본은 `RESOLUTION_ONLY=True`: 지난주 3차 완료 결과만 불러옵니다.**
상단에서 하단까지 실행하면 기존 자료 다운로드·10cm baseline/A/B/C 추론을 건너뛰고,
Drive에 있는 원본 RGB/CHM crop으로25cm·51cm RGB 실험만 추가합니다.
같은23개 조사구(GRSM15+ORNL8), 같은 매칭·공식·CHM을 유지합니다.
첫 설치 후 런타임 재시작 안내가 나오면 재시작하고 Drive 설정 셀부터 실행하세요.

`RESOLUTION_BASE_DIR=None`은 완료된3차 결과 폴더나 Drive의 `outputs/neon_results.zip`을 자동 탐색합니다.
없거나 서로 다른 결과가 여러 개면 중단하며, 새 영상 다운로드나 기존 방법 재계산으로 대체하지 않습니다.
필요하면02번 셀의 `RESOLUTION_BASE_DIR`에3차 결과 폴더 경로를 넣으세요.
새 결과는 별도 `outputs/resolution_…/`, RGB·탐지 캐시는 `data/cache/resolution/…/`에 저장합니다.
HF 캐시는 `/content/hf_cache`, 마지막 다운로드는 `/content/neon_results.zip`입니다.

원본 설명과 기존 분석 코드는 아래에 보존했으며 기본 실행에서는 정의와 import만 읽습니다.
''').splitlines(True) + ['\n'] + cells[0]['source']
    _set(cells[2], _source(cells[2]).replace('"deepforest>=1.5,<3"', '"deepforest==2.1.0"').replace('scikit-learn ', '"scikit-learn==1.6.1" '))
    setup = _source(cells[4])
    setup = setup.replace("CFG = dict(", "RESOLUTION_ONLY = True  # 기존3차 자료만 재사용; False로 바꾸지 마세요.\n"
        "RESOLUTION_BASE_DIR = None  # 예: ROOT/'outputs/f55cc4636e94'\nCFG = dict(", 1)
    setup = setup.replace("OUT = ROOT/'outputs'/RUN_ID", "OUT = ROOT/'outputs'/('resolution_setup' if RESOLUTION_ONLY else RUN_ID)")
    begin = setup.index('try:\n    from google.colab import userdata')
    end = setup.index("print('결과 폴더:'", begin)
    token = setup[begin:end]
    setup = setup[:begin] + "TOKEN = ''  # 해상도 실험에는 NEON 토큰/영상 다운로드 불필요\nif not RESOLUTION_ONLY:\n" + textwrap.indent(token, '    ') + setup[end:]
    _set(cells[4], setup)
    _set(cells[3], _source(cells[3]) + '\n\n**해상도 추가 실행:** RESOLUTION_ONLY=True를 유지하세요. '
         'NEON 토큰은 묻지 않습니다. CFG의 .15/400은 보존된 baseline용이며, 새 고정 실험은 별도 .2/600입니다.')
    # Keep common imports/utilities and the method definitions; disable old actions.
    for index in (8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 34, 36, 38, 40):
        _set(cells[index], _guard_legacy(_source(cells[index])))
        cells[index - 1]['source'] += ['\n\n해상도 실험 기본 실행에서는 기존 처리·다운로드를 건너뜁니다. 함수 정의만 재사용합니다.\n']
    runtime = (HERE / 'resolution_experiment.py').read_text(encoding='utf-8')
    runtime_hash = hashlib.sha256(runtime.encode()).hexdigest()
    def md(text):
        return dict(cell_type='markdown', metadata={}, source=text.strip().splitlines(True))
    def code(text):
        return dict(cell_type='code', metadata={}, execution_count=None, outputs=[], source=text.strip().splitlines(True))
    extra = [
        md('''## 17A. 해상도 실험 함수와 캐시 재개

3차 결과와 좌표 없는 나무를 포함한 참조표, 조사구 경계, 기존RGB·CHM crop을 불러옵니다.
원본10cm의 네 방법은 캐시 파일을 그대로 재사용합니다. 25/51cm는 면적평균 축소→bilinear 확대하며,
확대본의 CRS·범위·transform·shape는 원본과 정확히 같습니다. 정수 픽셀 크기에 따른 실제GSD를 출력합니다.
모델 입력은 항상10cm 격자이므로 patch600의 지상 범위는60m로 유지됩니다.
'''),
        code(f"RESOLUTION_CODE_SHA = {runtime_hash!r}\n" + runtime),
        md('''## 17B. 10·25·51cm × 고정·nested OOF 비교 실행

25/51cm에서만 추론합니다. patch600을 고정하고 score.1/.2/.3/.4/.5를 같은GroupKFold로 선택합니다.
10cm 고정 결과의F1·탄소 비율을 기존A와 자동 비교하고, 다르면 경고합니다.
실제 수치·18행 사이트 비교표·138행 조사구 비교표·추세도·요청한3개 나란히 그림을 저장합니다.
중간 탐지 후보까지 Drive에 저장하므로 재실행하면 가중치 로딩과GPU 추론도 건너뛸 수 있습니다.
'''),
        code('RESOLUTION_SUMMARY = run_resolution_experiment()'),
    ]
    summary = _source(cells[43])
    _set(cells[43], "if not RESOLUTION_ONLY:\n" + textwrap.indent(summary, '    ') + '''
print('해상도별 최종 요약 — ALL, 탄소 RMSE 단위:tC/ha')
resolution_summary_display = RESOLUTION_SUMMARY.loc[RESOLUTION_SUMMARY.siteID=='ALL',
    ['resolution_cm','setting','F1_visible','estimate_reference_ratio','carbon_RMSE','selected_score_threshold']].copy()
resolution_summary_display['setting'] = resolution_summary_display.setting.map({'fixed':'고정','retuned':'재선택(OOF)'})
resolution_summary_display = resolution_summary_display.rename(columns={
    'resolution_cm':'해상도(cm)','setting':'설정','F1_visible':'F1(보이는 나무)',
    'estimate_reference_ratio':'탄소 추정/참조','carbon_RMSE':'RMSE(tC/ha)','selected_score_threshold':'선택 score'})
display(resolution_summary_display.round(3))
print('fixed=.2/600; retuned=patch600, inner F1 선택 → outer OOF 평가')
print('10cm 자동 비교:', OUT/'resolution_10cm_validation.csv')
print('보고서:', OUT/'report.md')
''')
    _set(cells[42], _source(cells[42]) + '\n\n해상도별 F1_visible·탄소 추정/참조·RMSE를 마지막에 추가 출력합니다.')
    cells[39:39] = extra
    return cells
