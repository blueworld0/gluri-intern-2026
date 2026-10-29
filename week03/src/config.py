"""설정값을 한 곳에 모아 둔 파일. 값을 바꾸고 싶으면 여기만 고치면 됩니다."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"          # 내려받은 원자료 (깃허브에는 안 올림)
OUT = ROOT / "outputs"        # 계산 결과 (표, 지도)

# 작업 범위 (경도 서, 위도 남, 경도 동, 위도 북). 금정산~백양산을 넉넉히 포함.
BBOX = (128.98, 35.12, 129.14, 35.36)

# 공원 경계 파일. 공식 경계(국립공원공단)를 받으면 그 파일 경로로 바꾸세요.
# 없으면 01_boundary.py 가 만든 '임시 경계'를 씁니다.
OFFICIAL_BOUNDARY = None      # 예: DATA / "boundary" / "geumjeongsan_np_official.geojson"
PROXY_BOUNDARY = DATA / "boundary" / "geumjeongsan_proxy.geojson"
OFFICIAL_AREA_KM2 = 66.859    # 환경부 고시 면적 (2026-03-03 지정)

# 계산에 쓰는 좌표계: UTM 52N (미터 단위라 면적 계산이 쉬움)
CRS_M = "EPSG:32652"

# ESA CCI Biomass v4 (100 m, 단위 Mg/ha = t/ha). NASA MAAP 공개 저장소에서 받을 수 있는 연도.
CCI_YEARS = [2017, 2018, 2019, 2020]
CCI_URL = ("https://nasa-maap-data-store.s3.us-west-2.amazonaws.com/file-staging/nasa-map/"
           "ESACCI_Biomass_L4_AGB_V4_100m_{y}/N40E120_ESACCI-BIOMASS-L4-AGB{sd}-MERGED-100m-{y}-fv4.0.tif")

# ESA WorldCover (10 m 토지피복). 10 = 나무(Tree cover)
# 작업 범위가 동경 129도 선에 걸쳐 있어서 서쪽 타일(126E)과 동쪽 타일(129E)을 둘 다 씁니다.
_WC = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/{v}/{y}/map/ESA_WorldCover_10m_{y}_{v}_N33E{e}_Map.tif"
WORLDCOVER = {
    2020: [_WC.format(v="v100", y=2020, e=e) for e in ("126", "129")],
    2021: [_WC.format(v="v200", y=2021, e=e) for e in ("126", "129")],
}
TREE_CLASS = 10
_DEM = ("https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N35_00_E{e}_00_DEM/"
        "Copernicus_DSM_COG_10_N35_00_E{e}_00_DEM.tif")
DEM_URL = [_DEM.format(e=e) for e in ("128", "129")]

# 탄소 환산 계수 (설악산 논문 전청옥 외 2024 과 같은 IPCC 2006 기본값)
R_ROOT = 0.43     # 뿌리(지하부)/지상부 비율
CF = 0.5          # 바이오매스 중 탄소 비율
CO2_PER_C = 44 / 12

# ── 수종별 뿌리비율 R (week03 개선) ─────────────────────────────────
# "ipcc"    : 위의 R_ROOT 하나(0.43)를 모든 칸에 적용 (설악산 논문 방식)
# "species" : 임상도 수종별로 칸마다 다른 R 적용 (01b_species_r.py 가 만든 R 지도 사용)
R_MODE = "species"

# 임상도 1:5,000 (산림청, 공공데이터포털). 바이오매스 지도(2017–2020)와 시기를 맞추려고 2019년판 사용.
# zip 을 풀지 않고 그대로 읽습니다.  "zip://경로!zip안의.shp"
IMSANG = [
    "zip://" + str(DATA / "imsang" / "busan_2019.zip") + "!26.shp",      # 부산광역시
    "zip://" + str(DATA / "imsang" / "yangsan_2019.zip") + "!48330.shp", # 경남 양산시
]
R_MAP = DATA / "imsang" / "r_map_2019.tif"   # 01b 가 만드는 칸별 R 지도

# 수종 코드(KOFTR_GROU) → R.  값 출처: Lee et al. 2018 (Forests), 국립산림과학원 수종별 배출계수
#   표에 있는 수종은 그 값, 표에 없는 묶음은 아래 규칙으로 정함 (근거는 README 참고)
_R_CON = (0.25 + 0.29 + 0.36 + 0.29) / 4     # 표의 침엽수 4종 평균 (소나무·곰솔·리기다·낙엽송) ≈ 0.30
_R_OAK = (0.31 + 0.39 + 0.32) / 3            # 참나무 3종 평균 (상수리·신갈·굴참) = 0.34
_R_BRD = (0.31 + 0.39 + 0.32 + 0.48) / 4     # 표의 활엽수 4종 평균 (참나무 3종 + 아까시) = 0.375
SPECIES_R = {
    "11": ("소나무", 0.25), "15": ("곰솔", 0.29), "14": ("리기다소나무", 0.36), "13": ("낙엽송", 0.29),
    "31": ("상수리나무", 0.31), "32": ("신갈나무", 0.39), "33": ("굴참나무", 0.32), "49": ("아까시나무", 0.48),
    "34": ("기타참나무류", _R_OAK),
    "10": ("기타침엽수", _R_CON), "12": ("잣나무", _R_CON), "17": ("편백나무", _R_CON), "18": ("삼나무", _R_CON),
    "30": ("기타활엽수", _R_BRD), "35": ("오리나무", _R_BRD), "39": ("밤나무", _R_BRD), "41": ("서어나무", _R_BRD),
    "45": ("포플러", _R_BRD), "46": ("벚나무", _R_BRD), "47": ("느티나무", _R_BRD),
    # 침활혼효림: 금정산의 대표 침엽수(곰솔)와 대표 활엽수(기타참나무류)를 반반
    "77": ("침활혼효림", (0.29 + _R_OAK) / 2),
}
# 죽림·비산림 등 표에 없는 칸은 R 지도에서 0 으로 두고, 계산할 때 '산림 면적가중 평균 R'로 채움

# 민감도 분석용: 한국 수종별 국가계수 (Lee et al. 2018, 국립산림과학원 계열)
R_SENSITIVITY = {"IPCC 기본(설악산 논문)": 0.43, "소나무(중부)": 0.25, "곰솔": 0.29, "신갈나무": 0.39}
