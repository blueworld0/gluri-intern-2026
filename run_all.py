"""모든 단계를 순서대로 실행:  python run_all.py"""
import runpy, sys
from pathlib import Path
src = Path(__file__).parent / "src"
sys.path.insert(0, str(src))
for step in ["00_download.py", "01_boundary.py", "01b_species_r.py", "02_carbon.py", "03_maps.py"]:
    print(f"\n===== {step} =====")
    runpy.run_path(str(src / step), run_name="__main__")
