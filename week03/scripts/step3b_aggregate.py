"""Step 3b: does averaging GEDI shots into grid cells (less label noise) raise R²?"""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step3_diagnose import PROC, metrics, spatial_oof  # noqa: E402

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')


def main():
    names = json.loads((PROC / 'feature_names.json').read_text(encoding='utf-8'))
    df = pd.read_csv(PROC / 'training_all_samples.csv', dtype={'shot_number': 'uint64'})
    df = df[df.valid_target & df.agbd.le(600) & df.complete_features].reset_index(drop=True)
    rows = [{'cell_m': 'shot (25m)', 'min_shots': 1, **metrics(df.agbd, spatial_oof(df, names))}]
    for cell in (100, 250, 500):
        df['cell'] = (np.floor(df.easting / cell).astype(int).astype(str) + '_'
                      + np.floor(df.northing / cell).astype(int).astype(str))
        g = df.groupby('cell').agg(n=('agbd', 'size'), agbd=('agbd', 'mean'),
                                   easting=('easting', 'mean'), northing=('northing', 'mean'),
                                   **{f: (f, 'mean') for f in names}).reset_index()
        for min_shots in (3, 5):
            sub = g[g.n >= min_shots].copy()
            # Same 2 km spatial blocks as before, so results stay comparable.
            sub['block_id'] = (np.floor(sub.easting / 2000).astype(int).astype(str) + '_'
                               + np.floor(sub.northing / 2000).astype(int).astype(str))
            if len(sub) < 50 or sub.block_id.nunique() < 5:
                continue
            rows.append({'cell_m': cell, 'min_shots': min_shots,
                         **metrics(sub.agbd, spatial_oof(sub, names))})
            print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(PROC / 'step3_aggregation.csv', index=False, encoding='utf-8-sig')
    print(out.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
