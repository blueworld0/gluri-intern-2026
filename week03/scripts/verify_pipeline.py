"""Synthetic verification only. Never writes research output into data/ or figures/."""
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import h5py
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
import run_analysis as pipeline


class PipelineTest(unittest.TestCase):
    def test_uint64_filters_raster_and_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            h5 = root / 'GEDI04_A_synthetic.h5'
            with h5py.File(h5, 'w') as f:
                for i, beam in enumerate(pipeline.BEAMS):
                    g = f.create_group(beam)
                    values = {
                        'shot_number': np.array([2**63 + i*10 + j for j in range(6)], dtype='uint64'),
                        'lat_lowestmode': [37.7]*5+[38.0],
                        'lon_lowestmode': [128.1]*6,
                        'agbd': [100., -9999., 200., 300., 400., 500.],
                        'agbd_se': [10.]*6,
                        'l4_quality_flag': [1, 1, 0, 1, 1, 1],
                        'l2_quality_flag': [1]*6,
                        'degrade_flag': [0, 0, 0, 1, 0, 0],
                        'sensitivity': [.98]*6,
                        'solar_elevation': [-10., 1., 1., 1., 10., 1.],
                        'delta_time': [94694400.]*6,
                    }
                    for name, value in values.items():
                        g.create_dataset(name, data=value)
                    g['agbd'].attrs['_FillValue'] = -9999.
            all_points, log = pipeline.read_bbox([h5])
            self.assertEqual(len(log), 8)
            self.assertEqual(len(all_points), 40)
            self.assertEqual(str(all_points.shot_number.dtype), 'uint64')
            filtered, summary = pipeline.quality_filter(all_points)
            self.assertEqual(summary['count'].tolist(), [40, 32, 24, 16, 8])
            tile = root / pipeline.WC_NAME
            with rasterio.open(tile, 'w', driver='GTiff', height=4, width=4,
                               count=1, dtype='uint8', crs='EPSG:4326',
                               transform=from_origin(127.98, 37.85, .1, .1), nodata=0) as ds:
                ds.write(np.full((4, 4), 10, dtype='uint8'), 1)
            filtered = pipeline.worldcover_values(filtered, root)
            self.assertEqual(int(filtered.tree_cover.sum()), 16)
            with patch.object(pipeline, 'ROOT', root):
                pipeline.save_results(all_points, filtered, summary, root, root, 1)
            saved = pd.read_csv(root / 'gedi_l4a_2021_filtered.csv', dtype={'shot_number': 'uint64'})
            self.assertEqual(saved.shot_number.tolist(), filtered.shot_number.tolist())
            self.assertTrue((root / 'map_points.png').is_file())
            self.assertTrue((root / 'agbd_hist.png').is_file())
            self.assertIn('점 부족', (root / 'summary.md').read_text(encoding='utf-8'))
            # Zero retained points must still produce valid CSVs, plots, and report.
            empty = filtered.iloc[:0].copy()
            empty_filtered, empty_summary = pipeline.quality_filter(all_points.iloc[:0])
            empty_filtered = pipeline.worldcover_values(empty_filtered, root)
            with patch.object(pipeline, 'ROOT', root):
                pipeline.save_results(all_points.iloc[:0], empty_filtered, empty_summary, root, root, 1)
            self.assertEqual(empty_summary['count'].tolist(), [0]*5)


if __name__ == '__main__':
    unittest.main()
