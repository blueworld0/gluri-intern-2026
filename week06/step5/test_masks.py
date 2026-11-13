"""Offline regression check: valid zero-valued mask flags must remain zero."""
import unittest
import numpy as np
from rasterio.transform import from_origin
from rasterio.enums import Resampling
from process_candidate_a import warp


class MaskTests(unittest.TestCase):
    def test_integer_mask_zero_preserved(self):
        profile={"crs":"EPSG:32652","transform":from_origin(400000,3900000,10,10),"height":10,"width":10}
        a=np.zeros((10,10),dtype="uint8");a[4,5]=1
        result=warp(a,profile,profile,Resampling.nearest)
        np.testing.assert_array_equal(result,a)

    def test_invalid_support_stays_zero_under_bilinear(self):
        src={"crs":"EPSG:32652","transform":from_origin(400000,3900000,20,20),"height":5,"width":5}
        dst=dict(src,transform=from_origin(400000,3900000,10,10),height=10,width=10)
        a=np.ones((5,5),dtype="float32");a[2,2]=0
        result=warp(a,src,dst,Resampling.bilinear)
        self.assertLess(result[4,4],1)
        self.assertGreaterEqual(float(result.min()),0)
        self.assertLessEqual(float(result.max()),1)


if __name__=="__main__":unittest.main(verbosity=2)
