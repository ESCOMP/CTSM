#!/usr/bin/env python3

"""Unit tests for grid_one_variable.py"""

import tempfile
import shutil
import unittest
import os
import numpy as np
import xarray as xr

from ctsm import unit_testing
from ctsm.crop_calendars import grid_one_variable as g1v

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestGridOneVariable(unittest.TestCase):
    """Tests of grid_one_variable"""

    def setUp(self):
        self.lat_da = xr.DataArray(
            data=[-45.0, 45.0],
            dims=["lat"],
        )
        self.lon_da = xr.DataArray(
            data=[-145.0, 145.0],
            dims=["lon"],
        )
        self.result = None
        self.expected = None

    def _compare_arrays(self):
        print(f"Result:\n{self.result}")
        print(f"Expected:\n{self.expected}")
        self.assertTrue(np.array_equal(self.result, self.expected, equal_nan=True))

    def test_create_filled_array_fillnan(self):
        """Test create_filled_array() with NaN fill_value"""
        this_da = xr.DataArray(
            data=np.array([[1, 2], [3, 4]]),
            dims=["lat", "lon"],
            coords={"lat": self.lat_da, "lon": self.lon_da},
        )
        dummy_ds = xr.Dataset()

        fill_value = np.nan
        self.result = g1v.create_filled_array(dummy_ds, fill_value, this_da, ["lat"])

        self.expected = np.full_like(self.lat_da, fill_value)
        self._compare_arrays()

    def test_create_filled_array_fill6(self):
        """Test create_filled_array() with fill_value = 6.0"""
        this_da = xr.DataArray(
            data=np.array([[1, 2], [3, 4]]),
            dims=["lat", "lon"],
            coords={"lat": self.lat_da, "lon": self.lon_da},
        )
        dummy_ds = xr.Dataset()

        fill_value = 6.0
        self.result = g1v.create_filled_array(dummy_ds, fill_value, this_da, ["lat"])

        self.expected = np.full_like(self.lat_da, fill_value)
        self._compare_arrays()

    def test_create_filled_array_fillfalse(self):
        """Test create_filled_array() with false(y) fill_value"""
        this_da = xr.DataArray(
            data=np.array([[1, 2], [3, 4]]),
            dims=["lat", "lon"],
            coords={"lat": self.lat_da, "lon": self.lon_da},
        )
        dummy_ds = xr.Dataset()

        fill_value = False
        self.result = g1v.create_filled_array(dummy_ds, fill_value, this_da, ["lat"])

        self.expected = np.full_like(self.lat_da, np.nan)
        self._compare_arrays()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
