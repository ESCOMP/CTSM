#!/usr/bin/env python3
"""
Unit tests for netcdf_utils.py functions
"""

import os
import sys
import unittest
import shutil
import tempfile
import numpy as np
import pandas as pd
import xarray as xr

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
import ctsm.netcdf_utils as nu
from ctsm import unit_testing

# pylint: disable=invalid-name


class TestUnitGetNetcdfFormat(unittest.TestCase):
    """
    Unit tests for get_netcdf_format
    """

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.outfile = os.path.join(self.tempdir, "file.nc")
        da = xr.DataArray(data=[1, 2, 3])
        self.ds = xr.Dataset(data_vars={"var": da})

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_get_netcdf_format_classic(self):
        """
        Test that get_netcdf_format() gets "classic" format right
        """
        nc_format = "NETCDF3_CLASSIC"
        self.ds.to_netcdf(self.outfile, format=nc_format)
        self.assertEqual(nu.get_netcdf_format(self.outfile), nc_format)

    def test_get_netcdf_format_netcdf4(self):
        """
        Test that get_netcdf_format() gets "netCDF4" format right
        """
        nc_format = "NETCDF4"
        self.ds.to_netcdf(self.outfile, format=nc_format)
        self.assertEqual(nu.get_netcdf_format(self.outfile), nc_format)


class TestUnitAreXrDataArraysIdentical(unittest.TestCase):
    """
    Unit tests for are_xr_dataarrays_identical
    """

    # pylint: disable=too-many-public-methods

    def test_are_xr_dataarrays_identical_data_types_match(self):
        """Should be true if dtypes match"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_data_types_differ(self):
        """Should be false if dtypes differ"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(np.float64(1))
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_data_types_differ_precision(self):
        """Should be false if dtypes differ only in their precision"""
        da0 = xr.DataArray(np.float64(1))
        da1 = xr.DataArray(np.float32(1))
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_encodings_match(self):
        """Should be true if encodings are the exact same"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.encoding["dummy0"] = -999
        da1.encoding["dummy0"] = -999
        da0.encoding["dummy1"] = -999
        da1.encoding["dummy1"] = -999
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_encodings_match_except_order(self):
        """Should be true if encodings are the same as long as you don't care about order"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.encoding["dummy0"] = -999
        da0.encoding["dummy1"] = -999
        da1.encoding["dummy1"] = -999
        da1.encoding["dummy0"] = -999
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_encodings_differ_number(self):
        """Should be false if encodings have a different number of keys"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.encoding["dummy0"] = -999
        da1.encoding["dummy0"] = -999
        da1.encoding["dummy1"] = -999
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_encodings_differ(self):
        """Should be false if encodings have the same number and order of keys but not values"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.encoding["dummy0"] = -999
        da1.encoding["dummy0"] = -999
        da0.encoding["dummy1"] = -998
        da1.encoding["dummy1"] = -999
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_attributes_match(self):
        """Should be true if attributes are the exact same"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.attrs["dummy0"] = -999
        da1.attrs["dummy0"] = -999
        da0.attrs["dummy1"] = -999
        da1.attrs["dummy1"] = -999
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_attributes_match_nan(self):
        """Should be true if attributes are the exact same and NaN"""
        da0 = xr.DataArray(1.0)
        da1 = xr.DataArray(1.0)
        da0.attrs["_FillValue"] = np.float64(np.nan)
        da1.attrs["_FillValue"] = np.float64(np.nan)
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_attributes_match_except_order(self):
        """Should be true if attributes are the same as long as you don't care about order"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.attrs["dummy0"] = -999
        da0.attrs["dummy1"] = -999
        da1.attrs["dummy1"] = -999
        da1.attrs["dummy0"] = -999
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_attributes_differ_number(self):
        """Should be false if attributes have a different number of keys"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.attrs["dummy0"] = -999
        da1.attrs["dummy0"] = -999
        da1.attrs["dummy1"] = -999
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_attributes_differ(self):
        """Should be false if attributes have the same number and order of keys but not values"""
        da0 = xr.DataArray(int(1))
        da1 = xr.DataArray(int(1))
        da0.attrs["dummy0"] = -999
        da1.attrs["dummy0"] = -999
        da0.attrs["dummy1"] = -998
        da1.attrs["dummy1"] = -999
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_values_ndarrays_match(self):
        """Should be true if values match and they're both numpy arrays under the hood"""
        da0 = xr.DataArray(np.array(int(1)))
        da1 = xr.DataArray(np.array(int(1)))
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_values_ndarrays_differ(self):
        """Should be false if values differ and they're both numpy arrays under the hood"""
        da0 = xr.DataArray(np.array(int(1)))
        da1 = xr.DataArray(np.array(int(2)))
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_coords_match(self):
        """Should be true if coordinates match"""
        time = pd.date_range("2000-01-01", periods=3)
        da0 = xr.DataArray(
            dims=["time"],
            coords={"time": time},
        )
        da1 = xr.DataArray(
            dims=["time"],
            coords={"time": time},
        )
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

        # To make sure we don't need a separate "test_are_xr_dataarrays_identical_indexes_match"
        self.assertEqual(len(da0.indexes), len(da1.indexes))
        for key in da0.indexes:
            self.assertTrue(da0.indexes[key].equals(da1.indexes.get(key)))

    def test_are_xr_dataarrays_identical_coords_onemissing(self):
        """Should be false if only one has coords"""
        data = [1, 2, 3]
        time = pd.date_range("2000-01-01", periods=len(data))
        da0 = xr.DataArray(
            data=data,
            dims=["time"],
            coords={"time": time},
        )
        da1 = xr.DataArray(
            data=data,
            dims=["time"],
        )
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_coords_differ(self):
        """Should be false if coordinates differ"""
        time0 = pd.date_range("2000-01-01", periods=3)
        da0 = xr.DataArray(
            dims=["time"],
            coords={"time": time0},
        )
        time1 = pd.date_range("1987-01-01", periods=3)
        da1 = xr.DataArray(
            dims=["time"],
            coords={"time": time1},
        )
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_dims_match(self):
        """Should be true if dimensions match"""
        da0 = xr.DataArray(data=[1], dims=["dim"])
        da1 = xr.DataArray(data=[1], dims=["dim"])
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_dims_differ(self):
        """Should be false if dimensions differ"""
        da0 = xr.DataArray(data=[1], dims=["dim0"])
        da1 = xr.DataArray(data=[1], dims=["dim1"])
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_names_match(self):
        """Should be true if names match"""
        da0 = xr.DataArray(data=[1], name="da_name")
        da1 = xr.DataArray(data=[1], name="da_name")
        self.assertTrue(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_names_differ(self):
        """Should be false if names differ"""
        da0 = xr.DataArray(data=[1], name="da0_name")
        da1 = xr.DataArray(data=[1], name="da1_name")
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    def test_are_xr_dataarrays_identical_sizes_differ(self):
        """Should be false if sizes differ"""
        da0 = xr.DataArray(data=[1])
        da1 = xr.DataArray(data=[1, 1])
        self.assertFalse(nu.are_xr_dataarrays_identical(da0, da1))

    # Waiting on dask, sparse, or pint to be in ctsm_pylib:
    # TODO: False if data types don't match
    # TODO: NotImplementedError if data types match but aren't np.ndarray

    # Waiting on dask to be in ctsm_pylib:
    # TODO: True if the only difference is chunked or not
    # TODO: True if the only difference is chunk sizes


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
