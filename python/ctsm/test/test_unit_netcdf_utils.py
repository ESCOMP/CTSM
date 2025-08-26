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
# pylint: disable=protected-access


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

    def test_are_xr_dataarrays_identical_values_ndarrays_differ_nantypeerror(self):
        """
        Should be false if values differ, they're both numpy arrays under the hood, and they can't
        be coerced to a type capable of NaN
        """
        da0 = xr.DataArray(np.array(["a", "b", "c"]))
        da1 = xr.DataArray(np.array(["d", "e", "f"]))
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
    # TODO: NotImplementedError if data types match but aren't np.array

    # Waiting on dask to be in ctsm_pylib:
    # TODO: True if the only difference is chunked or not
    # TODO: True if the only difference is chunk sizes


class TestUnitAreDictsIdenticalNansEqual(unittest.TestCase):
    """
    Unit tests for _are_dicts_identical_nansequal
    """

    def test_are_dicts_identical_nansequal_yes_noignore_nonan(self):
        """
        Test two identical dicts with no keys being ignored and no nans
        """
        dict0 = {"a": 1, "b": 2}
        dict1 = {"a": 1, "b": 2}
        self.assertTrue(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_yes_ignorestr_nonan(self):
        """
        Test two identical dicts with a key being ignored (as str) and no nans
        """
        dict0 = {"a": 1, "b": 2}
        dict1 = {"a": 1, "b": 3}
        self.assertTrue(nu._are_dicts_identical_nansequal(dict0, dict1, keys_to_ignore="b"))

    def test_are_dicts_identical_nansequal_yes_ignorelist_nonan(self):
        """
        Test two identical dicts with a key being ignored (as list) and no nans
        """
        dict0 = {"a": 1, "b": 2}
        dict1 = {"a": 1, "b": 3}
        self.assertTrue(nu._are_dicts_identical_nansequal(dict0, dict1, keys_to_ignore=["b"]))

    def test_are_dicts_identical_nansequal_yes_ignorestr_nan(self):
        """
        Test two identical dicts with a key being ignored (as str) and one key with matching nan
        values
        """
        dict0 = {"a": np.nan, "b": 2}
        dict1 = {"a": np.nan, "b": 3}
        self.assertTrue(nu._are_dicts_identical_nansequal(dict0, dict1, keys_to_ignore="b"))

    def test_are_dicts_identical_nansequal_no_nonan(self):
        """
        Test two different dicts with no nans
        """
        dict0 = {"a": 1, "b": 2}
        dict1 = {"a": 1, "b": 4}
        self.assertFalse(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_no_nansmatch(self):
        """
        Test two different dicts with matching nans
        """
        dict0 = {"a": np.nan, "b": 2}
        dict1 = {"a": np.nan, "b": 4}
        self.assertFalse(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_no_nansdiffer(self):
        """
        Test two different dicts that only differ in that one has a nan where another doesn't
        """
        dict0 = {"a": 1, "b": 2}
        dict1 = {"a": 1, "b": np.nan}
        self.assertFalse(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_no_keysdiffer(self):
        """
        Test two different dicts that have different keys
        """
        dict0 = {"a": 1, "b": 2}
        dict1 = {"a": 1, "c": 2}
        self.assertFalse(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_no_lengthsdiffer(self):
        """
        Test two different dicts that have different lengths
        """
        dict0 = {"a": 1, "b": 2}
        dict1 = {"a": 1}
        self.assertFalse(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_no_both_nparrays(self):
        """
        Test two different dicts that have differing numpy arrays for one value
        """
        dict0 = {"a": 1, "b": np.array([1, 2])}
        dict1 = {"a": 1, "b": np.array([1, 3])}
        self.assertFalse(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_no_differ_nparrays(self):
        """
        Test two dicts where one has a value that's a numpy array and the other doesn't, but they're
        identical if you coerce them both to numpy arrays
        """
        dict0 = {"a": 1, "b": np.array([1, 2])}
        dict1 = {"a": 1, "b": [1, 2]}
        self.assertTrue(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_yes_both_nparrays(self):
        """
        Test two different dicts that have identical numpy arrays for one value
        """
        dict0 = {"a": 1, "b": np.array([1, 2])}
        dict1 = {"a": 1, "b": np.array([1, 2])}
        self.assertTrue(nu._are_dicts_identical_nansequal(dict0, dict1))

    def test_are_dicts_identical_nansequal_yes_both_nparrays_str(self):
        """
        Test two different dicts that have identical numpy arrays of strings for one value
        """
        dict0 = {"a": 1, "b": np.array(["1", "2"])}
        dict1 = {"a": 1, "b": np.array(["1", "2"])}
        self.assertTrue(nu._are_dicts_identical_nansequal(dict0, dict1))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
