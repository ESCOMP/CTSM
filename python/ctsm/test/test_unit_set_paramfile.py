#!/usr/bin/env python3

"""Unit tests for set_paramfile"""

import unittest
import os
import sys
import tempfile
import shutil
import numpy as np
import xarray as xr

from ctsm import unit_testing

from ctsm.netcdf_utils import get_netcdf_format
from ctsm.param_utils import set_paramfile as sp
from ctsm.param_utils.paramfile_shared import open_paramfile, are_paramfile_dataarrays_identical

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


def save_paramfile_with_integer_that_has_fillvalue(tempdir):
    """
    Convenience function for creating a parameter file that has an integer parameter with a fill
    value
    """
    input_path = os.path.join(tempdir, "input.nc")
    new_param_name = "new_param_abc123"
    fill_value = -999

    # Construct the Dataset
    data = np.array([1, 2, fill_value, 4, 5], dtype=np.int32)
    da = xr.DataArray(data)
    da.encoding["_FillValue"] = fill_value
    ds = xr.Dataset(data_vars={new_param_name: da})

    # Check it
    assert new_param_name in ds
    assert sp.is_integer(ds[new_param_name].values)

    # Save it
    sp.save_paramfile(ds, input_path)

    return input_path, new_param_name, fill_value


class TestUnitCheckCorrectNdims(unittest.TestCase):
    """Unit tests of check_correct_ndims"""

    def test_checkcorrectndims_0d_int(self):
        """Check True when given a standard int for a 0d parameter"""
        da = xr.DataArray(data=1)
        self.assertTrue(sp.check_correct_ndims(da, 1))

    def test_checkcorrectndims_0d_int_np(self):
        """Check True when given a numpy int for a 0d parameter"""
        da = xr.DataArray(data=1)
        self.assertTrue(sp.check_correct_ndims(da, np.int32(1)))

    def test_checkcorrectndims_1d_int(self):
        """Check False when given a standard int for a 0d parameter"""
        da = xr.DataArray(data=[1, 2])
        self.assertFalse(sp.check_correct_ndims(da, 1))

    def test_checkcorrectndims_1d_int_np(self):
        """Check False when given a numpy int for a 0d parameter"""
        da = xr.DataArray(data=[1, 2])
        self.assertFalse(sp.check_correct_ndims(da, np.int32(1)))

    def test_checkcorrectndims_0d_list(self):
        """Check False when given a list for a 0d parameter"""
        da = xr.DataArray(data=1)
        self.assertFalse(sp.check_correct_ndims(da, [1, 2]))

    def test_checkcorrectndims_0d_nparray(self):
        """Check False when given a numpy array for a 0d parameter"""
        da = xr.DataArray(data=1)
        self.assertFalse(sp.check_correct_ndims(da, np.array([1, 2])))

    def test_checkcorrectndims_0d_nparray_error(self):
        """Check for error when given a numpy array for a 0d parameter and requesting throw_error"""
        da = xr.DataArray(data=1)
        with self.assertRaisesRegex(RuntimeError, "Incorrect N dims: Expected 0, got 1"):
            sp.check_correct_ndims(da, np.array([1, 2]), throw_error=True)

    def test_checkcorrectndims_1d_list(self):
        """Check True when given a list for a 1d parameter"""
        da = xr.DataArray(data=[1, 2])
        self.assertTrue(sp.check_correct_ndims(da, [1, 2]))

    def test_checkcorrectndims_1d_nparray(self):
        """Check True when given a numpy array for a 1d parameter"""
        da = xr.DataArray(data=[1, 2])
        self.assertTrue(sp.check_correct_ndims(da, np.array([1, 2])))


class TestUnitIsInteger(unittest.TestCase):
    """Unit tests of is_integer"""

    def test_isinteger_obj_int(self):
        """Check True if given an object of type int"""
        self.assertTrue(sp.is_integer(int(1)))

    def test_isinteger_obj_int_np(self):
        """Check True if given an object of a numpy integer type"""
        self.assertTrue(sp.is_integer(np.int32(1)))

    def test_isinteger_obj_int_array0d_np(self):
        """Check True if given an object of a numpy array with integer dtype"""
        self.assertTrue(sp.is_integer(np.array(1, dtype=np.int32)))

    def test_isinteger_obj_int_0d_xr(self):
        """Check True if given a numpy scalar integer object via xarray"""
        da = xr.DataArray(data=1)
        self.assertTrue(sp.is_integer(da.values))

    def test_isinteger_obj_int_1d_xr(self):
        """Check True if given a numpy 1-d integer object via xarray"""
        da = xr.DataArray(data=[1])
        self.assertTrue(sp.is_integer(da.values))
        da = xr.DataArray(data=[1, 2])
        self.assertTrue(sp.is_integer(da.values))

    def test_isinteger_obj_float(self):
        """Check False if given an object of type float"""
        self.assertFalse(sp.is_integer(float(3.14)))

    def test_isinteger_obj_float_np(self):
        """Check False if given an object of a numpy float type"""
        self.assertFalse(sp.is_integer(np.float32(3.14)))

    def test_isinteger_type_int(self):
        """Check False if given a type"""
        self.assertFalse(sp.is_integer(int))


class TestUnitSetParamfile(unittest.TestCase):
    """Unit tests of set_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self.orig_argv

    def test_set_paramfile_args_short(self):
        """Test that all arguments can be set correctly with shortnames"""
        output_path = "/path/to/output.nc"
        sys.argv = [
            "get_arguments",
            "-i",
            PARAMFILE,
            "-p",
            "pft1,pft2",
            "-o",
            output_path,
            "param1=new_value1",
            "param2=new_value2",
        ]
        args = sp.get_arguments()
        self.assertEqual(PARAMFILE, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(output_path, args.output)
        self.assertEqual(["param1=new_value1", "param2=new_value2"], args.param_changes)

    def test_set_paramfile_args_long(self):
        """Test that all arguments can be set correctly with longnames"""
        output_path = "/path/to/output.nc"
        sys.argv = [
            "get_arguments",
            "--input",
            PARAMFILE,
            "--pft",
            "pft1,pft2",
            "--drop-other-pfts",
            "--output",
            output_path,
            "param1=new_value1",
            "param2=new_value2",
        ]
        args = sp.get_arguments()
        self.assertEqual(PARAMFILE, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(output_path, args.output)
        self.assertEqual(["param1=new_value1", "param2=new_value2"], args.param_changes)

    def test_set_paramfile_error_missing_input(self):
        """Test that it errors if input file doesn't exist"""
        output_path = "/path/to/output.nc"
        sys.argv = [
            "get_arguments",
            "--input",
            "nwuefweirbdfdiurbe",
            "--pft",
            "pft1,pft2",
            "--output",
            output_path,
        ]
        with self.assertRaises(FileNotFoundError):
            sp.get_arguments()

    def test_set_paramfile_error_existing_output(self):
        """Test that it errors if output file already exists"""
        sys.argv = [
            "get_arguments",
            "--input",
            PARAMFILE,
            "--pft",
            "pft1,pft2",
            "--output",
            PARAMFILE,
        ]
        with self.assertRaises(FileExistsError):
            sp.get_arguments()

    def test_set_paramfile_error_dropotherpfts_without_pft(self):
        """Test that it errors if given --drop-other-pfts without --pft"""
        output_path = "/path/to/output.nc"
        sys.argv = [
            "get_arguments",
            "--input",
            PARAMFILE,
            "--drop-other-pfts",
            "--output",
            output_path,
        ]
        with self.assertRaises(RuntimeError):
            sp.get_arguments()


class TestUnitSaveParamfile(unittest.TestCase):
    """Unit tests of save_paramfile"""

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.output_path = os.path.join(self.tempdir, "output.nc")

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_save_paramfile(self):
        """Test that save_paramfile can save our usual test file to a new file without changes"""
        input_path = PARAMFILE
        ds_in = open_paramfile(input_path)
        sp.save_paramfile(ds_in, self.output_path, nc_format=get_netcdf_format(input_path))
        ds_out = open_paramfile(self.output_path)
        self.assertTrue(ds_out.equals(ds_in))
        self.assertTrue(set(ds_in.variables) == set(ds_out.variables))
        for var in ds_in:
            self.assertTrue(are_paramfile_dataarrays_identical(ds_in[var], ds_out[var]))

    def test_save_paramfile_integer_with_fillvalue(self):
        """Test that save_paramfile can successfully save an integer parameter with a fill value"""

        # Create paramfile with a integer variable with fill value -999
        input_path, new_param_name, fill_value = save_paramfile_with_integer_that_has_fillvalue(
            self.tempdir
        )

        # Read it, checking that its fill value is what we asked for. Note: We need to mask, because
        # otherwise the fill value won't be read.
        ds = xr.open_dataset(input_path, mask_and_scale=True)
        self.assertTrue("_FillValue" in ds[new_param_name].encoding)
        self.assertEqual(fill_value, ds[new_param_name].encoding["_FillValue"])

        # Check that the saved variable is an integer type. Note: We need to NOT mask, because
        # masking converts fill values to NaN, which forces conversion to float.
        ds = xr.open_dataset(input_path, mask_and_scale=False)
        self.assertTrue(sp.is_integer(ds[new_param_name].values))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
