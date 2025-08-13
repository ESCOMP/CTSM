#!/usr/bin/env python3

"""Unit tests for set_paramfile"""

import unittest
import os
import sys
import numpy as np
import xarray as xr

from ctsm import unit_testing

from ctsm.param_utils import set_paramfile as sp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


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
            "-v",
            "var1,var2",
            "param1=new_value1",
            "param2=new_value2",
        ]
        args = sp.get_arguments()
        self.assertEqual(PARAMFILE, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(["var1", "var2"], args.variables)
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
            "--output",
            output_path,
            "--variables",
            "var1,var2",
            "param1=new_value1",
            "param2=new_value2",
        ]
        args = sp.get_arguments()
        self.assertEqual(PARAMFILE, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(["var1", "var2"], args.variables)
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


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
