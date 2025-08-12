#!/usr/bin/env python3

"""System tests for set_paramfile"""

import unittest
import os
import sys
import shutil
import tempfile
import xarray as xr
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm import unit_testing

from ctsm.netcdf_utils import get_netcdf_format
from ctsm.param_utils import set_paramfile as sp
from ctsm.param_utils.paramfile_shared import open_paramfile

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


class TestSysSetParamfile(unittest.TestCase):
    """System tests of set_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        sys.argv = self.orig_argv
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_set_paramfile_copyfile(self):
        """Test that set_paramfile can straight-up copy to a new file"""
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = ["set_paramfile", "-i", PARAMFILE, "-o", output_path]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that contents are functionally identical
        self.assertEqual(ds_in, ds_out)

        # Check that both are the same kind of netCDF
        self.assertEqual(get_netcdf_format(PARAMFILE), get_netcdf_format(output_path))

    def test_set_paramfile_extractvars(self):
        """Test that set_paramfile can copy to a new file with only some requested variables"""
        output_path = os.path.join(self.tempdir, "output.nc")
        vars_to_include = ["a_coef", "bgc_cn_s2"]
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-v",
            ",".join(vars_to_include),
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that output file variables/coords are what we expect
        # (use set() so order doesn't matter)
        expected_var_list = set(vars_to_include + list(ds_in.coords))
        self.assertEqual(set(ds_out.variables), expected_var_list)

        # Check that included variables/coords match
        for var in vars_to_include:
            self.assertEqual(ds_in[var], ds_out[var])

        # Check that both are the same kind of netCDF
        self.assertEqual(get_netcdf_format(PARAMFILE), get_netcdf_format(output_path))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
