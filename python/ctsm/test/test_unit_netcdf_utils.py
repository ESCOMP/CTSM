#!/usr/bin/env python3
"""
Unit tests for netcdf_utils.py functions
"""

import os
import sys
import unittest
import shutil
import tempfile
import xarray as xr

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm.netcdf_utils import get_netcdf_format
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
        self.assertEqual(get_netcdf_format(self.outfile), nc_format)

    def test_get_netcdf_format_netcdf4(self):
        """
        Test that get_netcdf_format() gets "netCDF4" format right
        """
        nc_format = "NETCDF4"
        self.ds.to_netcdf(self.outfile, format=nc_format)
        self.assertEqual(get_netcdf_format(self.outfile), nc_format)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
