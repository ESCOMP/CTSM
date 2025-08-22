#!/usr/bin/env python3

"""Unit tests for query_paramfile"""

import unittest
import sys
import io
from contextlib import redirect_stdout
import xarray as xr

from ctsm import unit_testing

from ctsm.param_utils import query_paramfile as qp
from ctsm.param_utils.paramfile_shared import PFTNAME_VAR

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


def _setup_pft_parameter_ds():
    """
    Set up a parameter Dataset with a PFT-dimensioned parameter
    """
    pft_dimname = "pft"
    pft_names_list = ["pft0", "pft1"]
    pft_names_da = xr.DataArray(
        name=PFTNAME_VAR,
        data=pft_names_list,
        dims=[pft_dimname],
        coords={pft_dimname: pft_names_list},
    )
    var_name = "pft_param"
    pft_param_da = xr.DataArray(
        data=[1986.0325, 1987.0724], dims=[pft_dimname], coords={PFTNAME_VAR: pft_names_da}
    )
    ds = xr.Dataset(data_vars={var_name: pft_param_da, pft_dimname: pft_names_da})
    return ds, var_name, pft_names_list


class TestUnitQueryParamfile(unittest.TestCase):
    """Unit tests of query_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self.orig_argv

    def test_query_paramfile_args_short(self):
        """Test that all arguments can be set correctly with shortnames"""
        input_path = "/path/to/input.nc"
        sys.argv = ["get_arguments", "-i", input_path, "-p", "pft1,pft2", "var1", "var2"]
        args = qp.get_arguments()
        self.assertEqual(input_path, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(["var1", "var2"], args.variables)

    def test_query_paramfile_args_long(self):
        """Test that all arguments can be set correctly with longnames"""
        input_path = "/path/to/input.nc"
        sys.argv = ["get_arguments", "--input", input_path, "--pft", "pft1,pft2", "var1", "var2"]
        args = qp.get_arguments()
        self.assertEqual(input_path, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(["var1", "var2"], args.variables)

    def test_query_paramfile_print_scalar(self):
        """Test that print_values works with a scalar parameter"""
        scalar_da = xr.DataArray(data=1987.0724)
        var_name = "scalar_param"
        ds = xr.Dataset(data_vars={var_name: scalar_da})

        f = io.StringIO()
        with redirect_stdout(f):
            qp.print_values(ds, var_name, selected_pfts=None, pft_names=None)
        out = f.getvalue()
        self.assertEqual("scalar_param: 1987.0724\n", out)

    def test_query_paramfile_print_pfts_selectnone(self):
        """Test that print_values works with PFT-dimensioned parameter, selecting no PFTs"""
        ds, var_name, pft_names_list = _setup_pft_parameter_ds()

        f = io.StringIO()
        with redirect_stdout(f):
            qp.print_values(ds, var_name, selected_pfts=None, pft_names=pft_names_list)
        out = f.getvalue()
        self.assertEqual("pft_param:\n   pft0: 1986.0325\n   pft1: 1987.0724\n", out)

    def test_query_paramfile_print_pfts_selectall(self):
        """Test that print_values works with PFT-dimensioned parameter, selecting all PFTs"""
        ds, var_name, pft_names_list = _setup_pft_parameter_ds()

        f = io.StringIO()
        with redirect_stdout(f):
            qp.print_values(ds, var_name, selected_pfts=pft_names_list, pft_names=pft_names_list)
        out = f.getvalue()
        self.assertEqual("pft_param:\n   pft0: 1986.0325\n   pft1: 1987.0724\n", out)

    def test_query_paramfile_print_pfts_selectone(self):
        """Test that print_values works with PFT-dimensioned parameter, selecting one PFT"""
        ds, var_name, pft_names_list = _setup_pft_parameter_ds()

        f = io.StringIO()
        with redirect_stdout(f):
            qp.print_values(ds, var_name, selected_pfts=["pft0"], pft_names=pft_names_list)
        out = f.getvalue()
        self.assertEqual("pft_param:\n   pft0: 1986.0325\n", out)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
