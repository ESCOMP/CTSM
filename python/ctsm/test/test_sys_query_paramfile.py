#!/usr/bin/env python3

"""System tests for query_paramfile"""

import unittest
import os
import sys
import io
from contextlib import redirect_stdout
import tempfile
import shutil
import xarray as xr

from ctsm import unit_testing

from ctsm.param_utils import query_paramfile as qp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


class TestSysQueryParamfile(unittest.TestCase):
    """System tests of query_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        sys.argv = self.orig_argv
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_query_paramfile_scalar_nopfts(self):
        """Test that print_values works with scalar parameter and no PFTs specified"""

        sys.argv = ["get_arguments", "-i", PARAMFILE, "phenology_soil_depth"]

        f = io.StringIO()
        with redirect_stdout(f):
            qp.main()
        out = f.getvalue()
        self.assertEqual("phenology_soil_depth: 0.08\n", out)

    def test_query_paramfile_scalar_ignorepfts(self):
        """Test that print_values works with scalar parameter and PFTs specified (ignored)"""

        sys.argv = ["get_arguments", "-i", PARAMFILE, "phenology_soil_depth", "--pft", "c3_crop"]

        f = io.StringIO()
        with redirect_stdout(f):
            qp.main()
        out = f.getvalue()
        self.assertEqual("phenology_soil_depth: 0.08\n", out)

    def test_query_paramfile_pft_selectnone(self):
        """Test that print_values works with PFT-dim parameter and no PFTs specified"""

        sys.argv = ["get_arguments", "-i", PARAMFILE, "rswf_min"]

        f = io.StringIO()
        with redirect_stdout(f):
            qp.main()
        out = f.getvalue()
        self.assertRegex(
            out,
            (
                r"rswf_min:\n"
                r"\s+not_vegetated\s*: 0\.25\n"
                r"\s+needleleaf_evergreen_temperate_tree\s*: 0\.25\n.*"
            ),
        )

    def test_query_paramfile_pft_select2(self):
        """Test that print_values works with PFT-dim parameter and two PFTs specified"""

        sys.argv = [
            "get_arguments",
            "-i",
            PARAMFILE,
            "--pft",
            "not_vegetated,needleleaf_evergreen_temperate_tree",
            "rswf_min",
        ]

        f = io.StringIO()
        with redirect_stdout(f):
            qp.main()
        out = f.getvalue()
        self.assertRegex(
            out,
            (
                r"rswf_min:\n"
                r"\s+not_vegetated\s*: 0\.25\n"
                r"\s+needleleaf_evergreen_temperate_tree\s*: 0\.25\n"
            ),
        )

    def test_query_paramfile_no_variables_fake(self):
        """
        Test that print_values prints every variable when no variables are given. Use a small fake
        paramfile so we can check that what gets printed is what we expect.
        """

        fake_da1 = xr.DataArray(data=[1, 2, 3])
        fake_da2 = xr.DataArray(data=[4, 5, 6])
        fake_ds = xr.Dataset(data_vars={"fake1": fake_da1, "fake2": fake_da2})
        fake_nc_path = os.path.join(self.tempdir, "fake_da.nc")
        fake_ds.to_netcdf(fake_nc_path)

        sys.argv = [
            "query_paramfile",
            "-i",
            fake_nc_path,
        ]

        f = io.StringIO()
        with redirect_stdout(f):
            qp.main()
        out = f.getvalue()
        self.assertRegex(out, (r"fake1: \[1 2 3\]\nfake2: \[4 5 6\]\n"))

    def test_query_paramfile_no_variables_real(self):
        """
        Test that query_paramfile doesn't error when trying to print every variable from a real
        paramfile. Don't actually check that it matches what we expect; that's done in
        test_query_paramfile_no_variables_fake.
        """

        sys.argv = [
            "query_paramfile",
            "-i",
            PARAMFILE,
        ]

        f = io.StringIO()
        with redirect_stdout(f):
            qp.main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
