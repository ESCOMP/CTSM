#!/usr/bin/env python3
"""
Unit tests for subset_data

You can run this by:
    python -m unittest test_unit_subset_data.py
"""

import unittest
import configparser
import argparse
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.subset_data import get_parser, setup_files, check_args
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=invalid-name


class TestSubsetData(unittest.TestCase):
    """
    Basic class for testing SubsetData class in subset_data.py.
    """

    def setUp(self):
        sys.argv = ["subset_data", "point", "--create-surface"]
        DEFAULTS_FILE = os.path.join(
            os.getcwd(), "../tools/site_and_regional/default_data_2000.cfg"
        )
        self.parser = get_parser()
        self.args = self.parser.parse_args()
        self.cesmroot = path_to_ctsm_root()
        self.defaults = configparser.ConfigParser()
        self.defaults.read(os.path.join(self.cesmroot, "tools/site_and_regional", DEFAULTS_FILE))

    def test_inputdata_setup_files_basic(self):
        """
        Test
        """
        check_args(self.args)
        files = setup_files(self.args, self.defaults, self.cesmroot)
        self.assertEqual(
            files["fsurf_in"],
            "surfdata_0.9x1.25_hist_2000_16pfts_c240216.nc",
            "fsurf_in filename not whats expected",
        )
        self.assertEqual(
            files["fsurf_out"],
            None,
            "fsurf_out filename not whats expected",
        )
        self.assertEqual(
            files["main_dir"],
            "/glade/campaign/cesm/cesmdata/cseg/inputdata",
            "main_dir directory not whats expected",
        )

    def test_inputdata_setup_files_inputdata_dne(self):
        """
        Test that inputdata directory does not exist
        """
        check_args(self.args)
        self.defaults.set("main", "clmforcingindir", "/zztop")
        with self.assertRaisesRegex(SystemExit, "inputdata directory does not exist"):
            setup_files(self.args, self.defaults, self.cesmroot)

    def test_check_args_nooutput(self):
        """
        Test that check args aborts when no-output is asked for
        """
        sys.argv = ["subset_data", "point"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(argparse.ArgumentError, "Must supply one of"):
            check_args(self.args)

    def test_check_args_notype(self):
        """
        Test that check args aborts when no type is asked for
        """
        sys.argv = ["subset_data"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(argparse.ArgumentError, "Must supply a positional argument:"):
            check_args(self.args)

    def test_check_args_badconfig(self):
        """
        Test that check args aborts when a config file is entered that doesn't exist
        """
        sys.argv = ["subset_data", "point", "--create-surface", "--cfg-file", "zztop"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError, "Entered default config file does not exist"
        ):
            check_args(self.args)

    def test_check_args_outsurfdat_provided(self):
        """
        Test that check args allows an output surface dataset to be specified
        when create-surface is on
        """
        sys.argv = ["subset_data", "point", "--create-surface", "--out-surface", "outputsurface.nc"]
        self.args = self.parser.parse_args()
        check_args(self.args)
        files = setup_files(self.args, self.defaults, self.cesmroot)
        self.assertEqual(
            files["fsurf_out"],
            "outputsurface.nc",
            "fsurf_out filename not whats expected",
        )

    def test_check_args_outsurfdat_fails_without_create_surface(self):
        """
        Test that check args does not allow an output surface dataset to be specified
        when create-surface is not on
        """
        sys.argv = ["subset_data", "point", "--create-datm", "--out-surface", "outputsurface.nc"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError,
            "out-surface option is given without the --create-surface option",
        ):
            check_args(self.args)

    def test_check_args_fails_for_timeseries_without_correct_surface_year(self):
        """
        Test that check args does not allow landuse-timeseries to be used
        without providing the correct start surface year
        """
        sys.argv = ["subset_data", "point", "--create-landuse", "--create-surface"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError,
            "--surf-year option is NOT set to 1850 and the --create-landuse option",
        ):
            check_args(self.args)

    def test_check_args_fails_for_surf_year_without_surface(self):
        """
        Test that check args does not allow surf_year to be set
        without the create-surface option
        """
        sys.argv = ["subset_data", "point", "--create-datm", "--surf-year", "1850"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError,
            "--surf-year option is set to something besides the default of 2000",
        ):
            check_args(self.args)

    def test_check_args_fails_for_landuse_without_surface(self):
        """
        Test that check args does not allow landuse to be set
        without the create-surface option
        """
        sys.argv = ["subset_data", "point", "--create-landuse"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError,
            "--create-landuse option requires the --create-surface option",
        ):
            check_args(self.args)

    def test_check_args_fails_bad_surface_year(self):
        """
        Test that check args does not allow --surf-year to be bad
        """
        sys.argv = ["subset_data", "point", "--create-surface", "--surf-year", "2305"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError, "--surf-year option can only be set to 1850 or 2000"
        ):
            check_args(self.args)

    def test_check_args_outsurfdat_fails_without_overwrite(self):
        """
        Test that check args does not allow an output surface dataset to be specified
        for an existing dataset without the overwrite option
        """
        outfile = os.path.join(
            os.getcwd(),
            "ctsm/test/testinputs/",
            "surfdata_1x1_mexicocityMEX_hist_16pfts_CMIP6_2000_c231103.nc",
        )
        self.assertTrue(os.path.exists(outfile), str(outfile) + " outfile should exist")

        sys.argv = ["subset_data", "point", "--create-surface", "--out-surface", outfile]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError,
            "out-surface filename exists and the overwrite option was not also selected",
        ):
            check_args(self.args)

    def test_inputdata_setup_files_bad_inputdata_arg(self):
        """
        Test that inputdata directory provided on command line does not exist if it's bad
        """
        check_args(self.args)
        self.args.inputdatadir = "/zztop"
        with self.assertRaisesRegex(SystemExit, "inputdata directory does not exist"):
            setup_files(self.args, self.defaults, self.cesmroot)

    def test_create_user_mods_without_create_mesh(self):
        """
        Test that you can't run create user mods without also doing create_mesh
        """
        sys.argv = ["subset_data", "region", "--create-user-mods", "--create-surface"]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError, "For regional cases, you can not create user_mods"
        ):
            check_args(self.args)

    def test_create_mesh_without_domain(self):
        """
        Test that you can't run create mesh without domain
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-user-mods",
            "--create-surface",
            "--create-mesh",
        ]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentError, "For regional cases, you can not create mesh files"
        ):
            check_args(self.args)

    def test_complex_option_works(self):
        """
        Test that check_args won't flag a set of complex options that is valid
        Do user-mods, surface and landuse-timeseries, as well as DATM, for verbose with crop
        """
        sys.argv = [
            "subset_data",
            "region",
            "--reg",
            "testname",
            "--create-user-mods",
            "--create-surface",
            "--create-landuse",
            "--surf-year",
            "1850",
            "--create-mesh",
            "--create-domain",
            "--create-datm",
            "--verbose",
            "--crop",
        ]
        self.args = self.parser.parse_args()
        check_args(self.args)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
