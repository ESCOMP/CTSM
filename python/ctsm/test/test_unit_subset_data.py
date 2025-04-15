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
from ctsm.subset_data import get_parser, setup_files, check_args, _set_up_regional_case
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=invalid-name,too-many-public-methods


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
        self.args = check_args(self.args)
        files = setup_files(self.args, self.defaults, self.cesmroot)
        self.assertEqual(
            files["fsurf_in"],
            "surfdata_0.9x1.25_hist_2000_16pfts_c240908.nc",
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
        self.args = check_args(self.args)
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
        self.args = check_args(self.args)
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
        self.args = check_args(self.args)
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

    # When CTSM issue #2110 is resolved, this test should be removed.
    def test_subset_region_errors_if_datm(self):
        """
        Test that you can't run subset_data for a region with --create-datm
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-datm",
        ]
        self.args = self.parser.parse_args()
        with self.assertRaisesRegex(
            NotImplementedError, "For regional cases, you can not subset datm data"
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
            # "--create-datm",  # Uncomment this when CTSM issue #2110 is resolved
            "--verbose",
            "--crop",
        ]
        args = self.parser.parse_args()
        args = check_args(args)
        _set_up_regional_case(args)

    def test_region_lon_type_360_ok(self):
        """
        In region mode, test that --lon-type 360 works with valid longitudes
        """
        lon_type = 360
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            str(lon_type),
            "--lon1",
            "320",
            "--lon2",
            "340",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        args = check_args(args)
        self.assertEqual(args.lon1.get(lon_type), 320)
        self.assertEqual(args.lon2.get(lon_type), 340)
        _set_up_regional_case(args)

    def test_region_lon_type_360_toolow(self):
        """
        In region mode, test that --lon-type 360 fails with a longitude value that's below [0, 360]
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            "360",
            "--lon1",
            "-1",
            "--lon2",
            "360",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        with self.assertRaisesRegex(
            ValueError,
            r"lon_in needs to be in the range \[0, 360\]",
        ):
            check_args(args)

    def test_region_lon_type_360_ok_at_360(self):
        """
        In region mode, test that --lon-type 360 works with a longitude value of 360
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            "360",
            "--lon1",
            "320",
            "--lon2",
            "360",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        args = check_args(args)
        _set_up_regional_case(args)

    def test_region_lon_type_360_crosses_pm_errors(self):
        """
        In region mode, test that --lon-type 360 errors if lon range crosses Prime Meridian
        """
        lon1 = 320
        lon2 = 300
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            "360",
            "--lon1",
            str(lon1),
            "--lon2",
            str(lon2),
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()

        expected_err_msg = rf"--lon1 \({lon1}[\.\d]*\) must be < --lon2 \({lon2}[\.\d]*\)"
        with self.assertRaisesRegex(ValueError, expected_err_msg):
            check_args(args)

    def test_region_lon_type_180_crosses_pm_errors(self):
        """
        In region mode, test that --lon-type 180 errors if lon range crosses Prime Meridian
        """
        lon1 = -5
        lon2 = 5
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            "180",
            "--lon1",
            str(lon1),
            "--lon2",
            str(lon2),
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()

        lon1_conv = lon1 % 360
        lon2_conv = lon2 % 360
        expected_err_msg = (
            "After converting to --lon-type 360, "
            + rf"--lon1 \({lon1_conv}[\.\d]*\) must be < --lon2 \({lon2_conv}[\.\d]*\)"
        )
        with self.assertRaisesRegex(ValueError, expected_err_msg):
            check_args(args)

    def test_region_lon_type_180_neg_ok(self):
        """
        In region mode, test that --lon-type 180 works with valid negative longitudes
        """
        lon1 = -87
        lon2 = -24
        lon_type = 180
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            str(lon_type),
            "--lon1",
            str(lon1),
            "--lon2",
            str(lon2),
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        args = check_args(args)
        self.assertEqual(args.lon1.get(lon_type), lon1)
        self.assertEqual(args.lon2.get(lon_type), lon2)
        _set_up_regional_case(args)

    def test_region_lon_type_180_pos_ok(self):
        """
        In region mode, test that --lon-type 180 works with valid positive longitudes
        """
        lon1 = 24
        lon2 = 87
        lon_type = 180
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            str(lon_type),
            "--lon1",
            str(lon1),
            "--lon2",
            str(lon2),
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        args = check_args(args)
        self.assertEqual(args.lon1.get(lon_type), lon1)
        self.assertEqual(args.lon2.get(lon_type), lon2)
        _set_up_regional_case(args)

    def test_point_ambiguous_lon_errors(self):
        """
        In point mode, test that an error is thrown if you give it an ambiguous longitude without
        also giving --lon-type
        """
        sys.argv = [
            "subset_data",
            "point",
            "--create-domain",
            "--verbose",
            "--lat",
            "0",
            "--lon",
            "87",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError,
            "When providing an ambiguous longitude, you must specify --lon-type 180 or 360",
        ):
            check_args(args)

    def test_point_unambiguous_lon_180_ok(self):
        """
        In point mode, test that no error is thrown if an unambiguous longitude is given without
        specifying --lon-type 180
        """
        sys.argv = [
            "subset_data",
            "point",
            "--create-domain",
            "--verbose",
            "--lat",
            "0",
            "--lon",
            "-87",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        check_args(args)

    def test_point_unambiguous_lon_360_ok(self):
        """
        In point mode, test that no error is thrown if an unambiguous longitude is given without
        specifying --lon-type 360
        """
        sys.argv = [
            "subset_data",
            "point",
            "--create-domain",
            "--verbose",
            "--lat",
            "0",
            "--lon",
            "194",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        check_args(args)

    def test_region_ambiguous_lon_errors(self):
        """
        In region mode, test that an error is thrown if you give it one ambiguous longitude without
        also giving --lon-type
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon1",
            "-24",
            "--lon2",
            "87",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError,
            "When providing an ambiguous longitude, you must specify --lon-type 180 or 360",
        ):
            check_args(args)

    def test_region_ambiguous_lons_errors(self):
        """
        In region mode, test that an error is thrown if you give two ambiguous longitudes without
        also giving --lon-type
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon1",
            "24",
            "--lon2",
            "87",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError,
            "When providing an ambiguous longitude, you must specify --lon-type 180 or 360",
        ):
            check_args(args)

    def test_region_unambiguous_lons_180_ok(self):
        """
        In region mode, test that no error is thrown if two unambiguous longitudes are given without
        specifying --lon-type 180
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon1",
            "-87",
            "--lon2",
            "-24",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        args = check_args(args)
        _set_up_regional_case(args)

    def test_region_unambiguous_lons_360_ok(self):
        """
        In region mode, test that no error is thrown if two unambiguous longitudes are given without
        specifying --lon-type 360
        """
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon1",
            "194",
            "--lon2",
            "287",
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        args = check_args(args)
        _set_up_regional_case(args)

    def test_region_lon_type_180_ok_at_180(self):
        """
        In region mode, test that --lon-type 180 passes at lon 180
        """
        lon1 = 24
        lon2 = 180
        lon_type = 180
        sys.argv = [
            "subset_data",
            "region",
            "--create-domain",
            "--verbose",
            "--lat1",
            "0",
            "--lat2",
            "40",
            "--lon-type",
            str(lon_type),
            "--lon1",
            str(lon1),
            "--lon2",
            str(lon2),
        ]
        self.parser = get_parser()
        args = self.parser.parse_args()
        args = check_args(args)
        self.assertEqual(args.lon1.get(lon_type), lon1)
        self.assertEqual(args.lon2.get(lon_type), lon2)
        _set_up_regional_case(args)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
