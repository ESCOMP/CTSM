#!/usr/bin/env python3
"""
System tests for subset_data

You can run this by:
    python -m unittest test_sys_subset_data.py
"""

import unittest
import os
import sys
import tempfile
import inspect
import xarray as xr
from CIME.scripts.create_newcase import _main_func as create_newcase  # pylint: disable=import-error

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm import subset_data
from ctsm.utils import find_one_file_matching_pattern


def _get_sitename_str_point(include_sitename, sitename, lon, lat):
    """
    Given a site, return the string to use in output filenames
    """
    if include_sitename:
        sitename_str = sitename
    else:
        sitename_str = f"{float(lon)}_{float(lat)}"
    return sitename_str


class TestSubsetDataSys(unittest.TestCase):
    """
    Basic class for testing subset_data.py.
    """

    def setUp(self):
        self.previous_dir = os.getcwd()
        self.temp_dir_out = tempfile.TemporaryDirectory()
        self.temp_dir_umd = tempfile.TemporaryDirectory()
        self.temp_dir_caseparent = tempfile.TemporaryDirectory()
        self.inputdata_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)

    def tearDown(self):
        self.temp_dir_out.cleanup()
        self.temp_dir_umd.cleanup()
        os.chdir(self.previous_dir)

    def _check_create_newcase(self):
        """
        Check that you can call create_newcase using the usermods from subset_data
        """
        case_dir = os.path.join(self.temp_dir_caseparent.name, "case")
        sys.argv = [
            "create_newcase",
            "--case",
            case_dir,
            "--res",
            "CLM_USRDAT",
            "--compset",
            "I2000Clm60Bgc",
            "--run-unsupported",
            "--user-mods-dir",
            self.temp_dir_umd.name,
        ]
        create_newcase()

    def _check_result_file_matches_expected(self, expected_output_files, caller_n):
        """
        Loop through a list of output files, making sure they match what we expect.

        caller_n should be an integer giving the number of levels above this function you need to
        traverse before you hit the actual test name. If the test is calling this function directly,
        caller_n = 1. If the test is calling a function that calls this function, caller_n = 2. Etc.
        """
        all_files_present_and_match = True
        result_file_found = True
        expected_file_found = True
        for basename in expected_output_files:

            # Check whether result (output) file exists. If not, note it but continue.
            result_file = os.path.join(self.temp_dir_out.name, basename)
            try:
                result_file = find_one_file_matching_pattern(result_file)
            except FileNotFoundError:
                result_file_found = False

            # Check whether expected file exists. If not, note it but continue.
            expected_file = os.path.join(
                os.path.dirname(__file__),
                "testinputs",
                "expected_result_files",
                inspect.stack()[caller_n][3],  # Name of calling function (i.e., test name)
                basename,
            )
            try:
                expected_file = find_one_file_matching_pattern(expected_file)
            except FileNotFoundError:
                expected_file_found = False

            # Raise an AssertionError if either file was not found
            if not (result_file_found and expected_file_found):
                msg = ""
                if not result_file_found:
                    this_dir = os.path.dirname(result_file)
                    msg += f"\nResult file '{result_file}' not found. "
                    msg += f"Contents of directory '{this_dir}':\n\t"
                    msg += "\n\t".join(os.listdir(this_dir))
                if not expected_file_found:
                    this_dir = os.path.dirname(expected_file)
                    msg += f"\nExpected file '{expected_file}' not found. "
                    msg += f"Contents of directory '{this_dir}':\n\t"
                    msg += "\n\t".join(os.listdir(this_dir))
                raise AssertionError(msg)

            # Compare the two files
            ds_result = xr.open_dataset(result_file)
            ds_expected = xr.open_dataset(expected_file)
            if not ds_result.equals(ds_expected):
                print("Result differs from expected: " + basename)
                print(ds_result)
                print(ds_expected)
                all_files_present_and_match = False
        return all_files_present_and_match

    def _do_test_subset_data_reg_amazon(self, include_regname=True):
        """
        Convenience function for multiple tests of subset_data region for the Amazon
        """
        regname = "TMP"
        lat1 = -12
        lat2 = -7
        lon1 = 291
        lon2 = 299
        cfg_file = os.path.join(
            self.inputdata_dir,
            "ctsm",
            "test",
            "testinputs",
            "subset_data_amazon.cfg",
        )
        print(cfg_file)
        sys.argv = [
            "subset_data",
            "region",
            "--lat1",
            str(lat1),
            "--lat2",
            str(lat2),
            "--lon1",
            str(lon1),
            "--lon2",
            str(lon2),
            "--create-mesh",
            "--create-domain",
            "--create-surface",
            "--surf-year",
            "2000",
            "--create-user-mods",
            "--outdir",
            self.temp_dir_out.name,
            "--user-mods-dir",
            self.temp_dir_umd.name,
            "--inputdata-dir",
            self.inputdata_dir,
            "--cfg-file",
            cfg_file,
            "--overwrite",
        ]
        if include_regname:
            sys.argv += ["--reg", regname]
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        if include_regname:
            regname_str = regname
        else:
            regname_str = f"{float(lon1)}-{float(lon2)}_{float(lat1)}-{float(lat2)}"
        expected_output_files = [
            f"domain.lnd.5x5pt-amazon_navy_{regname_str}_c{daystr}_ESMF_UNSTRUCTURED_MESH.nc",
            f"domain.lnd.5x5pt-amazon_navy_{regname_str}_c{daystr}.nc",
            f"surfdata_{regname_str}_amazon_hist_16pfts_CMIP6_2000_c{daystr}.nc",
        ]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 2))

        # Check that create_newcase works
        # SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        self._check_create_newcase()

    def test_subset_data_reg_amazon(self):
        """
        Test subset_data for Amazon region
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_reg_amazon()

    def test_subset_data_reg_amazon_noregname(self):
        """
        Test subset_data for Amazon region
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_reg_amazon(include_regname=False)

    def test_subset_data_reg_infile_detect360(self):
        """
        Test subset_data for region with ambiguous longitudes. We specify the longitude type for
        lon1 and lon2 but not for the input data files. This should still work as long as the input
        data file longitude type is detectable and matches --lon-type.
        """
        sys.argv = [
            "subset_data",
            "region",
            "--lat1",
            "-12",
            "--lat2",
            "-7",
            "--lon1",
            "15",
            "--lon2",
            "23",
            "--lon-type",
            "360",
            "--reg",
            "TMP",
            "--create-mesh",
            "--create-domain",
            "--create-surface",
            "--surf-year",
            "2000",
            "--create-user-mods",
            "--outdir",
            self.temp_dir_out.name,
            "--user-mods-dir",
            self.temp_dir_umd.name,
            "--overwrite",
        ]
        subset_data.main()

    def test_subset_data_reg_infile_detect180_error(self):
        """
        Specifying --lon-type 180 but an input file of type 360 should error
        """
        sys.argv = [
            "subset_data",
            "region",
            "--lat1",
            "-12",
            "--lat2",
            "-7",
            "--lon1",
            "15",
            "--lon2",
            "23",
            "--lon-type",
            "180",
            "--reg",
            "TMP",
            "--create-mesh",
            "--create-domain",
            "--create-surface",
            "--surf-year",
            "2000",
            "--create-user-mods",
            "--outdir",
            self.temp_dir_out.name,
            "--user-mods-dir",
            self.temp_dir_umd.name,
            "--overwrite",
        ]
        with self.assertRaisesRegex(
            RuntimeError, r"File lon type \(360\) doesn't match boundary lon type \(180\)"
        ):
            subset_data.main()

    def _do_test_subset_data_pt_surface(self, lon, include_sitename=True):
        """
        Given a longitude, test subset_data point --create-surface
        """
        lat = -12
        cfg_file = os.path.join(
            self.inputdata_dir,
            "ctsm",
            "test",
            "testinputs",
            "subset_data_amazon.cfg",
        )
        print(cfg_file)
        sys.argv = [
            "subset_data",
            "point",
            "--lat",
            str(lat),
            "--lon",
            str(lon),
            "--create-domain",
            "--create-surface",
            "--surf-year",
            "2000",
            "--create-user-mods",
            "--outdir",
            self.temp_dir_out.name,
            "--user-mods-dir",
            self.temp_dir_umd.name,
            "--inputdata-dir",
            self.inputdata_dir,
            "--cfg-file",
            cfg_file,
            "--overwrite",
        ]
        sitename = "TMP"
        if include_sitename:
            sys.argv += ["--site", sitename]
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        sitename_str = _get_sitename_str_point(include_sitename, sitename, lon, lat)
        expected_output_files = [
            f"surfdata_{sitename_str}_amazon_hist_16pfts_CMIP6_2000_c{daystr}.nc",
        ]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 2))

        # Check that create_newcase works
        # SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        self._check_create_newcase()

    def test_subset_data_pt_surface_amazon_type360(self):
        """
        Test subset_data --create-surface for Amazon point with longitude type 360
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_pt_surface(291)

    def test_subset_data_pt_surface_amazon_type180(self):
        """
        Test subset_data --create-surface for Amazon point with longitude type 180
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_pt_surface(-69)

    def test_subset_data_pt_surface_amazon_type180_nositename(self):
        """
        Test subset_data --create-surface for Amazon point with longitude type 180
        without specifying a site name
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_pt_surface(-69, include_sitename=False)

    def _do_test_subset_data_pt_landuse(self, lon, include_sitename=True):
        """
        Given a longitude, test subset_data point --create-landuse
        """
        lat = -12
        sitename = "TMP"
        cfg_file = os.path.join(
            self.inputdata_dir,
            "ctsm",
            "test",
            "testinputs",
            "subset_data_amazon_1850.cfg",
        )
        print(cfg_file)
        sys.argv = [
            "subset_data",
            "point",
            "--lat",
            str(lat),
            "--lon",
            str(lon),
            "--create-domain",
            "--create-surface",
            "--surf-year",
            "1850",
            "--create-landuse",
            "--create-user-mods",
            "--outdir",
            self.temp_dir_out.name,
            "--user-mods-dir",
            self.temp_dir_umd.name,
            "--inputdata-dir",
            self.inputdata_dir,
            "--cfg-file",
            cfg_file,
            "--overwrite",
        ]
        if include_sitename:
            sys.argv += ["--site", sitename]
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        sitename_str = _get_sitename_str_point(include_sitename, sitename, lon, lat)
        expected_output_files = [
            f"surfdata_{sitename_str}_amazon_hist_1850_78pfts_c{daystr}.nc",
            f"landuse.timeseries_{sitename_str}_amazon_hist_1850-1853_78pfts_c{daystr}.nc",
        ]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 2))

        # Check that create_newcase works
        # SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        self._check_create_newcase()

    def test_subset_data_pt_landuse_amazon_type360(self):
        """
        Test subset_data --create-landuse for Amazon point with longitude type 360
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_pt_landuse(291)

    def test_subset_data_pt_landuse_amazon_type360_nositename(self):
        """
        Test subset_data --create-landuse for Amazon point with longitude type 360 and no site name
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_pt_landuse(291, include_sitename=False)

    def test_subset_data_pt_landuse_amazon_type180(self):
        """
        Test subset_data --create-landuse for Amazon point with longitude type 180
        SHOULD WORK ONLY ON CESM-SUPPORTED MACHINES
        """
        self._do_test_subset_data_pt_landuse(-69)

    def _do_test_subset_data_pt_datm(self, lon, include_sitename=True):
        """
        Given a longitude, test subset_data point --create-datm
        """
        start_year = 1986
        end_year = 1988
        sitename = "TMP"
        lat = -12
        outdir = self.temp_dir_out.name
        sys.argv = [
            "subset_data",
            "point",
            "--lat",
            str(lat),
            "--lon",
            str(lon),
            "--create-datm",
            "--datm-syr",
            str(start_year),
            "--datm-eyr",
            str(end_year),
            "--create-user-mods",
            "--outdir",
            outdir,
            "--user-mods-dir",
            self.temp_dir_umd.name,
            "--overwrite",
        ]
        if include_sitename:
            sys.argv += ["--site", sitename]
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        sitename_str = _get_sitename_str_point(include_sitename, sitename, lon, lat)
        expected_output_files = [
            f"domain.crujra_v2.3_0.5x0.5_{sitename_str}_c{daystr}.nc",
        ]
        for year in list(range(start_year, end_year + 1)):
            for forcing in ["Solr", "Prec", "TPQWL"]:
                expected_output_files.append(
                    f"clmforc.CRUJRAv2.5_0.5x0.5.{forcing}.{sitename_str}.{year}.nc"
                )
        expected_output_files = [os.path.join("datmdata", x) for x in expected_output_files]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 2))

        # Check that create_newcase works
        self._check_create_newcase()

    def test_subset_data_pt_datm_amazon_type360(self):
        """
        Test subset_data --create-datm for Amazon point with longitude type 360
        FOR NOW CAN ONLY BE RUN ON DERECHO/CASPER
        """
        self._do_test_subset_data_pt_datm(291)

    def test_subset_data_pt_datm_amazon_type180(self):
        """
        Test subset_data --create-datm for Amazon point with longitude type 180
        FOR NOW CAN ONLY BE RUN ON DERECHO/CASPER
        """
        self._do_test_subset_data_pt_datm(-69)

    def test_subset_data_pt_datm_amazon_type180_nositename(self):
        """
        Test subset_data --create-datm for Amazon point with longitude type 180 without providing
        site name.
        FOR NOW CAN ONLY BE RUN ON DERECHO/CASPER
        """
        self._do_test_subset_data_pt_datm(-69, include_sitename=False)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
