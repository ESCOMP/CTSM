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

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm import subset_data
from ctsm.utils import find_one_file_matching_pattern


class TestSubsetDataSys(unittest.TestCase):
    """
    Basic class for testing subset_data.py.
    """

    def setUp(self):
        self.temp_dir_out = tempfile.TemporaryDirectory()
        self.temp_dir_umd = tempfile.TemporaryDirectory()
        self.inputdata_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)

    def tearDown(self):
        self.temp_dir_out.cleanup()
        self.temp_dir_umd.cleanup()

    def _check_result_file_matches_expected(self, expected_output_files, caller_n):
        """
        Loop through a list of output files, making sure they match what we expect.
        """
        all_files_present_and_match = True
        for basename in expected_output_files:
            result_file = os.path.join(self.temp_dir_out.name, basename)
            result_file = find_one_file_matching_pattern(result_file)
            expected_file = os.path.join(
                os.path.dirname(__file__),
                "testinputs",
                "expected_result_files",
                inspect.stack()[caller_n][3],  # Name of calling function (i.e., test name)
                basename,
            )
            expected_file = find_one_file_matching_pattern(expected_file)
            ds_result = xr.open_dataset(result_file)
            ds_expected = xr.open_dataset(expected_file)
            if not ds_result.equals(ds_expected):
                print("Result differs from expected: " + basename)
                print(ds_result)
                print(ds_expected)
                all_files_present_and_match = False
        return all_files_present_and_match

    def test_subset_data_reg_amazon(self):
        """
        Test subset_data for Amazon region
        """
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
            "-12",
            "--lat2",
            "-7",
            "--lon1",
            "291",
            "--lon2",
            "299",
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
            "--inputdata-dir",
            self.inputdata_dir,
            "--cfg-file",
            cfg_file,
            "--overwrite",
        ]
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        expected_output_files = [
            f"domain.lnd.5x5pt-amazon_navy_TMP_c{daystr}_ESMF_UNSTRUCTURED_MESH.nc",
            f"domain.lnd.5x5pt-amazon_navy_TMP_c{daystr}.nc",
            f"surfdata_TMP_amazon_hist_16pfts_CMIP6_2000_c{daystr}.nc",
        ]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 1))

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

    def _do_test_subset_data_pt_surface(self, lon):
        """
        Given a longitude, test subset_data point --create-surface
        """
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
            "-12",
            "--lon",
            str(lon),
            "--site",
            "TMP",
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
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        expected_output_files = [
            f"surfdata_TMP_amazon_hist_16pfts_CMIP6_2000_c{daystr}.nc",
        ]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 2))

    def test_subset_data_pt_surface_amazon_type360(self):
        """
        Test subset_data --create-surface for Amazon point with longitude type 360
        """
        self._do_test_subset_data_pt_surface(291)

    def test_subset_data_pt_surface_amazon_type180(self):
        """
        Test subset_data --create-surface for Amazon point with longitude type 180
        """
        self._do_test_subset_data_pt_surface(-69)

    def _do_test_subset_data_pt_landuse(self, lon):
        """
        Given a longitude, test subset_data point --create-landuse
        """
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
            "-12",
            "--lon",
            str(lon),
            "--site",
            "TMP",
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
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        expected_output_files = [
            f"surfdata_TMP_amazon_hist_1850_78pfts_c{daystr}.nc",
            f"landuse.timeseries_TMP_amazon_hist_1850-1853_78pfts_c{daystr}.nc",
        ]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 2))

    def test_subset_data_pt_landuse_amazon_type360(self):
        """
        Test subset_data --create-landuse for Amazon point with longitude type 360
        """
        self._do_test_subset_data_pt_landuse(291)

    def test_subset_data_pt_landuse_amazon_type180(self):
        """
        Test subset_data --create-landuse for Amazon point with longitude type 180
        """
        self._do_test_subset_data_pt_landuse(-69)

    def _do_test_subset_data_pt_datm(self, lon):
        """
        Given a longitude, test subset_data point --create-datm
        """
        start_year = 1986
        end_year = 1988
        sitename = "TMP"
        outdir = self.temp_dir_out.name
        sys.argv = [
            "subset_data",
            "point",
            "--lat",
            "-12",
            "--lon",
            str(lon),
            "--site",
            sitename,
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
        subset_data.main()

        # Loop through all the output files, making sure they match what we expect.
        daystr = "[0-9][0-9][0-9][0-9][0-9][0-9]"  # 6-digit day code, yymmdd
        expected_output_files = [
            f"domain.crujra_v2.3_0.5x0.5_{sitename}_c{daystr}.nc",
        ]
        for year in list(range(start_year, end_year + 1)):
            for forcing in ["Solr", "Prec", "TPQWL"]:
                expected_output_files.append(
                    f"clmforc.CRUJRAv2.5_0.5x0.5.{forcing}.{sitename}.{year}.nc"
                )
        expected_output_files = [os.path.join("datmdata", x) for x in expected_output_files]
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files, 2))

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


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
