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

    def _check_result_file_matches_expected(self, expected_output_files):
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
                inspect.stack()[1][3],  # Name of calling function (i.e., test name)
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
        self.assertTrue(self._check_result_file_matches_expected(expected_output_files))

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


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
