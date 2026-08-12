#!/usr/bin/env python3

"""System tests for regrid_ggcmi_shdates.py"""

import os

import unittest
import tempfile
import shutil
import sys

import xarray as xr
import numpy as np

# -- add python/ctsm  to path (needed if we want to run test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.crop_calendars.regrid_ggcmi_shdates import regrid_ggcmi_shdates
from ctsm.crop_calendars.regrid_ggcmi_shdates import regrid_ggcmi_shdates_arg_process

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestRegridGgcmiShdates(unittest.TestCase):
    """System tests for regrid_ggcmi_shdates"""

    def setUp(self):
        # Where in the /testinputs directory are the raw crop calendar file(s)?
        testinputs_path = os.path.join(path_to_ctsm_root(), "python", "ctsm", "test", "testinputs")
        testinputs_cc_path = os.path.join(testinputs_path, "cropcals")
        self._testinputs_cc_path = testinputs_cc_path

        # Make /_tempdir for use by these tests.
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()

        # Obtain path for the directory being created in /_tempdir
        self._regridded_cropcals = os.path.join(self._tempdir, "regridded_cropcals")

        # What extension do the raw crop calendar file(s) have?
        self._extension = ".nc4"

        # Which crop(s) should we test? (comma-separated string)
        self._crop_list = "swh_rf"

        # What is the complete set of input arguments (including script name)?
        regrid_template_file = os.path.join(
            testinputs_path, "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031.nc"
        )
        self._function_call_list = [
            "regrid_ggcmi_shdates",
            "-i",
            testinputs_cc_path,
            "-x",
            ".nc4",
            "-o",
            self._regridded_cropcals,
            "-rr",
            "5x5amazon",
            "-rt",
            regrid_template_file,
            "--crop-list",
            "swh_rf",
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_regrid_ggcmi_shdates(self):
        """
        Tests regrid_ggcmi_shdates
        """

        # Call script
        sys.argv = self._function_call_list
        args = regrid_ggcmi_shdates_arg_process()
        regrid_ggcmi_shdates(
            args.regrid_resolution,
            args.regrid_template_file,
            args.regrid_input_directory,
            args.regrid_output_directory,
            args.regrid_extension,
            args.crop_list,
        )

        # Read output file
        regrid_out_file = os.path.join(
            self._regridded_cropcals,
            "swh_rf_ggcmi_crop_calendar_phase3_v1.01_nninterp-5x5amazon.nc4",
        )
        regrid_out_ds = xr.open_dataset(regrid_out_file)

        # Check sowing dates
        expected_sow_dates = np.array(
            [
                [120, 120, 120, 120, 120],
                [120, 120, 120, 120, 120],
                [120, 120, 120, 120, 120],
                [330, 335, 335, 120, 120],
                [325, 335, 335, 335, 120],
            ]
        )
        np.testing.assert_array_equal(expected_sow_dates, regrid_out_ds["planting_day"].values)

        # Check maturity dates
        expected_mat_dates = np.array(
            [
                [221, 221, 221, 221, 221],
                [221, 221, 221, 221, 221],
                [221, 221, 221, 221, 221],
                [153, 128, 128, 221, 221],
                [163, 128, 128, 128, 221],
            ]
        )
        np.testing.assert_array_equal(expected_mat_dates, regrid_out_ds["maturity_day"].values)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
