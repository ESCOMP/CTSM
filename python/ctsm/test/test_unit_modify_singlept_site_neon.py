#!/usr/bin/env python3
"""
Unit tests for modify_singlept_site_neon

You can run this by:
    python -m unittest test_unit_modify_singlept_site_neon.py
"""

import os
import shutil
import sys
import tempfile
import unittest
from datetime import date
import xarray as xr

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.modify_singlept_site_neon import (
    get_neon,
    find_surffile,
    update_metadata,
    update_time_tag,
    check_neon_time,
)

# pylint: disable=invalid-name


class TestModifySingleptSiteNeon(unittest.TestCase):
    """
    Basic class for testing modify_singlept_site_neon.py.
    """

    def setUp(self):
        """
        Make /_tempdir for use by these tests.
        Check tempdir for history files
        """
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_get_neon(self):
        """
        Test to see if neon data for valid site name is found
        """
        site_name = "ABBY"
        neon_dir = self._tempdir
        file = get_neon(neon_dir, site_name)
        self.assertEqual(file.split("/")[-1][:4], "ABBY", "CSV file did not download as expected")

    def test_get_neon_false_site(self):
        """
        Test to see if neon data for invalid site name is found
        """
        site_name = "INVALID_SITE"
        neon_dir = self._tempdir
        with self.assertRaises(SystemExit):
            get_neon(neon_dir, site_name)

    def test_find_surffile(self):
        """
        Test that surface file does not exist in tempdir and raises system exit error
        """
        surf_dir = self._tempdir
        site_name = "BART"
        pft_16 = True
        with self.assertRaises(SystemExit):
            find_surffile(surf_dir, site_name, pft_16)

    def test_find_soil_structure(self):
        """
        Test to ensure that correct attributes are found for find_soil_structure.
        soil_texture_raw_data_file_name should be found, and test should go through sysexit.
        """
        surf_file_name = "surfdata_1x1_mexicocityMEX_hist_16pfts_Irrig_CMIP6_simyr2000_c221206.nc"
        surf_file = os.path.join(
            path_to_ctsm_root(),
            "python/ctsm/test/testinputs/",
            surf_file_name,
        )
        f1 = xr.open_dataset(surf_file)
        self.assertEqual(
            f1.attrs["Soil_texture_raw_data_file_name"],
            "mksrf_soitex.10level.c010119.nc",
            "did not retrieve expected surface soil texture filename from surf file",
        )

    def test_update_metadata(self):
        """
        Test to ensure that the file was updated today.
        """
        surf_file = "surfdata_1x1_mexicocityMEX_hist_16pfts_Irrig_CMIP6_simyr2000_c221206.nc"
        neon_file = "dummy_neon_file.nc"
        zb_flag = True
        f1 = xr.open_dataset(
            os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs/") + surf_file
        )
        f2 = update_metadata(f1, surf_file, neon_file, zb_flag)
        today = date.today()
        today_string = today.strftime("%Y-%m-%d")
        self.assertEqual(f2.attrs["Updated_on"], today_string, "File was not updated as expected")

    def test_update_time_tag(self):
        """
        Test that file ending is updated
        """
        self.assertEqual(
            update_time_tag("test_YYMMDD.nc")[-9:-3],
            date.today().strftime("%y%m%d"),
            "File ending not as expected",
        )

    def test_check_neon_time(self):
        """
        Test that dictionary containing last modified information is correctly downloaded
        """
        previous_dir = os.getcwd()
        os.chdir(self._tempdir)  # cd to tempdir
        last_abby_download = check_neon_time()[
            "https://storage.neonscience.org/neon-ncar/NEON/surf_files/v1/ABBY_surfaceData.csv"
        ]
        self.assertEqual(
            len(last_abby_download),
            19,
            "last ABBY download has unexpected date format or does not exist",
        )
        # Note: this checks that data is not pulled from before 2021;
        # we may want to update this occassionally,
        # but in any case it confirms that the oldest data is not found
        self.assertGreater(
            int(last_abby_download[:4]), 2021, "ABBY download is older than expected"
        )
        # change back to previous dir once listing.csv file is created in tempdir and test complete
        os.chdir(previous_dir)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
