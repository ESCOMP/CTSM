#!/usr/bin/env python3
"""
Unit tests for run_tower

You can run this by:
    python -m unittest test_unit_run_tower.py
"""

import unittest
from unittest import mock
import tempfile
import shutil
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
# Import add_cime_to_path for its side-effect of populating sys.path with CIME modules
import ctsm.add_cime_to_path  # pylint: disable=unused-import
from ctsm import unit_testing
from ctsm.site_and_regional.run_tower import check_neon_listing
from ctsm.site_and_regional.tower_site import TowerSite

# pylint: disable=invalid-name


class TestRunTower(unittest.TestCase):
    """
    Basic class for testing run_tower.py.
    """

    def setUp(self):
        """
        Make /_tempdir for use by these tests.
        """
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_check_neon_listing(self):
        """
        Test that neon listing is available for valid sites
        """
        valid_neon_sites = ["ABBY", "BART"]
        previous_dir = os.getcwd()
        os.chdir(self._tempdir)  # cd to tempdir
        available_list = check_neon_listing(valid_neon_sites)
        self.assertEqual(
            available_list[0].name, "ABBY", "available list of actual sites not as expected"
        )
        self.assertEqual(
            available_list[1].name, "BART", "available list of actual sites not as expected"
        )
        # change to previous dir once listing.csv file is created in tempdir and test complete
        os.chdir(previous_dir)

    def test_check_neon_listing_misspelled(self):
        """
        Test that neon listing is not available for invalid sites
        """
        valid_neon_sites = ["INVALID_SITE1", "INVALID_SITE2"]
        previous_dir = os.getcwd()
        os.chdir(self._tempdir)  # cd to tempdir
        available_list = check_neon_listing(valid_neon_sites)
        self.assertEqual(
            available_list, [], "available list of incorrect dummy site not as expected"
        )
        # change to previous dir once listing.csv file is created in tempdir and test complete
        os.chdir(previous_dir)

    def test_towerSite_cesmCheckoutLayout(self):
        """Test TowerSite sets cesmroot to top-level CESM root and user_mods_dirs to CTSM root
        in CESM checkout layout."""
        ctsm_path = os.path.join(self._tempdir, "components", "clm")
        cesm_path = self._tempdir
        os.makedirs(ctsm_path)
        os.makedirs(os.path.join(cesm_path, "cime"))

        # Create mock usermods_dirs under CTSM root
        usermod_dir = os.path.join(ctsm_path, "cime_config", "usermods_dirs", "clm", "NEON", "ABBY")
        os.makedirs(usermod_dir)

        with mock.patch(
            "ctsm.site_and_regional.tower_site.path_to_top_root", return_value=cesm_path
        ), mock.patch(
            "ctsm.site_and_regional.tower_site.path_to_ctsm_root", return_value=ctsm_path
        ):
            site = TowerSite(
                "NEON",
                name="ABBY",
                start_year=2018,
                end_year=2018,
                start_month=1,
                end_month=12,
                finidat=None,
            )

        self.assertEqual(site.cesmroot, cesm_path)
        self.assertEqual(site.user_mods_dirs[0], usermod_dir)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
