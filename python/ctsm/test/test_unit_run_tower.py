#!/usr/bin/env python3
"""
Unit tests for run_tower

You can run this by:
    python -m unittest test_unit_run_tower.py
"""

import unittest
import tempfile
import shutil
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.run_tower import check_neon_listing

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


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
