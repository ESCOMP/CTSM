#!/usr/bin/env python3
"""
Unit tests for NeonSite

You can run this by:
    python -m unittest test_unit_neon_site.py
"""

import unittest
import tempfile
import shutil
import os
import glob
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.neon_site import NeonSite

# pylint: disable=invalid-name


class TestNeonSite(unittest.TestCase):
    """
    Basic class for testing neon_site.py.
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

    def test_modify_user_nl_transient(self):
        """
        Test that modify_user_nl is correctly adding lines to namelist for transient cases
        """
        # NeonSite parameters:
        name = "ABBY"
        start_year = 2020
        end_year = 2021
        start_month = 1
        end_month = 12
        # finidat = None
        finidat = "dummy_finidat"

        # modify_user_nl parameters:
        case_root = self._tempdir
        run_type = "transient"
        rundir = ""

        # create NeonSite object and update namelist
        NeonSite(
            name=name,
            start_year=start_year,
            end_year=end_year,
            start_month=start_month,
            end_month=end_month,
            finidat=finidat,
        ).modify_user_nl(case_root, run_type, rundir)

        # gather file contents for test
        new_nl_file = open(glob.glob(case_root + "/*")[0], "r")
        lines_read = new_nl_file.readlines()[0]
        new_nl_file.close()

        # assertion
        self.assertEqual(
            lines_read,
            "finidat = '/inputdata/lnd/ctsm/initdata/dummy_finidat'\n",
            "transient case has unexpected nl",
        )

    def test_modify_user_nl_ad(self):
        """
        Test that modify_user_nl is correctly adding lines to namelist for ad cases
        """
        # NeonSite parameters:
        name = "ABBY"
        start_year = 2020
        end_year = 2021
        start_month = 1
        end_month = 12
        # finidat = None
        finidat = "dummy_finidat"

        # modify_user_nl parameters:
        case_root = self._tempdir
        run_type = "ad"
        rundir = ""

        # create NeonSite object and update namelist
        NeonSite(
            name=name,
            start_year=start_year,
            end_year=end_year,
            start_month=start_month,
            end_month=end_month,
            finidat=finidat,
        ).modify_user_nl(case_root, run_type, rundir)

        # gather file contents for test
        new_nl_file = open(glob.glob(case_root + "/*")[0], "r")
        lines_read = new_nl_file.readlines()[1]
        new_nl_file.close()

        # assertion
        self.assertEqual(lines_read, "hist_mfilt = 20\n", "ad case has unexpected nl")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
