#!/usr/bin/env python3
"""
Unit tests for tower_site that aren't specific to a single tower type

You can run this by:
    python -m unittest test_unit_tower_site.py
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
from ctsm.site_and_regional.tower_site import TowerSite, ALLOWED_SITE_TYPES

# pylint: disable=invalid-name


class Test_tower_site(unittest.TestCase):
    """
    Basic class for testing tower_site.py for things not specific to any tower type.
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

    def test_disallowed_tower_types(self):
        """
        Test that an invalid tower type throws an error
        """
        # NeonSite parameters (not that it matters)
        name = "ABBY"
        start_year = 2020
        end_year = 2021
        start_month = 1
        end_month = 12
        # finidat = None
        finidat = "dummy_finidat"

        # try to create site
        invalid_type = "INVALIDTYPE"
        err_msg = (
            f"tower_type '{invalid_type}' not allowed. "
            + f"Choose from: {','.join(ALLOWED_SITE_TYPES)}"
        )
        with self.assertRaisesRegex(ValueError, err_msg):
            TowerSite(invalid_type, name, start_year, end_year, start_month, end_month, finidat)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
