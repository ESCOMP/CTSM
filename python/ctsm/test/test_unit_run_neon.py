#!/usr/bin/env python3
"""
Unit tests for run_neon

You can run this by:
    python -m unittest test_unit_run_neon.py
"""

import unittest
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.run_neon import check_neon_listing

# pylint: disable=invalid-name


class TestRunNeon(unittest.TestCase):
    """
    Basic class for testing run_neon.py.
    """

    def test_check_neon_listing(self):
        """
        Test that neon listing is available for valid sites
        """
        valid_neon_sites = ["ABBY", "BART"]
        available_list = check_neon_listing(valid_neon_sites)
        self.assertEqual(
            available_list[0].name, "ABBY", "available list of actual sites not as expected"
        )
        self.assertEqual(
            available_list[1].name, "BART", "available list of actual sites not as expected"
        )

    def test_check_neon_listing_misspelled(self):
        """
        Test that neon listing is not available for invalid sites
        """
        valid_neon_sites = ["INVALID_SITE1", "INVALID_SITE2"]
        available_list = check_neon_listing(valid_neon_sites)
        self.assertEqual(
            available_list, [], "available list of incorrect dummy site not as expected"
        )


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
