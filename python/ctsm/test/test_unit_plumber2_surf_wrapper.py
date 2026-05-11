#!/usr/bin/env python3
"""
Unit tests for plumber2_surf_wrapper

You can run this by:
    python -m unittest test_unit_plumber2_surf_wrapper.py
"""

import unittest
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.plumber2_surf_wrapper import get_args

# pylint: disable=invalid-name


class TestPlumber2SurfWrapper(unittest.TestCase):
    """
    Basic class for testing plumber2_surf_wrapper.py.
    """

    def setUp(self):
        sys.argv = ["subset_data"]  # Could actually be anything

    def test_parser_default_csv_exists(self):
        """
        Test that default PLUMBER2 sites CSV file exists
        """

        args = get_args()
        self.assertTrue(os.path.exists(args.plumber2_sites_csv))

    def test_parser_custom_csv(self):
        """
        Test that script accepts custom CSV file path
        """

        custom_path = "path/to/custom.csv"
        sys.argv += ["--plumber2-sites-csv", custom_path]
        args = get_args()
        self.assertEqual(args.plumber2_sites_csv, custom_path)

    def test_parser_verbose_false_default(self):
        """
        Test that script is not verbose by default
        """

        args = get_args()
        self.assertFalse(args.verbose)

    def test_parser_verbose_true(self):
        """
        Test that --verbose sets verbose to True
        """

        sys.argv += ["--verbose"]
        args = get_args()
        self.assertTrue(args.verbose)

    def test_parser_78pft_false_default(self):
        """
        Test that script does not use 78pft mode by default
        """

        args = get_args()
        self.assertFalse(args.use_managed_crops)

    def test_parser_78pft_true(self):
        """
        Test that --crop sets use_managed_crops to True
        """

        sys.argv += ["--crop"]
        args = get_args()
        self.assertTrue(args.use_managed_crops)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
