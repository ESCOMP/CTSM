#!/usr/bin/env python3
"""
Unit tests for subset_data
"""

import unittest
import argparse

from ctsm.subset_data import plon_type
from ctsm import unit_testing

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestSubsetData(unittest.TestCase):
    """
    Tests for subset_data
    """
    def test_plonType_positive(self):
        """
        Test plot_type with 30
        """
        result = plon_type(30)
        self.assertEqual(result, 30.0)

    def test_plonType_negative(self):
        """
        Test plot_type with -30
        """
        result = plon_type(-30)
        self.assertEqual(result, 330.0)

    def test_plonType_outOfBounds(self):
        """
        Test plot_type with 361
        """
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "Latitude.*should be between"
        ):
            _ = plon_type(361)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
