#!/usr/bin/env python3

"""Unit tests for the iso functions in utils"""

import unittest

from ctsm import unit_testing
from ctsm.utils import parse_isoduration, get_isosplit

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestIsoUtils(unittest.TestCase):
    """Tests of iso functions in utils"""

    def test_iso_split_for_Year(self):
        """
        Tests the get_isosplit function for a strings with Years
        """
        iso_string = "0Y"
        self.assertEqual(get_isosplit(iso_string, "Y"), ("0", ""))
        iso_string = "1Y"
        self.assertEqual(get_isosplit(iso_string, "Y"), ("1", ""))
        iso_string = "4Y"
        self.assertEqual(get_isosplit(iso_string, "Y"), ("4", ""))
        iso_string = "100Y"
        self.assertEqual(get_isosplit(iso_string, "Y"), ("100", ""))
        iso_string = "999999Y"
        self.assertEqual(get_isosplit(iso_string, "Y"), ("999999", ""))

    def test_parse_isoduration_for_Years(self):
        """
        Tests the parse_isoduration function for iso strings with Years
        """
        days_in_year = 365
        iso_string = "0Y"
        self.assertEqual(parse_isoduration(iso_string), 0)
        iso_string = "1Y"
        self.assertEqual(parse_isoduration(iso_string), days_in_year)
        iso_string = "4Y"
        self.assertEqual(parse_isoduration(iso_string), 4 * days_in_year)
        iso_string = "100Y"
        self.assertEqual(parse_isoduration(iso_string), 100 * days_in_year)
        iso_string = "999999Y"
        self.assertEqual(parse_isoduration(iso_string), 999999 * days_in_year)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
