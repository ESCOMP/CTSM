#!/usr/bin/env python3

"""Unit tests for config_utils
"""

import unittest

from configparser import ConfigParser

from ctsm import unit_testing
from ctsm.config_utils import lon_range_0_to_360, get_config_value_or_array

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


class TestConfigUtils(unittest.TestCase):
    """Tests of config_utils"""

    # Allow these to be set outside of the __init__ method
    # pylint: disable=attribute-defined-outside-init
    def setUp(self):
        """Setup for testing"""
        self.config = ConfigParser()
        self.section = "main"
        self.file_path = "path_to_file"
        self.config[self.section] = {}

    def test_negative_lon(self):
        """Test lon_range_0_to_360 for a negative longitude"""
        lon = -180.0
        lon_new = lon_range_0_to_360(lon)
        self.assertEqual(lon_new, 180.0, "lon not as expected")

    def test_negative2_lon(self):
        """Test lon_range_0_to_360 for a negative longitude"""
        lon = -5.0
        lon_new = lon_range_0_to_360(lon)
        self.assertEqual(lon_new, 355.0, "lon not as expected")

    def test_regular_lon(self):
        """Test lon_range_0_to_360 for a regular longitude"""
        lon = 22.567
        lon_new = lon_range_0_to_360(lon)
        self.assertEqual(lon_new, lon, "lon not as expected")

    def test_lon_out_of_range(self):
        """Test lon_range_0_to_360 for longitude out of range"""
        lon = 361.0
        with self.assertRaisesRegex(SystemExit, "lon_in needs to be in the range 0 to 360"):
            lon_range_0_to_360(lon)

    def test_lon_out_of_range_negative(self):
        """Test lon_range_0_to_360 for longitude out of range"""
        lon = -181.0
        with self.assertRaisesRegex(SystemExit, "lon_in needs to be in the range 0 to 360"):
            lon_range_0_to_360(lon)

    def test_config_value_or_array_single_value(self):
        """Simple test of get_config_value_or_array"""
        item = "single_value_thing"
        # Test on a string, float and integer
        self.config.set(self.section, item, "one-thing")
        value = get_config_value_or_array(self.config, self.section, item)
        self.assertEqual(value, "one-thing", "Value as expected")
        self.config.set(self.section, item, "100.")
        value = get_config_value_or_array(self.config, self.section, item)
        self.assertEqual(value, "100.", "Value as expected")
        self.config.set(self.section, item, "100")
        value = get_config_value_or_array(self.config, self.section, item)
        self.assertEqual(value, "100", "Value as expected")
        # Run over again, with an explicit conversion
        self.config.set(self.section, item, "one-thing")
        value = get_config_value_or_array(self.config, self.section, item, convert_to_type=str)
        self.assertEqual(value, "one-thing", "Value as expected")
        self.config.set(self.section, item, "100.")
        value = get_config_value_or_array(self.config, self.section, item, convert_to_type=float)
        self.assertEqual(value, 100.0, "Value as expected")
        self.config.set(self.section, item, "100")
        value = get_config_value_or_array(self.config, self.section, item, convert_to_type=int)
        self.assertEqual(value, 100, "Value as expected")

    def test_config_value_or_array_for_list(self):
        """Simple test of get_config_value_or_array for a list"""
        item = "three_things"
        # Test on a string, float and integer
        mystr = "one two three"
        mystrlist = ["one", "two", "three"]
        myfloat = "1. 2. 3."
        myfloatlist = [1.0, 2.0, 3.0]
        myint = "1 2 3"
        myintlist = [1, 2, 3]
        self.config.set(self.section, item, mystr)
        value = get_config_value_or_array(self.config, self.section, item, convert_to_type=str)
        self.assertEqual(value, mystrlist, "List as expected")
        self.assertEqual(len(value), 3, "List size as expected")
        self.config.set(self.section, item, myfloat)
        value = get_config_value_or_array(self.config, self.section, item, convert_to_type=float)
        self.assertEqual(value, myfloatlist, "Value as expected")
        self.assertEqual(len(value), 3, "List size as expected")
        self.config.set(self.section, item, myint)
        value = get_config_value_or_array(self.config, self.section, item, convert_to_type=int)
        self.assertEqual(value, myintlist, "Value as expected")
        self.assertEqual(len(value), 3, "List size as expected")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
