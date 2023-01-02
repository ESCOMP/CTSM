#!/usr/bin/env python3

"""Unit tests for config_utils
"""

import unittest

from ctsm import unit_testing
from ctsm.config_utils import lon_range_0_to_360

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestConfigUtils(unittest.TestCase):
    """Tests of config_utils"""

    def test_negative_lon(self):
        """ Test lon_range_0_to_360 for a negative longitude"""
        lon = -180.
        lon_new = lon_range_0_to_360(lon)
        self.assertEqual( lon_new, 180., "lon not as expected" )

    def test_negative_lon(self):
        """ Test lon_range_0_to_360 for a negative longitude"""
        lon = -5.
        lon_new = lon_range_0_to_360(lon)
        self.assertEqual( lon_new, 355., "lon not as expected" )


    def test_regular_lon(self):
        """ Test lon_range_0_to_360 for a regular longitude"""
        lon = 22.567
        lon_new = lon_range_0_to_360(lon)
        self.assertEqual( lon_new, lon, "lon not as expected" )

    def test_lon_out_of_range(self):
        """ Test lon_range_0_to_360 for longitude out of range"""
        lon = 361.
        with self.assertRaisesRegex(SystemExit, "lon_in needs to be in the range 0 to 360"):
           lon_range_0_to_360(lon)

    def test_lon_out_of_range_negative(self):
        """ Test lon_range_0_to_360 for longitude out of range"""
        lon = -181.
        with self.assertRaisesRegex(SystemExit, "lon_in needs to be in the range 0 to 360"):
           lon_range_0_to_360(lon)
