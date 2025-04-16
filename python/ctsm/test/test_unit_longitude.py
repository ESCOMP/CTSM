#!/usr/bin/env python3

"""Unit tests for config_utils"""

import unittest

from ctsm import unit_testing
from ctsm.longitude import Longitude
from ctsm.longitude import _convert_lon_type_180_to_360, _convert_lon_type_360_to_180

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


class TestLongitude(unittest.TestCase):
    """Tests of Longitude class and helper functions"""

    # Converting between types 180 and 360

    def test_convert_lon_type_180_to_360_positive(self):
        """Test conversion 180→360 for a middle positive longitude"""
        lon = 80
        lon_new = _convert_lon_type_180_to_360(lon)
        self.assertEqual(lon_new, 80)

    def test_convert_lon_type_180_to_360_negative(self):
        """Test conversion 180→360 for a middle negative longitude"""
        lon = -80
        lon_new = _convert_lon_type_180_to_360(lon)
        self.assertEqual(lon_new, 280)

    def test_convert_lon_type_180_to_360_lowerbound(self):
        """Test conversion 180→360 at the lower bound of [-180, 180]"""
        lon = -180.0
        lon_new = _convert_lon_type_180_to_360(lon)
        self.assertEqual(lon_new, 180)

    def test_convert_lon_type_180_to_360_upperbound(self):
        """Test conversion 180→360 at the upper bound of [-180, 180]"""
        lon = 180.0
        lon_new = _convert_lon_type_180_to_360(lon)
        self.assertEqual(lon_new, 180)

    def test_convert_lon_type_180_to_360_toohigh(self):
        """Test conversion 180→360 for a value > 180: Should error"""
        lon = 555
        with self.assertRaisesRegex(ValueError, r"lon_in needs to be in the range \[-180, 180\]"):
            _convert_lon_type_180_to_360(lon)

    def test_convert_lon_type_180_to_360_toolow(self):
        """Test conversion 180→360 for a value < -180: Should error"""
        lon = -555

        with self.assertRaisesRegex(ValueError, r"lon_in needs to be in the range \[-180, 180\]"):
            _convert_lon_type_180_to_360(lon)

    def test_convert_lon_type_360_to_180_positive_low(self):
        """Test conversion 360→180 for a low positive longitude"""
        lon = 80
        lon_new = _convert_lon_type_360_to_180(lon)
        self.assertEqual(lon_new, 80)

    def test_convert_lon_type_360_to_180_positive_high(self):
        """Test conversion 360→180 for a high positive longitude"""
        lon = 270
        lon_new = _convert_lon_type_360_to_180(lon)
        self.assertEqual(lon_new, -90)

    def test_convert_lon_type_360_to_180_lowerbound(self):
        """Test conversion 360→180 at the lower bound of [0, 360]"""
        lon = 0
        lon_new = _convert_lon_type_360_to_180(lon)
        self.assertEqual(lon_new, 0)

    def test_convert_lon_type_360_to_180_upperbound(self):
        """Test conversion 360→180 at the upper bound of [0, 360]"""
        lon = 360
        lon_new = _convert_lon_type_360_to_180(lon)
        self.assertEqual(lon_new, 0)

    def test_convert_lon_type_360_to_180_toohigh(self):
        """Test conversion 360→180 for a value > 180: Should error"""
        lon = 555
        with self.assertRaisesRegex(ValueError, r"lon_in needs to be in the range \[0, 360\]"):
            _convert_lon_type_360_to_180(lon)

    def test_convert_lon_type_360_to_180_toolow(self):
        """Test conversion 360→180 for a value < -180: Should error"""
        lon = -555

        with self.assertRaisesRegex(ValueError, r"lon_in needs to be in the range \[0, 360\]"):
            _convert_lon_type_360_to_180(lon)

    # Initializing new Longitude objects

    def test_lon_obj_lon_errors(self):
        """Trying to initialize a Longitude object with a Longitude object should error"""
        lon_type = 360
        lon_obj = Longitude(0, lon_type)
        with self.assertRaises(TypeError):
            Longitude(lon_obj, lon_type)

    def test_lon_obj_type180_neg(self):
        """Test that creating an in-bounds negative Longitude of type 180 works"""
        lon_type = 180
        this_lon = -55
        lon_obj = Longitude(this_lon, lon_type)
        self.assertEqual(lon_obj.get(lon_type), this_lon)

    def test_lon_obj_type180_pos(self):
        """Test that creating an in-bounds positive Longitude of type 180 works"""
        lon_type = 180
        this_lon = 87
        lon_obj = Longitude(this_lon, lon_type)
        self.assertEqual(lon_obj.get(lon_type), this_lon)

    def test_lon_obj_type180_min(self):
        """Test that creating a lower-bound Longitude of type 180 works"""
        lon_type = 180
        this_lon = -180
        lon_obj = Longitude(this_lon, lon_type)
        self.assertEqual(lon_obj.get(lon_type), this_lon)

    def test_lon_obj_type180_max(self):
        """Test that creating an upper-bound Longitude of type 180 works"""
        lon_type = 180
        this_lon = 180
        lon_obj = Longitude(this_lon, lon_type)
        self.assertEqual(lon_obj.get(lon_type), this_lon)

    def test_lon_obj_type360_pos(self):
        """Test that creating an in-bounds Longitude of type 360 works"""
        lon_type = 360
        this_lon = 87
        lon_obj = Longitude(this_lon, lon_type)
        self.assertEqual(lon_obj.get(lon_type), this_lon)

    def test_lon_obj_type360_min(self):
        """Test that creating a lower-bound Longitude of type 360 works"""
        lon_type = 360
        this_lon = 0
        lon_obj = Longitude(this_lon, lon_type)
        self.assertEqual(lon_obj.get(lon_type), this_lon)

    def test_lon_obj_type360_max(self):
        """Test that creating an upper-bound Longitude of type 360 works"""
        lon_type = 360
        this_lon = 360
        lon_obj = Longitude(this_lon, lon_type)
        self.assertEqual(lon_obj.get(lon_type), this_lon)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
