#!/usr/bin/env python3

"""Unit tests for config_utils"""

import unittest
from argparse import ArgumentTypeError
import numpy as np

from ctsm import unit_testing
from ctsm.longitude import Longitude
from ctsm.longitude import _convert_lon_type_180_to_360, _convert_lon_type_360_to_180
from ctsm.longitude import _check_lon_type_180, _check_lon_type_360
from ctsm.longitude import detect_lon_type

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access,too-many-public-methods


class TestLongitude(unittest.TestCase):
    """Tests of Longitude class and helper functions"""

    # Checking longitude type

    def test_check_lon_type_180(self):
        """
        Check that a single value in [-180, 180] passes _check_lon_type_180()
        """
        _check_lon_type_180(-180)
        _check_lon_type_180(-55)
        _check_lon_type_180(55)
        _check_lon_type_180(155)
        _check_lon_type_180(180)

    def test_check_lon_type_180_list(self):
        """
        Check that a list with all values in [-180, 180] passes _check_lon_type_180()
        """
        _check_lon_type_180([-180, -55, 55, 155, 180])

    def test_check_lon_type_180_errors(self):
        """
        Check that a single value outside [-180, 180] fails _check_lon_type_180()
        """
        msg = r"lon_in must be in the range \[-180, 180\]"
        with self.assertRaisesRegex(ValueError, msg):
            _check_lon_type_180(-181)
        with self.assertRaisesRegex(ValueError, msg):
            _check_lon_type_180(181)

    def test_check_lon_type_180_list_errors(self):
        """
        Check that a list with one value outside [-180, 180] fails _check_lon_type_180()
        """
        with self.assertRaisesRegex(ValueError, r"lon_in must be in the range \[-180, 180\]"):
            _check_lon_type_180([-191, -180, -55, 55, 155, 180])

    def test_check_lon_type_360(self):
        """
        Check that a single value in [0, 360] passes _check_lon_type_360()
        """
        _check_lon_type_360(0)
        _check_lon_type_360(55)
        _check_lon_type_360(180)
        _check_lon_type_360(360)

    def test_check_lon_type_360_list(self):
        """
        Check that a list with all values in [0, 360] passes _check_lon_type_360()
        """
        _check_lon_type_360([0, 55, 180, 360])

    def test_check_lon_type_360_errors(self):
        """
        Check that a single value outside [0, 360] fails _check_lon_type_360()
        """
        msg = r"lon_in must be in the range \[0, 360\]"
        with self.assertRaisesRegex(ValueError, msg):
            _check_lon_type_360(-1)
        with self.assertRaisesRegex(ValueError, msg):
            _check_lon_type_360(361)

    def test_check_lon_type_360_list_errors(self):
        """
        Check that a list with one value outside [0, 360] fails _check_lon_type_360()
        """
        with self.assertRaisesRegex(ValueError, r"lon_in must be in the range \[0, 360\]"):
            _check_lon_type_360([-155, 0, 55, 180, 360])

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
        with self.assertRaisesRegex(ValueError, r"lon_in must be in the range \[-180, 180\]"):
            _convert_lon_type_180_to_360(lon)

    def test_convert_lon_type_180_to_360_toolow(self):
        """Test conversion 180→360 for a value < -180: Should error"""
        lon = -555

        with self.assertRaisesRegex(ValueError, r"lon_in must be in the range \[-180, 180\]"):
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
        with self.assertRaisesRegex(ValueError, r"lon_in must be in the range \[0, 360\]"):
            _convert_lon_type_360_to_180(lon)

    def test_convert_lon_type_360_to_180_toolow(self):
        """Test conversion 360→180 for a value < -180: Should error"""
        lon = -555

        with self.assertRaisesRegex(ValueError, r"lon_in must be in the range \[0, 360\]"):
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

    def test_lon_eq_both360(self):
        """Test that == works for two equal Longitudes both of type 360"""
        lon1 = Longitude(275, 360)
        lon2 = Longitude(275, 360)
        self.assertTrue(lon1 == lon2)

    def test_lon_eq_360num_error(self):
        """Test that == fails if RHS isn't Longitude"""
        lon1 = Longitude(275, 360)
        lon2 = 275
        with self.assertRaisesRegex(
            TypeError, "Comparison not supported between instances of 'Longitude' and "
        ):
            _ = lon1 == lon2

    def test_lon_eq_num360_error(self):
        """Test that == fails if LHS isn't Longitude"""
        lon1 = 275
        lon2 = Longitude(275, 360)
        with self.assertRaisesRegex(
            TypeError, "Comparison not supported between instances of 'Longitude' and "
        ):
            _ = lon1 == lon2

    def test_lon_eq_both180(self):
        """Test that == works for two equal Longitudes both of type 180"""
        lon1 = Longitude(-5, 180)
        lon2 = Longitude(-5, 180)
        self.assertTrue(lon1 == lon2)

    def test_lon_eq_180360(self):
        """Test that == works for two equal Longitudes of different types"""
        lon1 = Longitude(-5, 180)
        lon2 = Longitude(355, 360)
        self.assertTrue(lon1 == lon2)
        self.assertTrue(lon2 == lon1)

    def test_lon_eqfalse_180360(self):
        """Test that == works for two unequal Longitudes of different types"""
        lon1 = Longitude(-1, 180)
        lon2 = Longitude(355, 360)
        self.assertFalse(lon1 == lon2)
        self.assertFalse(lon2 == lon1)

    def test_lon_noteqtrue_180360(self):
        """Test that != works for two unequal Longitudes of different types"""
        lon1 = Longitude(-1, 180)
        lon2 = Longitude(355, 360)
        self.assertTrue(lon1 != lon2)
        self.assertTrue(lon2 != lon1)

    def test_lon_compare_both360(self):
        """
        Ensure that comparison operators work if both are type 360
        """
        lon1 = Longitude(155, 360)
        lon2 = Longitude(150, 360)
        self.assertTrue(lon1 > lon2)
        self.assertTrue(lon1 >= lon2)
        self.assertFalse(lon1 <= lon2)
        self.assertFalse(lon1 < lon2)

    def test_lon_compare_both180(self):
        """
        Ensure that comparison operators work if both are type 180
        """
        lon1 = Longitude(155, 180)
        lon2 = Longitude(150, 180)
        self.assertTrue(lon1 > lon2)
        self.assertTrue(lon1 >= lon2)
        self.assertFalse(lon1 <= lon2)
        self.assertFalse(lon1 < lon2)

    def test_lon_compare_arrays(self):
        """
        Ensure that comparison operators work elementwise for arrays
        """
        lon1 = Longitude(np.array([155, 156, 157]), 180)
        lon2 = Longitude(np.array([155, 157, 156]), 180)
        result = lon1 == lon2
        expected = [True, False, False]
        for i in np.arange(3):
            self.assertTrue(result[i] == expected[i])
        result = lon1 < lon2
        expected = [False, True, False]
        for i in np.arange(3):
            self.assertTrue(result[i] == expected[i])
        result = lon1 <= lon2
        expected = [True, True, False]
        for i in np.arange(3):
            self.assertTrue(result[i] == expected[i])
        result = lon1 > lon2
        expected = [False, False, True]
        for i in np.arange(3):
            self.assertTrue(result[i] == expected[i])
        result = lon1 >= lon2
        expected = [True, False, True]
        for i in np.arange(3):
            self.assertTrue(result[i] == expected[i])

    def test_lon_compare_lists(self):
        """
        Ensure that comparison operators work on the entire object for lists
        """
        lon1 = Longitude([155, 156, 157], 360)
        lon2 = Longitude([155, 157, 156], 360)
        self.assertFalse(lon1 == lon2)
        self.assertTrue(lon1 != lon2)
        lon2 = lon1
        self.assertTrue(lon1 == lon2)
        self.assertFalse(lon1 != lon2)

    def test_lon_compare_difftypes_error(self):
        """
        Ensure that comparison operators fail if Longitudes are different types
        """
        lon1 = Longitude(155, 360)
        lon2 = Longitude(150, 180)
        msg = "Comparison not supported between Longitudes of different types"
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 < lon2
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 > lon2
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 <= lon2
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 >= lon2

    def test_lon_compare_notlon_error(self):
        """
        Ensure that comparison operators fail if one isn't a Longitude
        """
        lon1 = Longitude(155, 360)
        lon2 = 255
        msg = "Comparison not supported between instances of 'Longitude' and"
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 < lon2
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 > lon2
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 <= lon2
        with self.assertRaisesRegex(TypeError, msg):
            _ = lon1 >= lon2

    def test_detect_lon_type_mid_180(self):
        """test that detect_lon_type works for an unambiguously 180 value"""
        self.assertEqual(detect_lon_type(-150), 180)

    def test_detect_lon_type_min_180(self):
        """test that detect_lon_type works at -180"""
        self.assertEqual(detect_lon_type(-180), 180)

    def test_detect_lon_type_mid_360(self):
        """test that detect_lon_type works for an unambiguously 360 value"""
        self.assertEqual(detect_lon_type(355), 360)

    def test_detect_lon_type_max_360(self):
        """test that detect_lon_type works at 360"""
        self.assertEqual(detect_lon_type(360), 360)

    def test_detect_lon_type_list_180(self):
        """test that detect_lon_type works for a list with just one unambiguously 180 value"""
        self.assertEqual(detect_lon_type([-150, 150]), 180)

    def test_detect_lon_type_list_360(self):
        """test that detect_lon_type works for a list with just one unambiguously 360 value"""
        self.assertEqual(detect_lon_type([256, 150]), 360)

    def test_detect_lon_type_ambig(self):
        """test that detect_lon_type fails if ambiguous"""
        with self.assertRaisesRegex(ArgumentTypeError, r"Longitude\(s\) ambiguous"):
            detect_lon_type(150)

    def test_detect_lon_type_list_ambig(self):
        """test that detect_lon_type fails for an ambiguous list"""
        with self.assertRaisesRegex(ArgumentTypeError, r"Longitude\(s\) ambiguous"):
            detect_lon_type([150, 170])

    def test_detect_lon_type_list_both(self):
        """test that detect_lon_type fails for a list with unambiguous members of both types"""
        with self.assertRaisesRegex(RuntimeError, r"Longitude array contains values of both types"):
            detect_lon_type([-150, 270])

    def test_detect_lon_type_ambig0(self):
        """test that detect_lon_type fails at 0"""
        with self.assertRaisesRegex(ArgumentTypeError, r"Longitude\(s\) ambiguous"):
            detect_lon_type(0)

    def test_detect_lon_type_oob_low(self):
        """test that detect_lon_type fails if out of bounds below min"""
        with self.assertRaisesRegex(ValueError, r"\(Minimum\) longitude < -180"):
            detect_lon_type(-300)

    def test_detect_lon_type_oob_high(self):
        """test that detect_lon_type fails if out of bounds above max"""
        with self.assertRaisesRegex(ValueError, r"\(Maximum\) longitude > 360"):
            detect_lon_type(500)

    def test_list_as_lon(self):
        """
        Test that converting a List to Longitude works
        """
        expected_lon = [1, 5]
        expected_type = 360
        result = Longitude(expected_lon, expected_type)
        self.assertEqual(expected_lon, result._lon)
        self.assertEqual(expected_type, result.lon_type())

    def test_array_as_lon(self):
        """
        Test that converting a Numpy array to Longitude works
        """
        expected_lon = np.array([1, 5])
        expected_type = 180
        result = Longitude(expected_lon, expected_type)
        self.assertTrue(np.array_equal(expected_lon, result._lon))
        self.assertEqual(expected_type, result.lon_type())

    def test_lon_type_getter(self):
        """Test lon_type() getter method"""
        lon = Longitude(55, 180)
        self.assertEqual(lon.lon_type(), 180)

    def test_no_implicit_string_conversion(self):
        """Ensure that implicit string conversion is disallowed"""
        lon = Longitude(55, 180)
        with self.assertRaisesRegex(
            NotImplementedError, r"Use Longitude\.get_str\(\) instead of implicit string conversion"
        ):
            _ = f"{lon}"
        with self.assertRaisesRegex(
            NotImplementedError, r"Use Longitude\.get_str\(\) instead of implicit string conversion"
        ):
            _ = str(lon)

    def test_get_str(self):
        """Ensure that explicit string conversion works as expected"""
        lon = Longitude(55, 180)
        self.assertEqual(lon.get_str(180), "55.0")
        self.assertEqual(lon.get_str(360), "55.0")
        lon = Longitude(-55, 180)
        self.assertEqual(lon.get_str(180), "-55.0")
        self.assertEqual(lon.get_str(360), "305.0")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
