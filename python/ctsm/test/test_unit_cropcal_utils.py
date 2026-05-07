#!/usr/bin/env python3

"""Unit tests for cropcal_utils.py"""

import unittest

from ctsm import unit_testing
from ctsm.crop_calendars import cropcal_utils as ccu

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestCropCalUtils(unittest.TestCase):
    """Tests of cropcal_utils.py"""

    def setUp(self):
        self.vegtype_mainlist = ["crop_1", "crop_2", "crop_3"]

    def test_vegtype_str2int_1string(self):
        """
        Tests vegtype_str2int() for a single string. Result should be an int.
        """
        result = ccu.vegtype_str2int("crop_1", vegtype_mainlist=self.vegtype_mainlist)
        self.assertEqual(result, 0)

    def test_vegtype_str2int_2strings(self):
        """
        Tests vegtype_str2int() for two strings. result should be a list of ints.
        """
        result = ccu.vegtype_str2int(["crop_1", "crop_3"], vegtype_mainlist=self.vegtype_mainlist)
        self.assertListEqual(result, [0, 2])


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
