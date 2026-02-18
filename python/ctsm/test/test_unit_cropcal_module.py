#!/usr/bin/env python3

"""Unit tests for cropcal_module.py"""

import unittest

import numpy as np

from ctsm import unit_testing
from ctsm.crop_calendars import cropcal_module as ccu

# Allow names that pylint doesn't like, because otherwise it's hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestCushionGsLength(unittest.TestCase):
    """Tests of cushion_gs_length()"""

    def setUp(self):
        self._mxmat_dict = {
            "cropA": 55,
            "cropB": 57,
            "cropC": 87,
            "cropD": 86,
        }

    def test_cushion_gs_length_1(self):
        """
        Test cushion_gs_length() for a cushion of 1.
        """
        cushion = 1
        result = ccu.cushion_gs_length(self._mxmat_dict, cushion)
        for crop, mxmat_orig in self._mxmat_dict.items():
            self.assertEqual(result[crop], mxmat_orig - cushion)

    def test_cushion_gs_length_neg1(self):
        """
        Test cushion_gs_length() for a cushion of -1.
        """
        cushion = -1
        result = ccu.cushion_gs_length(self._mxmat_dict, cushion)
        for crop, mxmat_orig in self._mxmat_dict.items():
            self.assertEqual(result[crop], mxmat_orig - cushion)

    def test_cushion_gs_length_inf(self):
        """
        As test_cushion_gs_length_1 but with an infinite value in original.
        """
        self._mxmat_dict["cropB"] = np.inf
        cushion = 1
        result = ccu.cushion_gs_length(self._mxmat_dict, cushion)
        for crop, mxmat_orig in self._mxmat_dict.items():
            if np.isinf(mxmat_orig):
                self.assertTrue(np.isinf(result[crop]))
            else:
                self.assertEqual(result[crop], mxmat_orig - cushion)

    def test_cushion_gs_length_limit_range_default_min1(self):
        """
        Test that cushion_gs_length() limits output range to a minimum of 1 by default
        """
        cushion = np.inf

        with self.assertWarnsRegex(RuntimeWarning, "increasing that to"):
            result = ccu.cushion_gs_length(self._mxmat_dict, cushion)

        self.assertTrue(all(x == 1 for x in result.values()))

    def test_cushion_gs_length_limit_range_default_max365(self):
        """
        Test that cushion_gs_length() limits output range to a maximum of 365 by default
        """
        cushion = -np.inf

        with self.assertWarnsRegex(RuntimeWarning, "decreasing that to"):
            result = ccu.cushion_gs_length(self._mxmat_dict, cushion)

        print(result)
        self.assertTrue(all(x == 365 for x in result.values()))

    def test_cushion_gs_length_limit_range_custom_min(self):
        """
        Test cushion_gs_length() with custom min
        """
        cushion = np.inf
        min_mxmat = -1

        with self.assertWarnsRegex(RuntimeWarning, "increasing that to"):
            result = ccu.cushion_gs_length(self._mxmat_dict, cushion, min_mxmat=min_mxmat)

        self.assertTrue(all(x == min_mxmat for x in result.values()))

    def test_cushion_gs_length_limit_range_custom_max(self):
        """
        Test cushion_gs_length() with custom max
        """
        cushion = -np.inf
        max_mxmat = 400

        with self.assertWarnsRegex(RuntimeWarning, "decreasing that to"):
            result = ccu.cushion_gs_length(self._mxmat_dict, cushion, max_mxmat=max_mxmat)

        self.assertTrue(all(x == max_mxmat for x in result.values()))

    def test_cushion_gs_length_assert_orig_too_low(self):
        """
        Test that cushion_gs_length() errors if original value is too low
        """
        self._mxmat_dict["cropC"] = -999
        with self.assertRaisesRegex(AssertionError, "is < min_mxmat"):
            ccu.cushion_gs_length(self._mxmat_dict, 1)

    def test_cushion_gs_length_assert_orig_too_high(self):
        """
        Test that cushion_gs_length() errors if original value is too high
        """
        self._mxmat_dict["cropC"] = 999
        with self.assertRaisesRegex(AssertionError, "is > max_mxmat"):
            ccu.cushion_gs_length(self._mxmat_dict, 1)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
