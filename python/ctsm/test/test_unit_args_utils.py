#!/usr/bin/env python3

import os
import sys
import unittest
import argparse

# -- add python/ctsm  to path
_CTSM_PYTHON = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
print (_CTSM_PYTHON)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm.args_utils import plon_type, plat_type
from ctsm import unit_testing

# pylint: disable=invalid-name

class TestArgsPlon(unittest.TestCase):

    # --between 0-360
    def test_plonType_positive(self):
        result = plon_type(30)
        self.assertEqual(result, 30.0)

    # --between -180-0
    def test_plonType_negative(self):
        result = plon_type(-30)
        self.assertEqual(result, 330.0)

    # -- > 360
    def test_plonType_outOfBounds_positive(self):
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "Longitude.*should be between"
        ):
            _ = plon_type(360.5)

    # -- < -180
    def test_plonType_outOfBounds_negative(self):
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "Longitude.*should be between"
        ):
            _ = plon_type(-200)

class TestArgsPlat(unittest.TestCase):
    def test_platType_outOfBounds_positive(self):
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "Latitude.*should be between"
        ):
            _ = plat_type(91)

    def test_platType_outOfBounds_negative(self):
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "Latitude.*should be between"
        ):
            _ = plat_type(-91)

if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
