#!/usr/bin/env python3
"""
Unit tests for arg_utils.py function and types.

You can run this by:
    python -m unittest test_unit_args_utils.py
"""

import os
import sys
import unittest
import argparse

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm.args_utils import plat_type
from ctsm import unit_testing

# pylint: disable=invalid-name


class TestArgsPlat(unittest.TestCase):
    """
    Tests for plat_type in args_util.py
    """

    def test_platType_outOfBounds_positive(self):
        """
        Test of plat_type bigger than 90
        """
        with self.assertRaisesRegex(argparse.ArgumentTypeError, "Latitude.*should be between"):
            _ = plat_type(91)

    def test_platType_outOfBounds_pos90(self):
        """
        Test of plat_type is 90
        """
        result = plat_type(90)
        self.assertEqual(result, 90.0)

    def test_platType_outOfBounds_neg90(self):
        """
        Test of plat_type is -90
        """
        result = plat_type(-90)
        self.assertEqual(result, -90.0)

    def test_platType_outOfBounds_negative(self):
        """
        Test of plat_type smaller than -90
        """
        with self.assertRaisesRegex(argparse.ArgumentTypeError, "Latitude.*should be between"):
            _ = plat_type(-91)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
