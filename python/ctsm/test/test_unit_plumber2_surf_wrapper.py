#!/usr/bin/env python3
"""
Unit tests for plumber2_surf_wrapper

You can run this by:
    python -m unittest test_unit_plumber2_surf_wrapper.py
"""

import unittest
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.plumber2_surf_wrapper import get_parser

# pylint: disable=invalid-name


class TestPlumber2SurfWrapper(unittest.TestCase):
    """
    Basic class for testing plumber2_surf_wrapper.py.
    """

    def test_parser(self):
        """
        Test that parser has same defaults as expected
        """

        self.assertEqual(get_parser().argument_default, None, "Parser not working as expected")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
