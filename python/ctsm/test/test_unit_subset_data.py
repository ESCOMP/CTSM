#!/usr/bin/env python3
"""
Unit tests for subset_data

You can run this by:
    python -m unittest test_unit_subset_data.py
"""

import unittest
import configparser
import argparse
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.subset_data import get_parser, setup_files
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=invalid-name


class TestSubsetData(unittest.TestCase):
    """
    Basic class for testing SubsetData class in subset_data.py.
    """

    def setUp(self):
        sys.argv = ["subset_data", "point"]
        DEFAULTS_FILE = "default_data.cfg"
        parser = get_parser()
        self.args = parser.parse_args()
        self.cesmroot = path_to_ctsm_root()
        self.defaults = configparser.ConfigParser()
        self.defaults.read(os.path.join(self.cesmroot, "tools/site_and_regional", DEFAULTS_FILE))

    def test_inputdata_setup_files_basic(self):
        """
        Test 
        """
        setup_files(self.args, self.defaults, self.cesmroot)

    def test_inputdata_setup_files_inputdata_dne(self):
        """
        Test that inputdata directory does not exist
        """
        self.defaults.set("main", "clmforcingindir", "/zztop")
        with self.assertRaisesRegex(SystemExit, "inputdata directory does not exist"):
            setup_files(self.args, self.defaults, self.cesmroot)

if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
