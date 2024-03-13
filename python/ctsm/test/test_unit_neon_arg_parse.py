#!/usr/bin/env python3
"""
Unit tests for neon_arg_parse

You can run this by:
    python -m unittest test_unit_neon_arg_parse.py
"""

import unittest
import tempfile
import shutil
import os
import sys
import glob

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.neon_arg_parse import get_parser
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=invalid-name


class Test_neon_arg_parse(unittest.TestCase):
    """
    Basic class for testing neon_arg_parse.py.
    """

    def setUp(self):
        """
        Make /_tempdir for use by these tests.
        """
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_function(self):
        """
        Test that neon_arg_parse is properly reading arguments...
        """
        sys.argv = [
            "neon_arg_parse",
            "--neon-sites",
            "ABBY",
            "--experiment",
            "test",
            "--run-type",
            "ad",
        ]
        description = ""
        cesmroot = path_to_ctsm_root()
        valid_neon_sites = glob.glob(
            os.path.join(cesmroot, "cime_config", "usermods_dirs", "NEON", "[!d]*")
        )
        valid_neon_sites = sorted([v.split("/")[-1] for v in valid_neon_sites])
        parsed_arguments = get_parser(sys.argv, description, valid_neon_sites)

        self.assertEqual(parsed_arguments[0][0], "ABBY", "arguments not processed as expected")
        self.assertEqual(parsed_arguments[3], "test", "arguments not processed as expected")
        self.assertEqual(parsed_arguments[4], False, "arguments not processed as expected")
        self.assertEqual(parsed_arguments[2], "ad", "arguments not processed as expected")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
