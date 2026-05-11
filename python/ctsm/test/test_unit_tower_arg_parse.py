#!/usr/bin/env python3
"""
Unit tests for tower_arg_parse

You can run this by:
    python -m unittest test_unit_tower_arg_parse.py
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
from ctsm.site_and_regional.tower_arg_parse import get_parser
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=invalid-name


class Test_tower_arg_parse(unittest.TestCase):
    """
    Basic class for testing tower_arg_parse.py.
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

    def test_neon(self):
        """
        Test that tower_arg_parse is properly reading arguments for neon...
        """
        sys.argv = [
            "tower_arg_parse",
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
            os.path.join(cesmroot, "cime_config", "usermods_dirs", "clm", "NEON", "[!Fd]*")
        )
        valid_neon_sites = sorted([v.split("/")[-1] for v in valid_neon_sites])

        valid_plumber_sites = glob.glob(
            os.path.join(cesmroot, "cime_config", "usermods_dirs", "clm", "PLUMBER2", "[!d]*")
        )
        valid_plumber_sites = sorted([v.split("/")[-1] for v in valid_plumber_sites])

        parsed_arguments = get_parser(sys.argv, description, valid_neon_sites, valid_plumber_sites)

        self.assertEqual(parsed_arguments[0][0], "ABBY", "arguments not processed as expected")
        self.assertEqual(parsed_arguments[4], "test", "arguments not processed as expected")
        self.assertEqual(parsed_arguments[5], False, "arguments not processed as expected")
        self.assertEqual(parsed_arguments[3], "ad", "arguments not processed as expected")

    def test_plumber(self):
        """
        Test that tower_arg_parse is properly reading arguments for plumber...
        """
        sys.argv = [
            "tower_arg_parse",
            "--plumber-sites",
            "AU-Emr",
        ]
        description = ""
        cesmroot = path_to_ctsm_root()

        valid_neon_sites = glob.glob(
            os.path.join(cesmroot, "cime_config", "usermods_dirs", "clm", "NEON", "[!Fd]*")
        )
        valid_neon_sites = sorted([v.split("/")[-1] for v in valid_neon_sites])

        valid_plumber_sites = glob.glob(
            os.path.join(cesmroot, "cime_config", "usermods_dirs", "clm", "PLUMBER2", "[!d]*")
        )
        valid_plumber_sites = sorted([v.split("/")[-1] for v in valid_plumber_sites])

        parsed_arguments = get_parser(sys.argv, description, valid_neon_sites, valid_plumber_sites)
        self.assertEqual(parsed_arguments[1][0], "AU-Emr", "arguments not processed as expected")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
