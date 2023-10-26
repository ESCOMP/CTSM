#!/usr/bin/env python3

"""
System tests for modify_singlept_site_neon.py
"""

import os
import re

import unittest
import tempfile
import shutil
import sys

import xarray as xr
import numpy as np

# THESE LINES ARE JUST HERE FOR TESTING
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.site_and_regional import modify_singlept_site_neon
from ctsm.site_and_regional.modify_singlept_site_neon import main, get_parser

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysModifySingleptSiteNeon(unittest.TestCase):
    """System tests for modify_singlept_site_neon"""

    def setUp(self):
        """
        Make /_tempdir for use by these tests.
        Check tempdir for history files
        """
        self._tempdir = tempfile.mkdtemp()
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._cfg_file_path = os.path.join(
            testinputs_path, "modify_singlept_site_neon_opt_sections.cfg"
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_modify_site(self):
        """
        Test modifying a singple point site
        The primary items to test here are the following:
        1) Fields are overwritten with site-specific data for neon sites
        2) Downloaded data is used in surface dataset
        3) Check specific fields listed in update_metadata for correct output
        4) Check that a netcdf with correct formatting is created
        """
        sys.argv = ["--neon_site", ["ABBY"]]  # self._cfg_file_path] #, "ABBY"]
        print("sys.argv:")
        print(sys.argv)
        parser = get_parser()
        print("parser:")
        print(parser)
        main()  # running into error here because main doesn't take in any arguments and doesn't know what arguments are needed.


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
