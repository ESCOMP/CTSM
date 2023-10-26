#!/usr/bin/env python3

"""System tests for run_neon

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

from ctsm.path_utils import path_to_ctsm_root, path_to_cime
from ctsm import unit_testing
from ctsm.site_and_regional import run_neon
from ctsm.site_and_regional.run_neon import main, get_parser

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysRunNeon(unittest.TestCase):
    """System tests for run_neon"""

    def setUp(self):
        """
        Make /_tempdir for use by these tests.
        Check tempdir for history files
        """
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_one_site(self):
        """
        This test specifies a site to run
        Run the tool, check history file output exists
        Could also check log files? Although this functionally doesn't change output,
        it might be good to ensure the log files are working as expected?
        Test running transient, ad and post ad cases.
        Test use of base case root.
        Test for using prism?
        Test setup_only? This should be encapsulated within the full run?
        """

        # run the run_neon tool
        # sys.argv = ["run_neon", "--neon_sites 'ABBY'"]
        # sys.argv = ["--neon_sites", ["ABBY"]]
        sys.argv = ["--neon-sites 'ABBY'"]
        print("sys.argv:")
        print(sys.argv)
        valid_neon_sites = ["ABBY", "OSBS", "BART"]  # ["all"]
        parser = get_parser(sys.argv, "description_for_parser", valid_neon_sites)
        print("parser:")
        print(parser)
        main("")
        # this seems to run OSBS (default site, instead of ABBY),
        #  but does create files! It takes a while though, should we do setup-only?
        # Could assert that dir is created with files; we should also move this into a tempdir?
        
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # fsurdat_in_data = xr.open_dataset(self._fsurdat_in)
        # fsurdat_out_data = xr.open_dataset(self._fsurdat_out)
        # assert that fsurdat_out equals fsurdat_in
        # self.assertTrue(fsurdat_out_data.equals(fsurdat_in_data))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
