#!/usr/bin/env python3

"""System tests for run_neon

"""

import glob
import os
import unittest
import tempfile
import shutil
import sys

# THESE LINES ARE JUST HERE FOR TESTING
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import unit_testing
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
        Run the tool, check that file structure is set up correctly
        """

        # run the run_neon tool
        sys.argv = [
            "run_neon",
            "--neon-sites",
            "BART",
            "--setup-only",
            "--output-root",
            self._tempdir,
        ]
        valid_neon_sites = ["ABBY", "OSBS", "BART"]  # ["all"]
        parser = get_parser(sys.argv, "description_for_parser", valid_neon_sites)
        main("")

        # assert that BART directories were created during setup
        self.assertTrue("BART" in glob.glob(self._tempdir + "/*"))

        # TODO: Would also be useful to test the following items:
        # It might be good to ensure the log files are working as expected?
        # Test running transient, ad and post ad cases.
        # Test use of base case root.
        # Test for using prism?


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
