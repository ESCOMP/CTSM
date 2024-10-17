#!/usr/bin/env python3

"""System tests for run_neon

"""

import glob
import os
import unittest
import tempfile
import shutil
import sys

from ctsm import unit_testing
from ctsm.site_and_regional.run_neon import main
from ctsm.path_utils import path_to_ctsm_root

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
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)  # cd to tempdir

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_one_site(self):
        """
        This test specifies a site to run
        Run the tool, check that file structure is set up correctly
        """

        # run the run_neon tool
        sys.argv = [
            os.path.join(path_to_ctsm_root(), "tools", "site_and_regional", "run_neon"),
            "--neon-sites",
            "BART",
            "--setup-only",
            "--output-root",
            self._tempdir,
        ]
        main("")

        # assert that BART directories were created during setup
        self.assertTrue("BART" in glob.glob(self._tempdir + "/BART*")[0])

        # TODO: Would also be useful to test the following items:
        # It might be good to ensure the log files are working as expected?
        # Test running transient, ad and post ad cases.
        # Test use of base case root.
        # Test for using prism?


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
