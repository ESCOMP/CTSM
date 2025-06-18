#!/usr/bin/env python3

"""System tests for plumber2_surf_wrapper"""

import glob
import os
import unittest
import tempfile
import shutil
import sys

from ctsm import unit_testing
from ctsm.site_and_regional.plumber2_surf_wrapper import main
from ctsm.site_and_regional.plumber2_shared import read_plumber2_sites_csv
from ctsm.path_utils import path_to_ctsm_root

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysPlumber2SurfWrapper(unittest.TestCase):
    """
    System tests for plumber2_surf_wrapper
    """

    def setUp(self):
        """
        Make tempdir for use by these tests.
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

    def test_plumber2_surf_wrapper(self):
        """
        Run the entire tool
        """

        tool_path = os.path.join(
            path_to_ctsm_root(),
            "tools",
            "site_and_regional",
            "plumber2_surf_wrapper",
        )
        sys.argv = [tool_path]
        main()

        # How many files do we expect?
        plumber2_csv = read_plumber2_sites_csv()
        n_files_expected = len(plumber2_csv)

        # How many files did we get?
        file_list = os.listdir("subset_data_single_point")
        n_files = len(file_list)

        # Check
        self.assertEqual(n_files_expected, n_files)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
