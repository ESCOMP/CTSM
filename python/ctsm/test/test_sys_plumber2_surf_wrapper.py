#!/usr/bin/env python3

"""System tests for plumber2_surf_wrapper"""

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

        # Path to script
        self.tool_path = os.path.join(
            path_to_ctsm_root(),
            "tools",
            "site_and_regional",
            "plumber2_surf_wrapper",
        )

        # Path to test inputs directory
        self.test_inputs = os.path.join(
            os.path.dirname(__file__), "testinputs", "plumber2_surf_wrapper"
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_plumber2_surf_wrapper(self):
        """
        Run the entire tool with default settings.
        CAN ONLY RUN ON SYSTEMS WITH INPUTDATA
        """

        sys.argv = [self.tool_path]
        main()

        # How many files do we expect?
        plumber2_csv = read_plumber2_sites_csv()
        n_files_expected = len(plumber2_csv)

        # How many files did we get?
        file_list = os.listdir("subset_data_single_point")
        n_files = len(file_list)

        # Check
        self.assertEqual(n_files_expected, n_files)

    def test_plumber2_surf_wrapper_16pft(self):
        """
        Run the entire tool with --16pft.
        CAN ONLY RUN ON SYSTEMS WITH INPUTDATA
        """

        sys.argv = [self.tool_path, "--16pft"]
        main()

        # How many files do we expect?
        plumber2_csv = read_plumber2_sites_csv()
        n_files_expected = len(plumber2_csv)

        # How many files did we get?
        file_list = os.listdir("subset_data_single_point")
        n_files = len(file_list)

        # Check
        self.assertEqual(n_files_expected, n_files)

    def test_plumber2_surf_wrapper_invalid_pft(self):
        """
        plumber2_surf_wrapper should error if invalid PFT is given
        """

        sys.argv = [
            self.tool_path,
            "--plumber2-sites-csv",
            os.path.join(self.test_inputs, "PLUMBER2_sites_invalid_pft.csv"),
        ]
        with self.assertRaisesRegex(RuntimeError, "must be a valid PFT"):
            main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
