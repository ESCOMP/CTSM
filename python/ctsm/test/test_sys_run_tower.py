#!/usr/bin/env python3

"""System tests for run_tower"""

import glob
import os
import unittest
import tempfile
import shutil
import sys
import pathlib

from ctsm import unit_testing
from ctsm.site_and_regional.run_tower import main
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=import-error
# pylint: disable=wrong-import-order
from CIME.case import Case

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysRunTower(unittest.TestCase):
    """
    System tests for run_tower

    TODO: Would also be useful to test the following items:
    Ensure the log files are working as expected?
    Test use of base case root.
    Test for using prism
    """

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
        This test specifies a site to run a default case with experiment label 'TEST'
        Run the tool, check that file structure is set up correctly
        """

        # run the run_tower tool
        print("about to run tower tool")
        sys.argv = [
            os.path.join(path_to_ctsm_root(), "tools", "site_and_regional", "run_tower"),
            "--neon-sites",
            "BART",
            "--no-input-data-check",
            "--experiment",
            "TEST",
            "--output-root",
            self._tempdir,
        ]
        main("")

        # assert that BART directories were created during setup
        self.assertTrue("BART" in glob.glob(self._tempdir + "/BART*")[0])
        print(glob.glob(self._tempdir))

    def test_ad_site(self):
        """
        This test specifies a site to run an 'ad' case for
        Run the tool, check that file structure is set up correctly
        """

        # run the run_tower tool
        sys.argv = [
            os.path.join(path_to_ctsm_root(), "tools", "site_and_regional", "run_tower"),
            "--neon-sites",
            "ABBY",
            "--no-input-data-check",
            "--run-type",
            "ad",
            "--output-root",
            self._tempdir,
        ]
        main("")

        # assert that ABBY directories were created during setup
        self.assertTrue("ABBY" in glob.glob(self._tempdir + "/ABBY*")[0])

    def test_plumber_site(self):
        """
        This test specifies a site to run a default plumber AD case with
        experiment label 'TEST'. Run the tool, check that file structure is set up.
        """

        # run the run_tower tool for plumber site
        sys.argv = [
            os.path.join(path_to_ctsm_root(), "tools", "site_and_regional", "run_tower"),
            "--plumber-sites",
            "AR-SLu",
            "--no-input-data-check",
            "--experiment",
            "TEST",
            "--output-root",
            self._tempdir,
        ]
        main("")

        # assert that AR-SLu directories were created during setup
        self.assertTrue("AR-SLu" in glob.glob(self._tempdir + "/AR-SLu*")[0])

    def test_xmlchange(self):
        """
        This test checks that the --xmlchange argument is obeyed.
        """

        # run the run_tower tool
        sys.argv = [
            os.path.join(path_to_ctsm_root(), "tools", "site_and_regional", "run_tower"),
            "--neon-sites",
            "BART",
            "--no-input-data-check",
            "--experiment",
            "TEST",
            "--xmlchange",
            "CLM_CO2_TYPE=constant,CCSM_CO2_PPMV=1987",
            "--output-root",
            self._tempdir,
        ]
        main("")

        # Check that the --xmlchange argument is obeyed
        case_dir = os.path.join(self._tempdir, "BART.TEST.transient")
        self.assertTrue(os.path.exists(case_dir))
        with Case(case_dir, read_only=True) as case:
            value = case.get_value("CLM_CO2_TYPE")
            print(f"CLM_CO2_TYPE = {value}")
            self.assertTrue(value == "constant")
            value = int(case.get_value("CCSM_CO2_PPMV"))
            print(f"CCSM_CO2_PPMV = {value}")
            self.assertTrue(int(value) == 1987)

    def test_setup_only(self):
        """
        This test checks that the --setup-only argument is obeyed
        """

        # run the run_tower tool
        site_name = "BART"
        sys.argv = [
            os.path.join(path_to_ctsm_root(), "tools", "site_and_regional", "run_tower"),
            "--neon-sites",
            site_name,
            "--setup-only",
            "--experiment",
            "TEST",
            "--output-root",
            self._tempdir,
        ]
        main("")

        # make sure that build didn't happen: this dir should be empty
        case_dir = os.path.join(self._tempdir, site_name)
        build_dir_to_check = os.path.join(case_dir, "bld", "cpl", "obj")
        self.assertTrue(os.path.exists(build_dir_to_check))
        self.assertTrue(len(os.listdir(build_dir_to_check)) == 0)

    def test_overwrite(self):
        """
        This test checks that the --overwrite argument is obeyed
        """

        # run the run_tower tool once
        site_name = "BART"
        sys.argv = [
            os.path.join(path_to_ctsm_root(), "tools", "site_and_regional", "run_tower"),
            "--neon-sites",
            site_name,
            "--no-input-data-check",
            "--experiment",
            "TEST",
            "--output-root",
            self._tempdir,
        ]
        main("")

        # create a file that should be erased during the upcoming overwrite
        case_dir = os.path.join(self._tempdir, site_name)
        test_file = os.path.join(case_dir, "test_file")
        pathlib.Path(test_file).touch()
        self.assertTrue(os.path.exists(test_file))

        # run the tool again, overwriting existing
        sys.argv += ["--overwrite"]
        main("")

        # ensure that file we created is gone
        self.assertFalse(os.path.exists(test_file))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
