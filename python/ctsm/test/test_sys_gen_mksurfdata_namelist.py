#!/usr/bin/env python3

"""System tests for gen_mksurfdata_namelist"""

import os

import unittest
import tempfile
import shutil
import sys

from ctsm.path_utils import path_to_ctsm_root
from ctsm.toolchain.gen_mksurfdata_namelist import main
from ctsm import unit_testing

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysGenMkSurfNML(unittest.TestCase):
    """System tests for gen_mksurfdata_namelist"""

    def setUp(self):
        """Setp temporary directory to make the files in"""
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._testinputs_path = testinputs_path
        self._original_wd = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)
        self.outfile = "surfdata.namelist"
        sys.argv = [
            "gen_mksurfdata_namelist",
            "--namelist",
            self.outfile,
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._original_wd)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_simple_namelist(self):
        """
        Test that a standard simple namelist works
        """
        sys.argv.extend(
            [
                "--start-year",
                "2000",
                "--end-year",
                "2000",
                "--res",
                "0.9x1.25",
            ]
        )
        main()
        self.assertTrue(os.path.exists(self.outfile), "Output surface dataset file should exist")

    def test_nocrop_inlandwet_glc_namelist(self):
        """
        Test a namelist with several options on
        """
        sys.argv.extend(
            [
                "--start-year",
                "1850",
                "--end-year",
                "1850",
                "--res",
                "1.9x2.5",
                "--nocrop",
                "--inlandwet",
                "--glc",
            ]
        )
        main()
        self.assertTrue(os.path.exists(self.outfile), "Output surface dataset file should exist")

    def test_hires_namelist(self):
        """
        Test that a high resolution namelist works
        """
        sys.argv.extend(
            [
                "--start-year",
                "1850",
                "--end-year",
                "1850",
                "--res",
                "mpasa15",
                "--glc-nec",
                "10",
                "--hires_soitex",
            ]
        )
        main()
        self.assertTrue(os.path.exists(self.outfile), "Output surface dataset file should exist")

    def test_ssp_transient_namelist(self):
        """
        Test that a SSP transient namelist works
        """
        sys.argv.extend(
            [
                "--start-year",
                "1850",
                "--end-year",
                "2100",
                "--res",
                "ne30np4.pg3",
                "--ssp-rcp",
                "SSP2-4.5",
                "--nosurfdata",
            ]
        )
        main()
        self.assertTrue(os.path.exists(self.outfile), "Output surface dataset file should exist")

    def test_potveg_namelist(self):
        """
        Test that a potential vegetation namelist works
        """
        sys.argv.extend(
            [
                "--start-year",
                "1850",
                "--end-year",
                "1850",
                "--res",
                "4x5",
                "--potveg_flag",
            ]
        )
        main()
        self.assertTrue(os.path.exists(self.outfile), "Output surface dataset file should exist")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
