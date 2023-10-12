#!/usr/bin/env python3

"""System tests for gen_mksurfdata_namelist

"""

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
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_simple_namelist(self):
        """
        Test that a standard simple namelist works
        """
        # pylint: disable=no-self-use
        sys.argv = [
            "gen_mksurfdata_namelist",
            "--start-year",
            "2000",
            "--end-year",
            "2000",
            "--res",
            "0.9x1.25",
        ]
        main()

    def test_vic_nocrop_inlandwet_glc_namelist(self):
        """
        Test that a namelist with several options on
        """
        # pylint: disable=no-self-use
        sys.argv = [
            "gen_mksurfdata_namelist",
            "--start-year",
            "1850",
            "--end-year",
            "1850",
            "--res",
            "1.9x2.5",
            "--vic",
            "--nocrop",
            "--inlandwet",
            "--glc",
        ]
        main()

    def test_hires_namelist(self):
        """
        Test that a high resolution namelist works
        """
        # pylint: disable=no-self-use
        sys.argv = [
            "gen_mksurfdata_namelist",
            "--start-year",
            "1850",
            "--end-year",
            "1850",
            "--res",
            "mpasa15",
            "--glc-nec",
            "10",
            "--hires_pft",
            "--hires_soitex",
        ]
        main()

    def test_ssp_transient_namelist(self):
        """
        Test that a SSP transient namelist works
        """
        # pylint: disable=no-self-use
        sys.argv = [
            "gen_mksurfdata_namelist",
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
        main()

    def test_potveg_namelist(self):
        """
        Test that a potential vegetation namelist works
        """
        # pylint: disable=no-self-use
        sys.argv = [
            "gen_mksurfdata_namelist",
            "--start-year",
            "1850",
            "--end-year",
            "1850",
            "--res",
            "4x5",
            "--potveg_flag",
        ]
        main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
