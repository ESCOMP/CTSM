#!/usr/bin/env python3

"""System tests for gen_mksurfdata_namelist

"""

import os
import re

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
        """ """
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


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
