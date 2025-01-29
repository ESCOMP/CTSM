#!/usr/bin/env python3
"""
Unit tests for the SystemTest sspmatrix.pty

You can run this by:
    python -m unittest test_unit_sspmatrix.py
"""

import unittest
import os
import sys
import tempfile
import shutil
import logging

from ctsm.path_utils import add_cime_lib_to_path
from ctsm import unit_testing

_CIME_PATH = add_cime_lib_to_path(standalone_only=True)

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
_CTSM_SYSTEST_DIR = os.path.join(_CIME_PATH, os.pardir, "cime_config")
sys.path.insert(1, _CTSM_PYTHON)
sys.path.insert(1, _CTSM_SYSTEST_DIR)
# pylint: disable=wrong-import-order
# pylint: disable=wrong-import-position
from SystemTests.sspmatrixcn import SSPMATRIXCN

# pylint: disable=wrong-import-order
# pylint: disable=wrong-import-position
from CIME.tests.case_fake import CaseFake

# pylint: disable=invalid-name

logger = logging.getLogger(__name__)


class TestSSPMatrix(unittest.TestCase):
    """
    Basic class for testing the sspmatrix.py SystemTest
    """

    def setUp(self):
        """
        Setup test directory
        """
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()

        # Set up the testing CaseFake
        caseroot = os.path.join(self._tempdir, "case")
        self._case = CaseFake(caseroot, create_case_root=True)
        os.chdir(caseroot)

        # Set XML variables that will be needed in the case
        self._case.set_value("DATM_YR_START", 2000)
        self._case.set_value("DATM_YR_END", 2001)
        self._case.set_value("COMP_LND", "clm")
        self._case.set_value("NINST", "1")

        self.ssp = SSPMATRIXCN(self._case)

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_logger(self):
        """
        Test the logger
        """
        stream_handler = logging.StreamHandler(sys.stdout)
        logger.addHandler(stream_handler)
        logger.level = logging.DEBUG
        logger.info("nyr_forcing = %s", self.ssp.nyr_forcing)
        for n in range(self.ssp.n_steps()):
            self.ssp.__logger__(n)
            if self.ssp.spin[n] == "sasu":
                logger.info("  SASU spinup is .true.")
                if self.ssp.sasu[n] != -999:
                    logger.info("  nyr_sasu = %s", self.ssp.sasu[n])
                if self.ssp.iloop[n] != -999:
                    logger.info("  iloop_avg = %s", self.ssp.iloop[n])

        logger.info("Total number of years %s", self.ssp.total_years())
        logger.removeHandler(stream_handler)

    def test_n_steps(self):
        """
        Test that n_steps is as expected
        """
        self.assertTrue(self.ssp.n_steps() == 3)

    def test_valid_n(self):
        """
        Test that check of n-step is good for the range it runs in
        """
        for n in range(self.ssp.n_steps()):
            self.ssp.check_n(n)

    def test_negative_n(self):
        """
        Test that fails when n-step is negative
        """
        self.assertRaises(SystemExit, self.ssp.check_n, -1)

    def test_n_too_big(self):
        """
        Test that fails when n-step is too big
        """
        self.assertRaises(SystemExit, self.ssp.check_n, self.ssp.n_steps())

    def test_append_user_nl_step2(self):
        """
        Test appending to user_nl_clm file for step 2
        """
        ufile = "user_nl_clm"
        if os.path.exists(ufile):
            os.remove(ufile)

        os.mknod(ufile)

        expect = "\nhist_nhtfrq = -8760, hist_mfilt = 2\n"
        self.ssp.append_user_nl(caseroot=".", n=2)
        log = open(ufile, "r").read()
        self.assertEqual(expect, log, "Append user_nl_clm file NOT as expected for step 2")
        os.remove(ufile)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
