#!/usr/bin/env python3
"""
Unit tests for the SystemTest sspmatrix.py
"""

import unittest
import os
import sys
import tempfile
import shutil
import logging
from pathlib import Path

from ctsm.path_utils import add_cime_lib_to_path, add_ctsm_systests_to_path
from ctsm import unit_testing

add_ctsm_systests_to_path()
add_cime_lib_to_path()

# pylint: disable=import-error
# pylint: disable=wrong-import-position
# pylint: disable=wrong-import-order
from SystemTests.sspmatrixcn import SSPMATRIXCN

# pylint: disable=import-error
# pylint: disable=wrong-import-position
# pylint: disable=wrong-import-order
from CIME.tests.case_fake import CaseFake

# pylint: disable=invalid-name

logger = logging.getLogger(__name__)


class SSPCaseFake(CaseFake):
    """
    Extend the CaseFake class with a couple things needed here
    """

    def __init__(self, case_root, tempdir, create_case_root=True):
        """
        Initialization handling the tempdir
        """
        super().__init__(case_root, create_case_root)
        self._tempdir = tempdir

    def create_clone(
        self,
        newcase,
        keepexe=False,
    ):
        """
        Extend to handle creation of user_nl_clm file
        """
        clone = super().create_clone(newcase, keepexe=keepexe)
        os.mknod(os.path.join(newcase, "user_nl_clm"))
        # Also make the needed case directories
        clone.make_case_dirs(self._tempdir)
        return clone

    def make_case_dirs(self, tempdir):
        """
        Create the directories needed for the CASE
        """
        casename = self.get_value("CASE")
        dout_s_root = Path(tempdir, "archive", casename)
        dout_s_root.mkdir(parents=True)
        self.set_value("DOUT_S_ROOT", str(dout_s_root))
        rest_root = Path(dout_s_root, "rest")
        rest_root.mkdir()
        rundir = Path(tempdir, casename, "run")
        rundir.mkdir()
        self.set_value("RUNDIR", rundir)

    def __str__(self):
        """
        String method
        """
        return "caseroot=%s" % (self.get_value("CASEROOT"))


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

        # Set up the testing SSPCaseFake
        caseroot = os.path.join(self._tempdir, "ssptestcase")
        self._case = SSPCaseFake(caseroot, tempdir=self._tempdir, create_case_root=True)
        os.chdir(caseroot)

        # Set XML variables that will be needed in the case
        self._case.set_value("DATM_YR_START", 2000)
        self._case.set_value("DATM_YR_END", 2001)
        self._case.set_value("COMP_LND", "clm")
        self._case.set_value("NINST", 1)

        self.ssp = SSPMATRIXCN(self._case)
        self._case.make_case_dirs(self._tempdir)

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

    def test_run_phase(self):
        """
        Test doing the standard run_phase, that does each step
        """
        self.ssp.run_phase()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
