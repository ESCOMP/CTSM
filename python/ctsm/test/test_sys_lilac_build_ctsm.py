#!/usr/bin/env python3

"""System tests for lilac_build_ctsm

These tests do a lot of work (interacting with cime, etc.), and thus take relatively long
to run.
"""

import unittest
import tempfile
import shutil
import os

from ctsm.path_utils import add_cime_lib_to_path
from ctsm import unit_testing
from ctsm.lilac_build_ctsm import build_ctsm

_CIME_PATH = add_cime_lib_to_path(standalone_only=True)

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestSysBuildCtsm(unittest.TestCase):
    """System tests for lilac_build_ctsm"""

    def setUp(self):
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        self.assertTrue(os.path.isdir(self._tempdir))

        # Hack around a check in CIME: As of https://github.com/ESMCI/cime/pull/4228, If
        # NCAR_HOST is in the environment, CIME checks if the machine you're running on is
        # the one you have set up the case for. This is a problem for trying to test the
        # user-defined machine infrastructure on an NCAR machine. So bypass this CIME
        # check by temporarily removing NCAR_HOST from the environment if it was present.
        # (Then we restore it in the tearDown method.)
        if "NCAR_HOST" in os.environ:
            self._ncarhost = os.environ["NCAR_HOST"]
            del os.environ["NCAR_HOST"]
        else:
            self._ncarhost = None

        self._original_wd = os.getcwd()
        os.chdir(self._tempdir)

    def tearDown(self):
        """tear down"""
        os.chdir(self._original_wd)
        shutil.rmtree(self._tempdir, ignore_errors=True)
        if self._ncarhost is not None:
            os.environ["NCAR_HOST"] = self._ncarhost

    def test_buildSetup_userDefinedMachine_minimalInfo(self):
        """Get through the case.setup phase with a user-defined machine

        This tests that the xml files are created successfully and that they are
        compatible with cime's xml schemas. It also ensures that the creation of
        various directories goes smoothly.

        This version specifies a minimal amount of information
        """
        build_dir = os.path.join(self._tempdir, "ctsm_build")
        build_ctsm(
            cime_path=_CIME_PATH,
            build_dir=build_dir,
            compiler="gnu",
            no_build=True,
            os_type="linux",
            netcdf_path="/path/to/netcdf",
            esmf_mkfile_path="/path/to/esmf/lib/esmf.mk",
            max_mpitasks_per_node=16,
            gmake="gmake",
            gmake_j=8,
            no_pnetcdf=True,
        )
        self.assertTrue(os.path.isdir(build_dir))
        # the critical piece of this test is that the above command doesn't generate any
        # errors; however we also do some assertions below

        # ensure that inputdata directory was created
        inputdata = os.path.join(build_dir, "inputdata")
        self.assertTrue(os.path.isdir(inputdata))

    def test_buildSetup_userDefinedMachine_allInfo(self):
        """Get through the case.setup phase with a user-defined machine

        This tests that the xml files are created successfully and that they are
        compatible with cime's xml schemas. It also ensures that the creation of
        various directories goes smoothly.

        This version specifies all possible information
        """
        build_dir = os.path.join(self._tempdir, "ctsm_build")
        inputdata_path = os.path.realpath(os.path.join(self._tempdir, "my_inputdata"))
        os.makedirs(inputdata_path)
        self.assertTrue(os.path.isdir(inputdata_path))
        build_ctsm(
            cime_path=_CIME_PATH,
            build_dir=build_dir,
            compiler="gnu",
            no_build=True,
            os_type="linux",
            netcdf_path="/path/to/netcdf",
            esmf_mkfile_path="/path/to/esmf/lib/esmf.mk",
            max_mpitasks_per_node=16,
            gmake="gmake",
            gmake_j=8,
            pnetcdf_path="/path/to/pnetcdf",
            pio_filesystem_hints="gpfs",
            gptl_nano_timers=True,
            extra_fflags="-foo",
            extra_cflags="-bar",
            build_debug=True,
            build_with_openmp=True,
            inputdata_path=os.path.join(self._tempdir, "my_inputdata"),
        )
        self.assertTrue(os.path.isdir(build_dir))
        # the critical piece of this test is that the above command doesn't generate any
        # errors; however we also do some assertions below

        # ensure that inputdata directory is NOT created
        inputdata = os.path.join(build_dir, "inputdata")
        self.assertFalse(os.path.exists(inputdata))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
