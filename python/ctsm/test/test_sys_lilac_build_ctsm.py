#!/usr/bin/env python

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
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_buildSetup_userDefinedMachine_minimalInfo(self):
        """Get through the case.setup phase with a user-defined machine

        This tests that the xml files are created successfully and that they are
        compatible with cime's xml schemas. It also ensures that the creation of
        various directories goes smoothly.

        This version specifies a minimal amount of information
        """
        build_dir = os.path.join(self._tempdir, 'ctsm_build')
        build_ctsm(cime_path=_CIME_PATH,
                   build_dir=build_dir,
                   compiler='gnu',
                   no_build=True,
                   os_type='linux',
                   netcdf_path='/path/to/netcdf',
                   esmf_lib_path='/path/to/esmf/lib',
                   max_mpitasks_per_node=16,
                   gmake='gmake',
                   gmake_j=8,
                   no_pnetcdf=True)
        # the critical piece of this test is that the above command doesn't generate any
        # errors; however we also do some assertions below

        # ensure that inputdata directory was created
        inputdata = os.path.join(build_dir, 'inputdata')
        self.assertTrue(os.path.isdir(inputdata))

    def test_buildSetup_userDefinedMachine_allInfo(self):
        """Get through the case.setup phase with a user-defined machine

        This tests that the xml files are created successfully and that they are
        compatible with cime's xml schemas. It also ensures that the creation of
        various directories goes smoothly.

        This version specifies all possible information
        """
        build_dir = os.path.join(self._tempdir, 'ctsm_build')
        inputdata_path = os.path.realpath(os.path.join(self._tempdir, 'my_inputdata'))
        os.makedirs(inputdata_path)
        build_ctsm(cime_path=_CIME_PATH,
                   build_dir=build_dir,
                   compiler='gnu',
                   no_build=True,
                   os_type='linux',
                   netcdf_path='/path/to/netcdf',
                   esmf_lib_path='/path/to/esmf/lib',
                   max_mpitasks_per_node=16,
                   gmake='gmake',
                   gmake_j=8,
                   pnetcdf_path='/path/to/pnetcdf',
                   pio_filesystem_hints='gpfs',
                   gptl_nano_timers=True,
                   extra_fflags='-foo',
                   extra_cflags='-bar',
                   build_debug=True,
                   build_with_openmp=True,
                   inputdata_path=os.path.join(self._tempdir, 'my_inputdata'))
        # the critical piece of this test is that the above command doesn't generate any
        # errors; however we also do some assertions below

        # ensure that inputdata directory is NOT created
        inputdata = os.path.join(build_dir, 'inputdata')
        self.assertFalse(os.path.exists(inputdata))

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
