#!/usr/bin/env python

"""System tests for build_ctsm

These tests do a lot of work (interacting with cime, etc.), and thus take relatively long
to run.
"""

import unittest
import tempfile
import shutil
import os

from ctsm.path_utils import add_cime_lib_to_path
from ctsm import unit_testing
from ctsm.build_ctsm import build_ctsm

_CIME_PATH = add_cime_lib_to_path(standalone_only=True)

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

class TestSysBuildCtsm(unittest.TestCase):
    """System tests for build_ctsm"""

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
        build_ctsm(cime_path=_CIME_PATH,
                   build_dir=os.path.join(self._tempdir, 'ctsm_build'),
                   compiler='gnu',
                   no_build=True,
                   os_type='linux',
                   netcdf_path='/path/to/netcdf',
                   esmf_lib_path='/path/to/esmf/lib',
                   gmake='gmake',
                   gmake_j=8)
        # no assertions: test passes as long as the command doesn't generate any errors

    def test_buildSetup_userDefinedMachine_allInfo(self):
        """Get through the case.setup phase with a user-defined machine

        This tests that the xml files are created successfully and that they are
        compatible with cime's xml schemas. It also ensures that the creation of
        various directories goes smoothly.

        This version specifies all possible information
        """
        build_ctsm(cime_path=_CIME_PATH,
                   build_dir=os.path.join(self._tempdir, 'ctsm_build'),
                   compiler='gnu',
                   no_build=True,
                   os_type='linux',
                   netcdf_path='/path/to/netcdf',
                   esmf_lib_path='/path/to/esmf/lib',
                   gmake='gmake',
                   gmake_j=8,
                   pnetcdf_path='/path/to/pnetcdf',
                   pio_filesystem_hints='gpfs',
                   gptl_nano_timers=True,
                   extra_fflags='-foo',
                   extra_cflags='-bar')
        # no assertions: test passes as long as the command doesn't generate any errors

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
