"""
Implementation of the CTSM FUNIT test.

This "system" test runs CTSM's Fortran unit tests. We're abusing the system test
infrastructure to run these, so that a run of the test suite can result in the unit tests
being run as well.

Grid and compset are irrelevant for this test type.
"""

import os
from CIME.SystemTests.funit import FUNIT
from CIME.XML.standard_module_setup import *
from CIME.utils import append_testlog, get_cime_root

logger = logging.getLogger(__name__)

class FUNITCTSM(FUNIT):

    def __init__(self, case):
        FUNIT.__init__(self, case)

    def get_test_spec_dir(self):
        lnd_root = self._case.get_value("COMP_ROOT_DIR_LND")
        return os.path.join(lnd_root, "src")

    # TODO(wjs, 2018-10-28) For now, we're overriding run_phase with a slightly modified
    # version that has a fix for the path to the unit_test_tool. Once we can point to an
    # updated cime that has this fix in place, we should remove run_phase, relying on
    # FUNIT's run_phase. At that point, we should also remove the import of append_testlog
    # and get_cime_root, above.
    def run_phase(self):

        rundir   = self._case.get_value("RUNDIR")
        exeroot  = self._case.get_value("EXEROOT")
        mach     = self._case.get_value("MACH")

        log = os.path.join(rundir, "funit.log")
        if os.path.exists(log):
            os.remove(log)

        test_spec_dir = self.get_test_spec_dir()
        unit_test_tool = os.path.abspath(os.path.join(get_cime_root(),"scripts","fortran_unit_testing","run_tests.py"))
        args = "--build-dir {} --test-spec-dir {} --machine {}".format(exeroot, test_spec_dir, mach)
        stat = run_cmd("{} {} >& funit.log".format(unit_test_tool, args), from_dir=rundir)[0]

        append_testlog(open(os.path.join(rundir, "funit.log"), "r").read())

        expect(stat == 0, "RUN FAIL for FUNIT")

