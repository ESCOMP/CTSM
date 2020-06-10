"""
Implementation of the CIME LILACSMOKE (LILAC smoke) test.

This is a CTSM-specific test. It tests the building and running of CTSM via LILAC. Grid
and compset are ignored.
"""

import os

from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import run_cmd_no_fail, append_testlog
from CIME.build import post_build
from CIME.test_status import GENERATE_PHASE, BASELINE_PHASE, TEST_PASS_STATUS
from CIME.XML.standard_module_setup import *

logger = logging.getLogger(__name__)

class LILACSMOKE(SystemTestsCommon):

    def __init__(self, case):
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        if not sharedlib_only:
            caseroot = self._case.get_value('CASEROOT')
            lndroot = self._case.get_value('COMP_ROOT_DIR_LND')
            exeroot = self._case.get_value('EXEROOT')
            build_dir = os.path.join(caseroot, 'lilac_build')
            script_path = os.path.abspath(os.path.join(lndroot, 'lilac', 'build_ctsm'))
            logs = []

            # We only run the initial build command if the build_dir doesn't exist
            # yet. This is to support rebuilding the test case. (The first time through,
            # the build_dir won't exist yet; subsequent times, it will already exist, so
            # we skip to the rebuild command.)
            if not os.path.isdir(build_dir):
                machine = self._case.get_value('MACH')
                compiler = self._case.get_value('COMPILER')
                # TODO(wjs, 2020-06-10) Add --build-debug if the test is a debug test
                cmd = '{script_path} {build_dir} --machine {machine} --compiler {compiler}'.format(
                    script_path=script_path,
                    build_dir=build_dir,
                    machine=machine,
                    compiler=compiler)
                append_testlog(cmd)
                run_cmd_no_fail(cmd, arg_stdout='build_ctsm.bldlog', combine_output=True, from_dir=exeroot)
                logfile = os.path.join(exeroot, 'build_ctsm.bldlog')
                logs.append(logfile)
                with open(logfile) as lf:
                    append_testlog(lf.read())

            # We call the build script with --rebuild even for an initial build. This is
            # so we make sure to test the code path for --rebuild. (This is also needed if
            # the user rebuilds the test case, in which case this will be the only command
            # run, since the build_dir will already exist.)
            cmd = '{script_path} {build_dir} --rebuild'.format(
                script_path=script_path,
                build_dir=build_dir)
            append_testlog(cmd)
            run_cmd_no_fail(cmd, arg_stdout='rebuild_ctsm.bldlog', combine_output=True, from_dir=exeroot)
            logfile = os.path.join(exeroot, 'rebuild_ctsm.bldlog')
            logs.append(logfile)
            with open(logfile) as lf:
                append_testlog(lf.read())

            post_build(self._case, logs, build_complete=True)

    def run_phase(self):
        # TODO(wjs, 2020-06-10) Fill this in
        pass
