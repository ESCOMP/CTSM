"""
Implementation of the FATES unit tests.

This "system" test runs FATES's Fortran unit tests. We're abusing the system test
infrastructure to run these, so that a run of the test suite can result in the unit tests
being run as well.

Grid and compset are irrelevant for this test type.
"""

import os
import sys
from CIME.SystemTests.funit import FUNIT
from CIME.XML.standard_module_setup import *
import systemtest_utils as stu

_FATES_TESTING_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "src", "fates", "testing"
)

logger = logging.getLogger(__name__)


class FUNITFATES(FUNIT):
    def __init__(self, case):
        FUNIT.__init__(self, case)

    def run_phase(self):
        tool_path = os.path.join(_FATES_TESTING_PYTHON, "run_unit_tests.py")
        build_dir = os.path.join("bld", "fates_unit_tests")
        cmd = f"{tool_path} -b {build_dir}"

        stu.run_python_script(
            self._get_caseroot(),
            "ctsm_pylib",
            cmd,
            tool_path,
        )
