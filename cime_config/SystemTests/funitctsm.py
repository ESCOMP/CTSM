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

logger = logging.getLogger(__name__)

class FUNITCTSM(FUNIT):

    def __init__(self, case):
        FUNIT.__init__(self, case)

    def get_test_spec_dir(self):
        lnd_root = self._case.get_value("COMP_ROOT_DIR_LND")
        return os.path.join(lnd_root, "src")
