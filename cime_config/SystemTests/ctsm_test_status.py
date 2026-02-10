
from CIME.test_status import *

CUPID_PHASE = "CUPID"

ALL_PHASES = [
    CREATE_NEWCASE_PHASE,
    XML_PHASE,
    SETUP_PHASE,
    NAMELIST_PHASE,
    SHAREDLIB_BUILD_PHASE,
    MODEL_BUILD_PHASE,
    SUBMIT_PHASE,
    RUN_PHASE,
    COMPARE_PHASE,
    BASELINE_PHASE,
    THROUGHPUT_PHASE,
    MEMCOMP_PHASE,
    MEMLEAK_PHASE,
    STARCHIVE_PHASE,
    CUPID_PHASE,
    GENERATE_PHASE,
]

# Extend the TestStatus class to include the CUPID phase and have CUPID phase within ALL_PHASES
class CTSMTestStatus(TestStatus):

    def __init__(self, case):
        super().__init__(case, ALL_PHASES)
        self._test_status = CTSMTestStatus(test_dir=case.get_value("CASEROOT"), test_name=self.case.get_value("CASEBASEID") )