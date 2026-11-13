"""
Implementation of the CIME LVG (Land Virtual Glacier) test.

This is a CLM specific test:
Verifies that adding virtual glacier columns doesn't change answers
(1) do a run with the standard set of virtual columns (suffix base)
(2) add virtual columns over Antarctica (suffix more_virtual)

This will only pass if there are no column or patch-level outputs in CLM's
history files.

In PR #4234, which updated from gx (POP) to tx (MOM) grids, this test failed in COMPARE_base_more_virtual until we added the fexcl statement below for PCT_GLC_MEC in both case1 and case2. Bill Sacks posted in this PR that the MOM grid had a grid cell on the south tip of South America with 0% glacier that likely was newly present in the land mask, and with extra virtual glacier (case2) it was given non-zero percentages in each glacier column due to falling within the Antarctica mask. The PCT_GLC_MEC difference was valid and as such could be excluded from the test.
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)


class LVG(SystemTestsCompareTwo):
    def __init__(self, case):
        SystemTestsCompareTwo.__init__(
            self,
            case,
            separate_builds=False,
            run_two_suffix="more_virtual",
            run_one_description="standard set of virtual columns",
            run_two_description="add virtual columns over Antarctica",
        )

    def _case_one_setup(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="glacier_region_behavior = 'single_at_atm_topo', 'UNSET', 'virtual', 'multiple'",
        )
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="hist_fexcl1 = 'PCT_GLC_MEC'",
        )

    def _case_two_setup(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="glacier_region_behavior = 'single_at_atm_topo', 'UNSET', 'virtual', 'virtual'",
        )
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="hist_fexcl1 = 'PCT_GLC_MEC'",
        )
