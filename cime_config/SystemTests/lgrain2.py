"""
Implementation of the CIME LGRAIN2 (Land 2-grain-pool) test.

This is a CTSM specific test: Verifies that we can allocate and use a second grain pool,
and that answers are identical to a run with a single grain pool. This is useful until we
have tests that include AgSys, which will exercise this capability more completely.

(1) do a run with a second grain pool
(2) do a standard run only a single grain pool

The reason for having the second grain pool in the first case (rather than the other way
around) is that this results in having a few extra fields present in the baselines.

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)


class LGRAIN2(SystemTestsCompareTwo):
    def __init__(self, case):
        SystemTestsCompareTwo.__init__(
            self,
            case,
            separate_builds=False,
            run_two_suffix="grain1",
            run_one_description="use a second grain pool",
            run_two_description="use a single grain pool",
            ignore_fieldlist_diffs=True,
        )

    def _case_one_setup(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="for_testing_use_second_grain_pool=.true.",
        )

    def _case_two_setup(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="for_testing_use_second_grain_pool=.false.",
        )
