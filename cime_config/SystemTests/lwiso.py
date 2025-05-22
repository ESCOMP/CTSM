"""
Implementation of the CIME LWISO (Land Water Isotope) test.

This is a CTSM specific test:
Verifies turning on water isotopes doesn't change answers
(1) do a run with water isotopes on
(2) do a run with water isotopes off

The reason for having water isotopes on in the first case (rather than the other way
around) is that this results in having water isotope quantities present in the baselines.

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)


class LWISO(SystemTestsCompareTwo):
    def __init__(self, case):
        SystemTestsCompareTwo.__init__(
            self,
            case,
            separate_builds=False,
            run_two_suffix="nowiso",
            run_one_description="water isotopes on",
            run_two_description="water isotopes off",
            ignore_fieldlist_diffs=True,
        )

    def _case_one_setup(self):
        # BUG(wjs, 2019-07-30, ESCOMP/ctsm#495) We currently can't turn on actual water
        # isotopes in a multi-timestep test, so we're setting
        # enable_water_tracer_consistency_checks rather than enable_water_isotopes;
        # eventually, though, we should change this to the latter. (See
        # <https://github.com/ESCOMP/ctsm/issues/495#issuecomment-516619853>.)
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="enable_water_tracer_consistency_checks=.true.",
        )

    def _case_two_setup(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="enable_water_tracer_consistency_checks=.false.",
        )
