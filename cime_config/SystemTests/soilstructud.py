"""
Implementation of the CIME SOILSTRUCTUD test.

This is a CLM specific test:
Verifies that a simulation that points to user_nl_ctsm containing
soil_layerstruct_userdefined_nlevsoi = 4
soil_layerstruct_userdefined = 0.1d0,0.3d0,0.6d0,1.0d0,1.0d0
gives bfb same results as one that points to user_nl_ctsm containing
soil_layerstruct_predefined = '4SL_2m'

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)


class SOILSTRUCTUD(SystemTestsCompareTwo):
    def __init__(self, case):
        SystemTestsCompareTwo.__init__(
            self,
            case,
            separate_builds=False,
            run_two_suffix="4SL_2m",
            run_one_description="soil_layerstruct_userdefined",
            run_two_description="soil_layerstruct_predefined",
        )

    def _case_one_setup(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="soil_layerstruct_userdefined_nlevsoi = 4,soil_layerstruct_userdefined = 0.1d0,0.3d0,0.6d0,1.0d0,1.0d0",
        )

    def _case_two_setup(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="soil_layerstruct_predefined = '4SL_2m'",
        )
