"""
Implementation of the CIME LCISO (Land Carbon Isotope) test.
This is a CTSM specific test:
Verifies turning on carbon isotopes doesn't change answers
(1) do a run with Carbon isotopes off (suffix base)
(2) add C13 and C14 carbon isotopes on with their time-series (suffix cisoallon)
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)

class LCISO(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'cisoallon',
                                       run_one_description = 'carbon isotopes off',
                                       run_two_description = 'c13 and c14 isotopes on as well as C isotope time series')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_c13=F, use_c14=F")

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_c13=.true.,use_c14=.true.,use_c13_timeseries=.true.,use_c14_bombspike=.true.," + \
                                "hist_fexcl1='C13_AR','C13_GPP','C13_HR','C13_NBP','C13_SOILC_vr','C13_TOTECOSYSC'," + \
                                "'C13_TOTLITC','C13_TOTSOMC','C13_TOTVEGC','C14_AR','C14_GPP','C14_HR','C14_NBP'," + \
                                "'C14_SOILC_vr','C14_TOTECOSYSC','C14_TOTLITC','C14_TOTSOMC','C14_TOTVEGC'")


