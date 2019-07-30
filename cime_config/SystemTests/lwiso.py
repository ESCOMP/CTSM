"""
Implementation of the CIME LWISO (Land Water Isotope) test.

This is a CTSM specific test:
Verifies turning on water isotopes doesn't change answers
(1) do a run with water isotopes off (suffix base)
(2) add water isotopes

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)

class LWISO(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'wiso',
                                       run_one_description = 'water isotopes off',
                                       run_two_description = 'water isotopes on',
                                       ignore_fieldlist_diffs = True)

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "enable_water_isotopes=.false.")

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "enable_water_isotopes=.true.")

