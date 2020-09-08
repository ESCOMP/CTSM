"""
Implementation of the CIME LVG (Land Virtual Glacier) test.

This is a CLM specific test:
Verifies that adding virtual glacier columns doesn't change answers
(1) do a run with the standard set of virtual columns (suffix base)
(2) add virtual columns over Antarctica (suffix more_virtual)

This will only pass if there are no column or patch-level outputs in CLM's
history files.
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)

class LVG(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'more_virtual',
                                       run_one_description = 'standard set of virtual columns',
                                       run_two_description = 'add virtual columns over Antarctica')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "glacier_region_behavior = 'single_at_atm_topo', 'virtual', 'virtual', 'multiple'")

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "glacier_region_behavior = 'single_at_atm_topo', 'virtual', 'virtual', 'virtual'")

