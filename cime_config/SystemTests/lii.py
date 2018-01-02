"""
Implementation of the CIME LII test.

This is a CLM specific test:
Verifies that interpolation of initial conditions onto an identical
configuration gives identical results:
(1) do a run with use_init_interp true (suffix base)
(2) do a run with use_init_interp false (suffix no_interp)

It is more intuitive to think of the no_interp test as the "base". However, we
do the use_init_interp=true test first to facilitate updating initial conditions
whenever this is necessary, as documented below.

The LII test needs to point to an initial conditions file that is compatible
with the given model configuration. Thus, the pointed-to initial conditions file
needs to be updated whenever surface datasets are changed, or the land-mask is
changed, or an imporant change is made to model physics (for example where new
fields are added to the restart file). The procedure for updating the initial
conditions files used by the LII test is as follows:

(1) Run the LII test; the 'base' case should run to completion even if the
no_interp test fails.

(2) Copy the finidat_interp_dest.nc file from the 'base' case to the inputdata
space. Rename this to be similar to the name of the file pointed to in this
test's user_nl_clm file, but with a new creation date.

(3) Update this test's user_nl_clm file (in the appropriate testmods directory)
to point to the new finidat file.
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)

class LII(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'no_interp',
                                       run_one_description = 'use_init_interp set to true',
                                       run_two_description = 'use_init_interp set to false')

    def _case_one_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_init_interp = .true.")

    def _case_two_setup(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "use_init_interp = .false.")

