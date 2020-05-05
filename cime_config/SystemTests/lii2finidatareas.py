"""
Implementation of the LII2FINIDATAREAS test.

This is similar to the LII test, but tests init_interp with mode
'use_finidat_areas'.

As with the standard LII test, this must be used in a configuration for
which we have a compatible out-of-the-box finidat file (so that the run
with use_init_interp = .false. runs successfully). In constrast to our
standard LII test (which uses glcMEC_spunup_1way), this one can use a
standard CISM2%NOEVOLVE configuration: we do *not* need to set
GLC_TWO_WAY_COUPLING=FALSE; in fact, it's a better test if we have
GLC_TWO_WAY_COUPLING=TRUE: with this mode of operation, areas should
match between the two runs.

It may not be totally necessary to have a system test covering this
init_interp_method: between unit tests and inline checks, we have quite
a bit of testing of this code. It's important to have at least one LII
test to make sure we're interpolating all the fields we're supposed to
interpolate, but we don't necessarily need a LII test with this mode of
operation as long as we have one with the general mode of operation.

However, this test is a good check that every point has a unique type
(because this mode of operation will fail if that's not the case); this
is something needed even for the 'general' init_interp method to work
correctly, but the 'general' mode won't catch problems in this regard.

One other reason why this test is useful is to cover the threading
directives in the set_single_match routine, since those aren't covered
by unit tests. So this test mod should be used in a test with threading.

To update the initial conditions (finidat) file for this test:

(1) Run the test; the 'base' case should run to completion even if the
no_interp test fails. (If the 'base' case fails, you may need to retry
with init_interp_method='general'.)

(2) Copy the finidat_interp_dest.nc file from the 'base' case to the inputdata
space. Rename this to be similar to the out-of-the-box finidat file
currently used by this test, but with a new creation date.

(3) Update namelist defaults to point to the new finidat file. If
updating the out-of-the-box file is not desired, then you could instead
point to this new finidat file with a user_nl_clm file in this testmod.
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

# We can import lii directly because the SystemTests directory has been
# added to sys.path.
#
# A cleaner and more general way to accomplish this would be: In cime:
# Rather than adding the SystemTests directory to sys.path and then
# importing an individual module from there: instead, allow each
# component to have a COMPNAME_pylib directory within its cime_config,
# and then have cime add each component's cime_config directory to
# sys.path and then do something like:
# import_module("COMPNAME_pylib.SystemTests.TESTNAME"). Then, for
# example, ctsm could access its own modules via "import
# ctsm_pylib.foo", or (in this case) "from ctsm_pylib.SystemTests.lii
# import LII".
from lii import LII

logger = logging.getLogger(__name__)

class LII2FINIDATAREAS(LII):

    def __init__(self, case):
        super(LII2FINIDATAREAS, self).__init__(case)

    def _case_one_setup(self):
        super(LII2FINIDATAREAS, self)._case_one_setup()
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "init_interp_method = 'use_finidat_areas'")
