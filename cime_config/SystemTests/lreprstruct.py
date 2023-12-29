"""
Implementation of the CIME LREPRSTRUCT (Land Reproductive Structure) test.

This is a CTSM specific test: Verifies that we can allocate and use a crop reproductive
structure pool, and that answers are identical to a run with a single reproductive grain
pool. This is useful until we have tests that include AgSys, which will exercise this
capability more completely.

(1) do a run with a second grain pool and two reproductive structure pools, with all
    allocation going to the second reproductive structure pool
(2) do a standard run only a single grain pool, but with no crop seed replenishment

The reason for having the new reproductive structure pools in the first case (rather than
the other way around) is that this results in having a few extra fields present in the
baselines.

"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)


class LREPRSTRUCT(SystemTestsCompareTwo):
    def __init__(self, case):
        SystemTestsCompareTwo.__init__(
            self,
            case,
            separate_builds=False,
            run_two_suffix="grain1",
            run_one_description="use a reproductive structure pool",
            run_two_description="use a single grain pool",
            ignore_fieldlist_diffs=True,
        )

    def _case_one_setup(self):
        # We don't really need a second grain pool for this test, but we might as well do
        # this to further exercise the looping over different reproductive components.
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="for_testing_use_second_grain_pool=.true.",
        )
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="for_testing_use_repr_structure_pool=.true.",
        )

        # Replace any GRAIN outputs with the same outputs for REPRODUCTIVE1 and REPRODUCTIVE2
        user_nl_clm_path = os.path.join(self._get_caseroot(), "user_nl_clm")
        with open(user_nl_clm_path) as f:
            user_nl_clm_text = f.read()
        for grain_output in re.findall("GRAIN\w*", user_nl_clm_text):
            user_nl_clm_text = user_nl_clm_text.replace(
                grain_output,
                grain_output.replace("GRAIN", "REPRODUCTIVE1")
                + "', '"
                + grain_output.replace("GRAIN", "REPRODUCTIVE2"),
            )
        with open(user_nl_clm_path, "w") as f:
            f.write(user_nl_clm_text)

    def _case_two_setup(self):
        # This is needed in the nearly-standard case to prevent grain from being used to
        # replenish crop seed deficits, thus making grain act like the reproductive
        # structure pools. (It wouldn't hurt to do this in case one as well, but it
        # shouldn't be needed there, since we shouldn't have any grain there anyway.)
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="for_testing_no_crop_seed_replenishment=.true.",
        )
