from CIME.SystemTests.sms import SMS
from ctsm_test_status import *


class SMSAR(SMS):
    def setup_indv(
        self, clean=False, test_mode=False, reset=False, keep=False, disable_git=False
    ):
        """
        Perform an individual setup
        """
        super().setup_indv(
            clean=clean,
            test_mode=test_mode,
            reset=reset,
            keep=keep,
            disable_git=disable_git,
        )
        # Make sure the st_archiver is turned on, as CUPiD only runs after it runs
        self._case.set_value("DOUT_S", True)

    def run_indv(
        self,
        suffix="base",
        st_archive=True,
        cupid=True,
        submit_resubmits=None,
        keep_init_generated_files=False,
    ):
        """
        Perform an individual run. Raises an EXCEPTION on fail.

        Just add the CUPiD phase after the standard run.
        """
        super().run_indv(
            suffix=suffix,
            st_archive=st_archive,
            submit_resubmits=submit_resubmits,
            keep_init_generated_files=keep_init_generated_files,
        )
