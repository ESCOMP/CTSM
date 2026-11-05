from CIME.SystemTests.sms import SMS

class SMSPOSTPROC(SMS):
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
        self._case.set_value("RUN_POSTPROCESSING", True)
