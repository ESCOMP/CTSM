from CIME.SystemTests.sms import SMS
from ctsm_test_status import *


class SMSCUPID(SMS):
    def __init__(self, case):
        super().__init__(case)
        self._test_status = CTSMTestStatus(
            test_dir=self._case.get_value("CASEROOT"),
            test_name=self._case.get_value("CASEBASEID"),
        )

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
        self._phase_modifying_call(CUPID_PHASE, self._cupid_case_test)

    def case_test_cupid(self, testdir="cupid_test"):
        # create the run directory testdir
        if os.path.exists(testdir):
            logger.info("Removing existing test directory {}".format(testdir))
            shutil.rmtree(testdir)
        # Check that the CUPid postprocessing directories and config file exist
        cupid_dir = os.path.join(self._case.get_value("CASEROOT"), "cupid-postprocessing" )
        notebooks_dir = os.path.join(cupid_dir, "compute_notebooks")
        data_dir = os.path.join(cupid_dir, "temp_data")
        for dir in [cupid_dir, notebooks_dir, data_dir]:
           expect( os.isdir(dir),
                   "CUPiD postprocessing directory {} does not exist".format(dir) )
        cupid_config = os.path.join(cupid_dir, "config.yml")
        expect( os.isfile(cupid_config),
                "CUPiD config file {} does not exist".format(cupid_config) )
        # TODO: Add more checks about more files that should exist


        # TODO: Populate the testdir with data files, config files and notebooks from the cupid-postprocessing directory

        # TODO: Save various files to the baseline directory to use for BASELINE comparison

        return True

    def _cupid_case_test(self):
        # For the st_archiver this test is under the case object
        # Here we create it in this object, but probably should be moved to under the case object in CIME
        result = self.test_cupid()
        with self._test_status:
            if result:
                self._test_status.set_status(CUPID_PHASE, TEST_PASS_STATUS)
            else:
                self._test_status.set_status(CUPID_PHASE, TEST_FAIL_STATUS)
