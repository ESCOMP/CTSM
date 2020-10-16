"""
Implementation of the CIME SSP test.  This class inherits from SystemTestsCommon

This is a CLM specific test:
Verifies that spinup works correctly
this test is only valid for CLM compsets with CLM45 or CLM50
(1) do an initial spin test
    - set CLM_ACCELERATED_SPINUP to on
    - write restarts at the end of the run, turn on short term archiving
    - turn MOSART off
(2) do a hybrid non-spinup simulation run
    - start from the restart files generated in (1)
    - turn MOSART off
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
import shutil, glob, os

logger = logging.getLogger(__name__)

class SSP(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SSP system test
        """
        SystemTestsCommon.__init__(self, case)
        rof = self._case.get_value("COMP_ROF")
        expect(rof == "mosart", "ERROR: SSP test requires that ROF component be mosart")

    def build_phase(self, sharedlib_only=False, model_only=False):
        self._case.set_value("MOSART_MODE", "NULL")
        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def run_phase(self):
        caseroot = self._case.get_value("CASEROOT")
        orig_case = self._case
        orig_casevar = self._case.get_value("CASE")

        # clone the main case to create ref1
        clone_path = "{}.ref1".format(caseroot)
        if os.path.exists(clone_path):
            shutil.rmtree(clone_path)
        clone = self._case.create_clone(clone_path, keepexe=True)

        # determine run lengths needed below
        stop_nf = self._case.get_value("STOP_N")
        stop_n1 = int(stop_nf / 2)
        stop_n2 = stop_nf - stop_n1

        #-------------------------------------------------------------------
        # (1) do a spinup run in the main case in the cloned ref case
        #      (short term archiving is on)
        #-------------------------------------------------------------------
        os.chdir(clone_path)
        self._set_active_case(clone)

        logger.info("startup: doing a {} {} 00000 seconds startup run".format(stop_n1, stop_nf))
        logger.info("  writing restarts at end of run")
        logger.info("  short term archiving is on ")

        with clone:
            clone.set_value("CLM_ACCELERATED_SPINUP", "on")
            clone.set_value("STOP_N",stop_n1)

        dout_sr = clone.get_value("DOUT_S_ROOT")
        # No history files expected, set suffix=None to avoid compare error
        self._skip_pnl = False
        self.run_indv(suffix=None, st_archive=True)

        #-------------------------------------------------------------------
        # (2) do a hybrid, non-spinup run in orig_case
        #-------------------------------------------------------------------
        os.chdir(caseroot)
        self._set_active_case(orig_case)

        refdate = run_cmd_no_fail(r'ls -1dt {}/rest/*-00000* | head -1 | sed "s/-00000.*//" | sed "s/^.*rest\///"'.format(dout_sr))
        refsec = "00000"

        # obtain rpointer files and necessary restart files from short term archiving directory
        rundir = self._case.get_value("RUNDIR")

        rest_path = os.path.join(dout_sr, "rest", "{}-{}".format(refdate, refsec))

        for item in glob.glob("{}/*{}*".format(rest_path, refdate)):
            link_name = os.path.join(rundir, os.path.basename(item))
            if os.path.islink(link_name) and os.readlink(link_name) == item:
                # Link is already set up correctly: do nothing
                # (os.symlink raises an exception if you try to replace an
                # existing file)
                pass
            else:
                os.symlink(item, link_name)

        for item in glob.glob("{}/*rpointer*".format(rest_path)):
            shutil.copy(item, rundir)

        self._case.set_value("CLM_ACCELERATED_SPINUP", "off")
        self._case.set_value("RUN_TYPE", "hybrid")
        self._case.set_value("GET_REFCASE", False)
        self._case.set_value("RUN_REFCASE", "{}.ref1".format(orig_casevar))

        self._case.set_value("RUN_REFDATE", refdate)
        self._case.set_value("STOP_N", stop_n2)
        self._case.set_value("DOUT_S", False)
        self._case.flush()

        # do the restart run (short term archiving is off)
        self.run_indv()
