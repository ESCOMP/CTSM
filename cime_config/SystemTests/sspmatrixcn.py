"""

CTSM only test to do the CN-matrix spinup procedure

This is a CLM specific test:
Verifies that spinup works correctly
this test is only valid for CLM compsets

Step 0: Run a AD cold-start with matrix and matrix spinup off
        Fast mode and fast-mode 2-loop spinup steps are now skipped
        These were labeled as Step 1 and Step 2.
Step 3: Run a slow-mode spinup
Step 4: matrix Spinup off
"""

import shutil, glob, os, sys
from datetime import datetime

if __name__ == "__main__":
    CIMEROOT = os.environ.get("CIMEROOT")
    if CIMEROOT is None:
        CIMEROOT = "../../cime"

    sys.path.append(os.path.join(CIMEROOT, "scripts", "lib"))
    sys.path.append(os.path.join(CIMEROOT, "scripts"))
else:
    from CIME.status import append_testlog

from CIME.case import Case
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.SystemTests.test_utils import user_nl_utils


logger = logging.getLogger(__name__)


class SSPMATRIXCN(SystemTestsCommon):

    # Class data
    nyr_forcing = 2
    # Get different integer multiples of the number of forcing years
    full = nyr_forcing
    twice = 2 * nyr_forcing
    thrice = 3 * nyr_forcing
    # Define the settings that will be used for each step
    steps = ["0-AD", "1-SASU", "2-norm"]
    desc = [
        "Accell-Decomp(AD)-coldstart",
        "slow-mode Semi-Analytic SpinUp(SASU)",
        "normal",
    ]
    runtyp = ["startup", "hybrid", "branch"]
    spin = ["on", "sasu", "off"]
    stop_n = [5, thrice, thrice]
    cold = [True, False, False]
    iloop = [-999, -999, -999]
    sasu = [-999, -999, -999]

    def __init__(self, case=None):
        """
        initialize an object interface to the SSPMATRIXCN system test
        """
        expect(
            len(self.steps) == len(self.sasu),
            "length of steps must be the same as sasu",
        )
        expect(
            len(self.steps) == len(self.spin),
            "length of steps must be the same as spin",
        )
        expect(
            len(self.steps) == len(self.desc),
            "length of steps must be the same as desc",
        )
        expect(
            len(self.steps) == len(self.cold),
            "length of steps must be the same as cold",
        )
        expect(
            len(self.steps) == len(self.runtyp),
            "length of steps must be the same as runtyp",
        )
        expect(
            len(self.steps) == len(self.iloop),
            "length of steps must be the same as iloop",
        )
        expect(
            len(self.steps) == len(self.stop_n),
            "length of steps must be the same as stop_n",
        )

        if __name__ != "__main__":
            SystemTestsCommon.__init__(self, case)
            ystart = int(self._case.get_value("DATM_YR_START"))
            yend = int(self._case.get_value("DATM_YR_END"))
            self.comp = self._case.get_value("COMP_LND")
        else:
            self._case = None
            self.comp = "clm"
            ystart = 2000
            yend = 2001

        for n in range(len(self.steps)):
            if n == 0:
                expect(self.cold[n] == True, "First step MUST be a cold-start")
                expect(self.runtyp[n] == "startup", "First step MUST be a startup")
            else:
                expect(self.cold[n] == False, "Other steps must NOT be a cold-start")
                expect(self.runtyp[n] != "startup", "Other steps MUST NOT be a startup")

            if self.spin[n] == "sasu":
                expect(self.cold[n] == False, "SASU step should NOT be a cold-start")
                if self.sasu[n] != -999:
                    expect(self.sasu[n] > 0, "SASU steps must set SASU cycle")
                    expect(
                        self.sasu[n] <= self.nyr_forcing,
                        "SASU cycles can't be greater than a full forcing cycle",
                    )

        expect(
            yend - ystart + 1 == self.nyr_forcing,
            "Number of years run over MUST correspond to nyr_forcing",
        )
        self._testname = "SSPMATRIX"

    def check_n(self, n):
        "Check if n is within range"
        expect(
            ((n >= 0) and (n < self.n_steps())),
            "Step number is out of range = " + str(n),
        )

    def __logger__(self, n=0):
        "Log info on this step"

        self.check_n(n)
        msg = "Step {}: {}: doing a {} run for {} years".format(
            self.steps[n], self.runtyp[n], self.desc[n], self.stop_n[n]
        )
        logger.info(msg)
        logger.info("  spinup type: {}".format(self.spin[n]))
        if __name__ != "__main__":
            append_testlog(msg)
        if n + 1 < self.n_steps():
            logger.info("  writing restarts at end of run")
            logger.info("  short term archiving is on ")

    def n_steps(self):
        "Total number of steps"

        return len(self.steps)

    def total_years(self):
        "Total number of years needed to do the full spinup"

        ysum = 0
        for nyr in self.stop_n:
            ysum = ysum + nyr

        return ysum

    def append_user_nl(self, caseroot, n=0):
        "Append needed settings to the user_nl files"

        self.check_n(n)
        # For all set output to yearly
        contents_to_append = "hist_nhtfrq = -8760"
        contents_to_append = contents_to_append + ", hist_mfilt = " + str(self.nyr_forcing)
        # For all but last step turn extra matrix output to off
        b4last = self.n_steps() - 1
        if n < b4last:
            contents_to_append = contents_to_append + ", hist_wrt_matrixcn_diag = .False."
        # For matrix spinup steps, set the matrix spinup and other variables associated with it
        if self.spin[n] == "sasu":
            contents_to_append = contents_to_append + ", nyr_forcing = " + str(self.nyr_forcing)
            if self.sasu[n] != -999:
                contents_to_append = contents_to_append + ", nyr_sasu = " + str(self.sasu[n])
            if self.iloop[n] != -999:
                contents_to_append = contents_to_append + ", iloop_avg = " + str(self.iloop[n])

        # For cold start, run with matrix off
        if self.cold[n]:
            contents_to_append = contents_to_append + ", use_matrixcn = .False."
            contents_to_append = contents_to_append + ", use_soil_matrixcn = .False."

        # Always append to the end
        user_nl_utils.append_to_user_nl_files(
            caseroot=caseroot, component=self.comp, contents=contents_to_append
        )

    def run_indv(self, nstep, st_archive=True):
        """
        Individual run of a given step
        """
        suffix = "step{}".format(self.steps[nstep])
        if isinstance(self._case, Case):
            super().run_indv(suffix, st_archive=True)
        else:
            caseroot = self._case.get_value("CASEROOT")
            dout_sr = self._case.get_value("DOUT_S_ROOT")
            rest_r = os.path.join(dout_sr, "rest")
            nyear = 1851 + nstep
            rundate = "%s-01-01-00000" % nyear
            restdir = os.path.join(rest_r, rundate)
            os.mkdir(restdir)
            rpoint = os.path.join(restdir, "rpointer.clm." + rundate)
            os.mknod(rpoint)
            rpoint = os.path.join(restdir, "rpointer.cpl." + rundate)
            os.mknod(rpoint)

    def run_phase(self):
        "Run phase"

        caseroot = self._case.get_value("CASEROOT")
        orig_case = self._case
        orig_casevar = self._case.get_value("CASE")

        # Get a clone of each step except the last one
        b4last = self.n_steps() - 1
        for n in range(b4last):
            #
            # Clone the main case, and get it setup for the next step
            #
            clone_path = "{}.step{}".format(caseroot, self.steps[n])
            if os.path.exists(clone_path):
                shutil.rmtree(clone_path)
            if n > 0:
                del clone
            self._set_active_case(orig_case)
            clone = self._case.create_clone(clone_path, keepexe=True)
            os.chdir(clone_path)
            self._set_active_case(clone)

            self.__logger__(n)

            with clone:
                clone.set_value("RUN_TYPE", self.runtyp[n])
                clone.set_value("STOP_N", self.stop_n[n])

                clone.set_value("CLM_ACCELERATED_SPINUP", self.spin[n])

                if self.cold[n]:
                    clone.set_value("CLM_FORCE_COLDSTART", "on")
                else:
                    clone.set_value("CLM_FORCE_COLDSTART", "off")

                self.append_user_nl(clone_path, n)

            dout_sr = clone.get_value("DOUT_S_ROOT")
            ninst = self._case.get_value("NINST")

            self._skip_pnl = False
            #
            # Start up from the previous case
            #
            rundir = clone.get_value("RUNDIR")
            with clone:
                if n > 0:
                    clone.set_value("GET_REFCASE", False)
                    expect("refcase" in locals(), "refcase was NOT previously set")
                    clone.set_value("RUN_REFCASE", refcase)
                    expect("refdate" in locals(), "refdate was NOT previously set")
                    clone.set_value("RUN_STARTDATE", refdate)
                    clone.set_value("RUN_REFDATE", refdate)
                    for item in glob.glob("{}/*{}*".format(rest_path, refdate)):
                        linkfile = os.path.join(rundir, os.path.basename(item))
                        if os.path.exists(linkfile):
                            os.remove(linkfile)
                        if not os.path.isdir(rundir):
                            os.makedirs(rundir)
                        os.symlink(item, linkfile)

                    # For a branch the cpl rpointer file needs to be handled
                    if self.runtyp[n] == "branch":

                        drvrest = "rpointer.cpl"
                        if ninst > 1:
                            drvrest += "_0001"
                        drvrest += rest_time

                        self._set_drv_restart_pointer(drvrest)
                        try:
                            shutil.copy(drvrest, rundir)
                        except shutil.SameFileError:
                            pass
            #
            # Run the case (Archiving on)
            #
            self._case.flush()
            self.run_indv(nstep=n, st_archive=True)

            #
            # Get the reference case from this step for the next step
            #
            refcase = clone.get_value("CASE")
            refdate = run_cmd_no_fail(
                r'ls -1dt {}/rest/*-00000* | head -1 | sed "s/-00000.*//" | sed "s/^.*rest\///"'.format(
                    dout_sr
                )
            )
            refsec = "00000"
            rest_path = os.path.join(dout_sr, "rest", "{}-{}".format(refdate, refsec))
            rest_time = "." + refdate + "-" + refsec

        #
        # Last step in original case
        #
        n = self.n_steps() - 1
        #
        # Setup the case to run from the previous clone step
        #
        os.chdir(caseroot)
        self._set_active_case(orig_case)
        self.__logger__(n)
        self._case.set_value("DOUT_S", False)
        self._case.set_value("RUN_TYPE", self.runtyp[n])
        self._case.set_value("STOP_N", self.stop_n[n])
        rundir = self._case.get_value("RUNDIR")
        self._case.set_value("GET_REFCASE", False)
        expect("refcase" in locals(), "refcase was NOT previously set")
        self._case.set_value("RUN_REFCASE", refcase)
        expect("refdate" in locals(), "refdate was NOT previously set")
        self._case.set_value("RUN_REFDATE", refdate)
        self._case.set_value("RUN_STARTDATE", refdate)
        for item in glob.glob("{}/*{}*".format(rest_path, refdate)):
            linkfile = os.path.join(rundir, os.path.basename(item))
            if os.path.exists(linkfile):
                os.remove(linkfile)
            expect(True, os.path.exists(item), "expected file does NOT exist = " + item)
            os.symlink(item, linkfile)

        # For a branch the cpl rpointer file needs to be handled
        if self.runtyp[n] == "branch":

            drvrest = "rpointer.cpl"
            if ninst > 1:
                drvrest += "_0001"
            drvrest += rest_time

            self._set_drv_restart_pointer(drvrest)
            try:
                shutil.copy(os.path.join(rest_path, drvrest), rundir)
            except shutil.SameFileError:
                pass

        self.append_user_nl(clone_path, n)
        #
        # Don't need to set COLDSTART or ACCEL_SPINUP
        #

        #
        # Run the case (short term archiving is off)
        #
        self._case.flush()
        self.run_indv(nstep=n, st_archive=False)
