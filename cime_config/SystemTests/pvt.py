"""
FATES land use potential vegetation spin up + transient test

This is a FATES specific test:

1) conduct a spinup with use_fates_potentialveg on
    - write restart file
    - use CLM_ACCELERATED_SPINUP?
2) run a transient landuse case with use_fates_lupft
    - start from the restart file generated in (1)
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
import shutil, glob, os

logger = logging.getLogger(__name__)


class PVT(SystemTestsCommon):
    def __init__(self, case):
        SystemTestsCommon.__init__(self, case)

        # Do not allow PVT to be run with certain testmods
        # Should this be targeted to a specific testmod for simplicity for now?
        # Technically this could be run with the luh fates_harvest_modes
        error_message = None
        casebaseid = self._case.get_value("CASEBASEID")
        casebaseid = casebaseid.split("-")[-1]
        if casebaseid[0:10] != "FatesLUPFT":
            error_message = f"Only call PVT with testmod FatesLUPFT. {casebaseid} selected."

        # Only allow to run if resolution is 4x5 for now
        # Other grid resolutions will be pre-processed and included in the namelist defaults at a future date.
        # Potentially we could generate these on the fly although doing so would result in increased build time
        lnd_grid = self._case.get_value("LND_GRID")
        if lnd_grid != "4x5":
            error_message = (
                f"PVT can currently only be run with 4x5 resolution. {lnd_grid} selected."
            )

        if error_message is not None:
            logger.error(error_message)
            raise RuntimeError(error_message)

    def run_phase(self):
        # -------------------------------------------------------------------
        # (1) Run FATES spin-up case in potential vegetation mode
        # -------------------------------------------------------------------
        orig_case = self._case
        orig_casevar = self._case.get_value("CASE")
        caseroot = self._case.get_value("CASEROOT")

        # Set the run start date based on the desired starting reference case year
        refcase_year = 1700
        stop_n_pveg = 5
        startyear_pveg = refcase_year - stop_n_pveg

        # clone the main case to create spinup case
        logger.info("PVT log:  cloning setup")
        clone_path = "{}.potveg".format(caseroot)
        if os.path.exists(clone_path):
            shutil.rmtree(clone_path)
        logger.info("PVT log:  cloning")
        clone = self._case.create_clone(clone_path, keepexe=True)
        logger.info("PVT log: cloning complete")

        # setup the clone case
        os.chdir(clone_path)
        self._set_active_case(clone)

        # set the clone case values
        with clone:
            clone.set_value("CLM_ACCELERATED_SPINUP", "off")
            clone.set_value("STOP_N", stop_n_pveg)
            clone.set_value("STOP_OPTION", "nyears")
            clone.set_value("RUN_STARTDATE", "{}-01-01".format(startyear_pveg))

        # Modify the spin up case to use the potential vegetation mode.
        # Checks for incompatible cases and necessary mapping files are
        # handled in the build case.
        # Turn off fates_harvest_mode for the spin up.

        logger.info("PVT log:  modify user_nl_clm file for spin up run")
        added_content = ["use_fates_potentialveg = .true.", "fates_harvest_mode = 'no_harvest'"]
        append_to_user_nl_files(clone_path, "clm", added_content)

        # Run the spin up case
        # As per SSP test:
        # "No history files expected, set suffix=None to avoid compare error"
        logger.info("PVT log:  starting spin-up run")
        dout_sr = clone.get_value("DOUT_S_ROOT")
        self._skip_pnl = False
        self.run_indv(suffix=None, st_archive=True)

        # -------------------------------------------------------------------
        # (2) Run FATES transient case using restart file from spin-up
        # -------------------------------------------------------------------
        os.chdir(caseroot)
        self._set_active_case(orig_case)

        # Copy restart files from spin up to the transient case run directory
        # obtain rpointer files and necessary restart files from short term archiving directory
        rundir = self._case.get_value("RUNDIR")

        refdate = str(refcase_year) + "-01-01-00000"
        rest_path = os.path.join(dout_sr, "rest", "{}".format(refdate))

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

        # Update run case settings
        self._case.set_value("CLM_ACCELERATED_SPINUP", "off")
        self._case.set_value("RUN_TYPE", "hybrid")
        self._case.set_value("GET_REFCASE", False)
        self._case.set_value("RUN_REFCASE", "{}.potveg".format(orig_casevar))
        self._case.set_value("RUN_REFDATE", "{}-01-01".format(refcase_year))
        self._case.set_value("RUN_STARTDATE", "{}-01-01".format(refcase_year))
        self._case.set_value("DOUT_S", False)
        self._case.flush()

        # do the restart run (short term archiving is off)
        self.run_indv()
