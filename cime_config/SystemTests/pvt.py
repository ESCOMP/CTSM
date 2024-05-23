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
from systemtest_utils import read_user_nl_clm
import shutil, glob, os

logger = logging.getLogger(__name__)

class PVT(SystemTestsCommon):
    def __init__(self,case):
        SystemTestsCommon.__init__(self,case)

        # Do not allow PVT to be run with certain testmods
        # Should this be targeted to a specific testmod for simplicity for now?
        # Technically this could be run with the luh fates_harvest_modes
        error_message = None
        casebaseid = self._case.get_value("CASEBASEID")
        casebaseid = casebaseid.split("-")[-1]
        if casebaseid[0:10] != "FatesLUPFT":
            error_message = (f"Only call PVT with testmod FatesLUPFT. {casebaseid} selected.")

        # Only allow to run if resolution is 4x5 for now
        # Eventually we could set this up to generate the necessary land use x pft mapping
        # on the fly, although this would also require generating the land use timeseries
        # regridding on the fly which is a more time consuming endevour currently
        lnd_grid = self._case.get_value("LND_GRID")
        if lnd_grid != "4x5":
            error_message = (f"PVT can currently only be run with 4x5 resolution. {lnd_grid} selected.")

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
        self._setup_all()

        # Modify the spin up case to use the potential vegetation mode.
        # Checks for incompatible cases and necessary mapping files are
        # handled in the build case.
        # use_fates_lupft should be set to true in the testmod
        # save off the harvest mode for reinstating later
        found, hmode = self.read_user_nl_clm('fates_harvest_mode')
        logger.info("PVT log:  modify user_nl_clm file for spin up run")
        # TODO: remove fates_harvest_mode if found
        self._append_to_user_nl_clm(
            [
                "use_fates_potentialveg = .true.",
            ]
        )

        stop_n_pveg = 1
        with clone:
            # clone.set_value("CLM_ACCELERATED_SPINUP", "on")
            clone.set_value("STOP_N", stop_n_pveg)

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

        # logger.info("PVT log:  modify user_nl_clm file for transient run")
        # self._append_to_user_nl_clm(
        #     [
        #         "use_fates_potentialveg = .true.",
        #     ]
        # )


        self._case.set_value("CLM_ACCELERATED_SPINUP", "off")
        self._case.set_value("RUN_TYPE", "hybrid")
        self._case.set_value("GET_REFCASE", False)
        self._case.set_value("RUN_REFCASE", "{}.potveg".format(orig_casevar))
        self._case.set_value("RUN_REFDATE", "1700-01-01")
        self._case.set_value("DOUT_S", False)
        self._case.flush()

        # do the restart run (short term archiving is off)
        self.run_indv()
