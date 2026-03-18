"""
Parent class for CTSM-specific tests that first run the subset_data tool and then ensure
that CTSM does not fail using the just-generated input files
"""

import os
import sys
import logging
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.user_mod_support import apply_user_mods
from CIME.XML.standard_module_setup import *

# In case we need to import subset_data later
_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)

logger = logging.getLogger(__name__)


class SUBSETDATASHARED(SystemTestsCommon):
    def __init__(self, case, subset_data_cmd):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

        # Check the test setup
        if not self._case.get_value("LND_GRID") == "CLM_USRDAT":
            raise RuntimeError("SUBSETDATA tests require resolution CLM_USRDAT")
        if "serial" not in self._case.get_value("MPILIB"):
            raise RuntimeError("SUBSETDATA tests require a serial MPILIB")
        if "BGC-CROP" not in self._case.get_value("COMPSET"):
            raise RuntimeError("SUBSETDATA tests require a BGC-CROP compset")

        # Add standard subset_data arguments
        out_dir = os.path.join(self._get_caseroot(), "subset_data_output")
        self.usermods_dir = os.path.join(out_dir, "user_mods")
        self.subset_data_cmd = subset_data_cmd + [
            "--create-user-mods",
            "--outdir",
            out_dir,
            "--user-mods-dir",
            self.usermods_dir,
            "--overwrite",
        ]

    def setup_phase(self, clean=False, test_mode=False, reset=False, keep=False, disable_git=False):
        """
        override default SMS system test setup phase
        """

        # Default SMS behavior
        self.setup_indv(
            clean=clean,
            test_mode=test_mode,
            reset=reset,
            keep=keep,
            disable_git=disable_git,
        )

        # Run subset_data, if needed.
        # It's needed during SETUP and/or NLCOMP phases if comparing/generating a baseline because
        # the namelist comparison will require generating a namelist, and that will fail if we
        # haven't specified our custom fsurdat and other stuff. By calling self._run_subset_data()
        # only if the usermods directory doesn't yet exist, we avoid it being called every time the
        # test class is initialized (which happens, e.g., in RUN phase).
        if not os.path.exists(self.usermods_dir):
            self._run_subset_data()

    def _run_subset_data(self):
        """
        Run subset_data
        """
        # Import subset_data. Do it here rather than at top because otherwise the import will
        # be attempted even during RUN phase.
        # pylint: disable=wrong-import-position,import-outside-toplevel
        from ctsm.subset_data import main as subset_data

        # Run subset_data
        sys.argv = self.subset_data_cmd
        subset_data()

        # Required so that CTSM doesn't fail
        user_nl_clm_path = os.path.join(self.usermods_dir, "user_nl_clm")
        with open(user_nl_clm_path, "a", encoding="utf-8") as user_nl_clm:
            user_nl_clm.write("\ncheck_dynpft_consistency = .false.\n")

        # Apply the user mods
        self._case.flush(flushall=True)
        apply_user_mods(self._get_caseroot(), self.usermods_dir)
        self._case.read_xml()
