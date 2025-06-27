"""
CTSM-specific test that first runs the subset_data point tool and then ensures
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


class SUBSETDATAPOINT(SystemTestsCommon):
    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

        # Check the test setup
        if not self._case.get_value("LND_GRID") == "CLM_USRDAT":
            raise RuntimeError("SUBSETDATAPOINT tests require resolution CLM_USRDAT")
        if "serial" not in self._case.get_value("MPILIB"):
            raise RuntimeError("SUBSETDATAPOINT tests require a serial MPILIB")
        if "BGC-CROP" not in self._case.get_value("COMPSET"):
            raise RuntimeError("SUBSETDATAPOINT tests require a BGC-CROP compset")

    def build_phase(self, sharedlib_only=False, model_only=False):

        # Where the output files will be saved
        out_dir = os.path.join(self._get_caseroot(), "subset_data_output")
        usermods_dir = os.path.join(out_dir, "user_mods")

        # Import subset_data. Do it here rather than at top because otherwise the import will
        # be attempted even during RUN phase.
        # pylint: disable=wrong-import-position,import-outside-toplevel
        from ctsm.subset_data import main as subset_data

        # Run the tool
        lat = 45.402252
        lon = -92.798085
        sys.argv = [
            "tools/site_and_regional/subset_data",
            "point",
            "--lat",
            str(lat),
            "--lon",
            str(lon),
            "--create-surface",
            "--crop",
            "--create-landuse",
            "--surf-year",
            "1850",
            "--create-datm",
            "--datm-syr",
            "1901",
            "--datm-eyr",
            "1901",
            "--create-user-mods",
            "--outdir",
            out_dir,
            "--user-mods-dir",
            usermods_dir,
            "--overwrite",
        ]
        subset_data()

        # Required so that CTSM doesn't fail
        user_nl_clm_path = os.path.join(usermods_dir, "user_nl_clm")
        with open(user_nl_clm_path, "a", encoding="utf-8") as user_nl_clm:
            user_nl_clm.write("\ncheck_dynpft_consistency = .false.")

        # Apply the user mods
        self._case.flush(flushall=True)
        apply_user_mods(self._get_caseroot(), usermods_dir)
        self._case.read_xml()

        # Do the build
        super().build_phase(sharedlib_only, model_only)
