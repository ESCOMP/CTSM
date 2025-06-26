"""
CTSM-specific test that first runs the subset_data point tool and then ensures
that CTSM does not fail using the just-generated input files
"""

import os
import sys
import logging
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.XML.standard_module_setup import *

_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)
from ctsm.subset_data import (  # pylint: disable=wrong-import-position
    main as subset_data,
)

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

    def build_phase(self, sharedlib_only=False, model_only=False):

        # Where the output files will be saved
        out_dir = os.path.join(self._get_caseroot(), "subset_data_output")
        usermods_dir = os.path.join(out_dir, "user_mods")

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
            "--create-datm",
            "--datm-syr",
            "2000",
            "--datm-eyr",
            "2000",
            "--create-user-mods",
            "--outdir",
            out_dir,
            "--user-mods-dir",
            usermods_dir,
            "--overwrite",
        ]
        subset_data()

        self._case.apply_user_mods([usermods_dir])
        super().build_phase(sharedlib_only, model_only)
