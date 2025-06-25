"""
CTSM-specific test that first runs the subset_data point tool and then ensures
that CTSM does not fail using the just-generated input files
"""

import os
import sys
import logging
import subprocess
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.XML.standard_module_setup import *

_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)
from ctsm.subset_data import (
    main as subset_data,
)  # pylint: disable=wrong-import-position

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

        # Where the output files will be saved
        self.out_dir = os.path.join(self._get_caseroot(), "subset_data_output")
        self.usermods_dir = os.path.join(self.out_dir, "user_mods")

        # Run the tool
        self.lat = 45.402252
        self.lon = -92.798085
        sys.argv = [
            "tools/site_and_regional/subset_data",
            "point",
            "--lat",
            str(self.lat),
            "--lon",
            str(self.lon),
            "--create-surface",
            "--create-datm",
            "--datm-syr",
            "2000",
            "--datm-eyr",
            "2000",
            "--create-user-mods",
            "--outdir",
            self.out_dir,
            "--user-mods-dir",
            self.usermods_dir,
        ]
        subset_data()

        # Append files in usermods_dir to equivalents in top level of case
        for file in os.listdir(self.usermods_dir):
            # Get the full path to source and destination files
            source_file_path = os.path.join(self.usermods_dir, file)
            dest_file_path = os.path.join(self._get_caseroot(), file)

            # Read the content of the source file
            with open(source_file_path, "r", encoding="utf-8") as source_file:
                content_to_append = source_file.read()

            # Append the content to the destination file
            with open(dest_file_path, "a", encoding="utf-8") as dest_file:
                dest_file.write("\n\n\n" + content_to_append)


    def build_phase(self, sharedlib_only=False, model_only=False):

        # Run shell_commands
        subprocess.check_call("./shell_commands", shell=True)

        super().build_phase(sharedlib_only, model_only)
