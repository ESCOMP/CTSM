"""
This test passes if mksurfdata_esmf generates an fsurdat (surface dataset)
and the CTSM completes a simulation with this fsurdat file.

We test res = '10x15' because it uses a lower-res topography file instead of
the 1-km topography raw dataset. The 1-km file causes the test to run out of
memory on cheyenne.

Currently casper complains that `git -C` is not a valid option.
I added -C to the `git describe` in gen_mksurfdata_namelist for this
system test to work.
"""

import os
import sys
import subprocess
from datetime import datetime
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)


class MKSURFDATAESMF(SystemTestsCommon):
    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

        # Paths and strings needed throughout
        ctsm_root = self._case.get_value("COMP_ROOT_DIR_LND")
        self._tool_path = os.path.join(ctsm_root, "tools/mksurfdata_esmf")
        self._tool_bld = os.path.join(self._get_caseroot(), "tool_bld")
        time_stamp = datetime.today().strftime("%y%m%d")
        self._res = "10x15"  # see important comment in script's docstring
        self._model_yr = "1850"
        self._jobscript = os.path.join(
            self._get_caseroot(), "mksurfdataesmf_test_jobscript_single.sh"
        )
        self._fsurdat_namelist = os.path.join(
            self._get_caseroot(),
            f"surfdata_{self._res}_hist_{self._model_yr}_78pfts_c{time_stamp}.namelist",
        )
        self._fsurdat_nc = os.path.join(
            self._get_caseroot(),
            f"surfdata_{self._res}_hist_{self._model_yr}_78pfts_c{time_stamp}.nc",
        )
        self._TestStatus_log_path = os.path.join(self._get_caseroot(), "TestStatus.log")

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Build executable that will generate fsurdat
        Generate namelist for generating fsurdat
        Generate jobscript that runs executable
        Modify user_nl_clm to point to the generated fsurdat
        """
        # build_phase gets called twice:
        # - once with sharedlib_only = True and
        # - once with model_only = True
        # Call the following steps only once during the test but do not skip
        # if the test stops and gets restarted.
        if sharedlib_only:
            # Paths and strings
            build_script_path = os.path.join(self._tool_path, "gen_mksurfdata_build")
            nml_script_path = os.path.join(self._tool_path, "gen_mksurfdata_namelist")
            gen_jobscript_path = os.path.join(self._tool_path, "gen_mksurfdata_jobscript_single")
            gen_mksurfdata_namelist = f"{nml_script_path} --res {self._res} --start-year {self._model_yr} --end-year {self._model_yr}"

            if not os.path.exists(nml_script_path):
                sys.exit(f"ERROR The build naemlist script {nml_script_path} does NOT exist")

            if not os.path.exists(gen_jobscript_path):
                sys.exit(f"ERROR The jobscript script {gen_jobscript_path} does NOT exist")

            gen_mksurfdata_jobscript = (
                f"{gen_jobscript_path} --number-of-nodes 1 --tasks-per-node 64 --namelist-file "
                + f"{self._fsurdat_namelist} --bld-path {self._tool_bld} --jobscript-file {self._jobscript}"
            )
            if not os.path.exists(build_script_path):
                sys.exit(f"ERROR The build script {build_script_path} does NOT exist")

            # Rm tool_bld and build executable that will generate fsurdat
            try:
                subprocess.check_call(f"rm -rf {self._tool_bld}", shell=True)
            except subprocess.CalledProcessError as e:
                sys.exit(
                    f"{e} ERROR REMOVING {self._tool_bld}. DETAILS IN {self._TestStatus_log_path}"
                )
            try:
                subprocess.check_call(f"{build_script_path} --blddir {self._tool_bld}", shell=True)
            except subprocess.CalledProcessError as e:
                print(f"build directory = {self._tool_bld}\n")
                sys.exit(
                    f"{e} ERROR RUNNING {build_script_path} DETAILS IN {self._TestStatus_log_path}"
                )

            # Generate namelist for generating fsurdat (rm namelist if exists)
            if os.path.exists(self._fsurdat_namelist):
                os.remove(self._fsurdat_namelist)
            try:
                subprocess.check_call(gen_mksurfdata_namelist, shell=True)
            except subprocess.CalledProcessError as e:
                sys.exit(
                    f"{e} ERROR RUNNING {gen_mksurfdata_namelist}. DETAILS IN {self._TestStatus_log_path}"
                )

            # Generate jobscript that will run the executable
            if os.path.exists(self._jobscript):
                os.remove(self._jobscript)
            try:
                subprocess.check_call(gen_mksurfdata_jobscript, shell=True)
            except subprocess.CalledProcessError as e:
                sys.exit(
                    f"{e} ERROR RUNNING {gen_mksurfdata_jobscript}. DETAILS IN {self._TestStatus_log_path}"
                )
            # Change self._jobscript to an executable file
            subprocess.check_call(f"chmod a+x {self._jobscript}", shell=True)

        # Call this step only once even if the test stops and gets restarted.
        if not os.path.exists(os.path.join(self._get_caseroot(), "done_MKSURFDATAESMF_setup.txt")):
            # Modify user_nl_clm to point to the generated fsurdat
            self._modify_user_nl()
            with open("done_MKSURFDATAESMF_setup.txt", "w") as fp:
                pass

        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def run_phase(self):
        """
        Run executable to generate fsurdat
        Submit CTSM run that uses fsurdat just generated
        """

        # Run executable to generate fsurdat (rm fsurdat if exists)
        if os.path.exists(self._fsurdat_nc):
            os.remove(self._fsurdat_nc)
        try:
            subprocess.check_call(self._jobscript, shell=True)
        except subprocess.CalledProcessError as e:
            sys.exit(f"{e} ERROR RUNNING {self._jobscript}; details in {self._TestStatus_log_path}")

        # Submit CTSM run that uses fsurdat just generated
        self.run_indv()

    def _modify_user_nl(self):
        """
        Modify user_nl_clm to point to the generated fsurdat
        """
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="fsurdat = '{}'".format(self._fsurdat_nc)
            + "\n"
            + "convert_ocean_to_land = .true.",
        )
