"""
This test passes if mksurfdata_esmf generates an fsurdat (surface dataset)
and the CTSM completes a simulation with this fsurdat file.
"""

import os
import subprocess
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

        # Paths and command strings
        self._ctsm_root = self._case.get_value('COMP_ROOT_DIR_LND')
        self._tool_path = os.path.join(self._ctsm_root, 'tools/mksurfdata_esmf')

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Build executable that will generate fsurdat
        Generate namelist for generating fsurdat
        SKIP Generate jobscript
        Submit jobscript to generate fsurdat
        Modify user_nl_clm to point to the generated fsurdat
        """
        # These steps need only happen once
        if not os.path.exists(os.path.join(self._get_caseroot(),
            'done_MKSURFDATAESMF_setup.txt')):
            # Paths and command strings
            self._rm_bld_dir = f"rm -rf {self._tool_path}/bld"
            self._fsurdat_out_prefix = os.path.join(self._get_caseroot(), 'surfdata_10x15_hist_78pfts_CMIP6_1850_c220517.')
            self._build_script_path = os.path.join(self._tool_path,
                'gen_mksurfdata_build.sh')
            self._nml_script_path = os.path.join(self._tool_path,
                'gen_mksurfdata_namelist.py')
            self._gen_mksurfdata_namelist = f"{self._nml_script_path} --res 10x15 --start-year 1850 --end-year 1850"

            # Build executable
            subprocess.check_call(self._rm_bld_dir, shell=True)
            subprocess.check_call(self._build_script_path, shell=True)

            # Generate namelist self._nml_script_path
            # TODO Error: No such file or directory: './gen_mksurfdata_namelist.xml' ...probably need explicit path for this to work
            # subprocess.check_call(self._gen_mksurfdata_namelist, shell=True)

            # Modify user_nl_clm to point to the generated fsurdat
            self._modify_user_nl()
            with open('done_MKSURFDATAESMF_setup.txt', 'w') as fp:
                pass

        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def run_phase(self):
        """
        Submit jobscript to generate fsurdat
        Submit CTSM run that uses fsurdat just generated
        """
        # Paths and command strings
        self._executable_path = os.path.join(self._tool_path, 'bld/mksurfdata')
        # self._mpiexec_mpt_cmd = f"mpiexec_mpt -np 144 {self._executable_path} < {self._fsurdat_out_prefix}namelist"  # TODO Rm nxt line when this is ready
        self._mpiexec_mpt_cmd = f"mpiexec_mpt -np 144 {self._executable_path} < {self._tool_path}/surfdata_10x15_hist_78pfts_CMIP6_1850_c220517.namelist"

        # Submit jobscript to generate fsurdat
        subprocess.check_call(self._mpiexec_mpt_cmd, shell=True)

        # Submit CTSM run that uses fsurdat just generated
        self.run_indv()

    def _modify_user_nl(self):
        """
        Modify user_nl_clm to point to the generated fsurdat
        """
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "fsurdat = '{}nc'".format(self._fsurdat_out_prefix))
