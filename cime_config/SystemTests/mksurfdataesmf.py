"""
"""

import os
import re
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

        # Steps for here or build phase (e.g. see ssp.py)
        # 1) build executable self._build_script_path
        # 2) generate namelist self._nml_script_path
        # 3) generate jobscript self._job_script_path
        # 4) change user_nl_clm as done in fsurdatmodify...

    def run_phase(self):
        """
        Submit jobscript
        Submit CTSM run that uses fsurdat just generated
        """
        # Submit jobscript
        # MKSURFDATAESMF_P144x1.f19_g17.I1850Clm50BgcCrop.cheyenne_intel.clm-default
        # with -np 144 stopped; presumably ran out of memory
        # MKSURFDATAESMF_P48x3.f19_g17.I1850Clm50BgcCrop.cheyenne_intel.clm-default
        # with -np 48 (144 gave error) stopped; presumably ran out of memory
        # MKSURFDATAESMF_P144x1.f10_f10_mg37.I1850Clm50BgcCrop.cheyenne_intel.clm-default
        # with -np 144 worked!
        subprocess.check_call("mpiexec_mpt -np 144 /glade/work/slevis/git/mksurfdata_toolchain/tools/mksurfdata_esmf/bld/mksurfdata < /glade/work/slevis/git/mksurfdata_toolchain/tools/mksurfdata_esmf/surfdata_10x15_hist_78pfts_CMIP6_1850_c220517.namelist", shell=True)
        # Submit CTSM run that uses fsurdat just generated
        self.run_indv()
