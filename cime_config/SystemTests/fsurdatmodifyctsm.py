"""
CTSM-specific test that first runs the fsurdat_modifier tool and then ensures
that the CTSM does not fail using the just-generated modified fsurdat file
"""

import os
import re
import subprocess
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)

class FSURDATMODIFYCTSM(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

        if not os.path.exists(os.path.join(self._get_caseroot(),
            'done_FSURDATMODIFYCTSM_setup.txt')):
            # Create out-of-the-box lnd_in to obtain fsurdat_in
            case.create_namelists(component='lnd')
            # If fsurdat_in does not exist, download it from the server
            case.check_all_input_data()

            lnd_in_path = os.path.join(self._get_caseroot(), 'CaseDocs/lnd_in')
            with open (lnd_in_path,'r') as lnd_in:
                for line in lnd_in:
                    fsurdat_in = re.match(r" *fsurdat *= *'(.*)'", line)
                    if fsurdat_in:
                        self._fsurdat_in = fsurdat_in.group(1)
                        break

            self._fsurdat_out = os.path.join(self._get_caseroot(), 'fsurdat.nc')
            self._ctsm_root = self._case.get_value( 'COMP_ROOT_DIR_LND')
            self._cfg_file_path = os.path.join(self._get_caseroot(),
                                               'modify_fsurdat.cfg')

            logger.info("  create config file to modify")
            self._create_config_file()
            logger.info("  run modify_fsurdat")
            self._run_modify_fsurdat()
            logger.info("  modify user_nl files")
            self._modify_user_nl()
            with open('done_FSURDATMODIFYCTSM_setup.txt', 'w') as fp:
                pass

    def _create_config_file(self):
        cfg_template_path = os.path.join(self._ctsm_root,
            'tools/modify_input_files/modify_fsurdat_template.cfg')

        with open (self._cfg_file_path,'w') as cfg_out:
            with open (cfg_template_path,'r') as cfg_in:
                for line in cfg_in:
                    if re.match(r' *fsurdat_in *=', line):
                        line = 'fsurdat_in = {}'.format(self._fsurdat_in)
                    elif re.match(r' *fsurdat_out *=', line):
                        line = 'fsurdat_out = {}'.format(self._fsurdat_out)
                    elif re.match(r' *idealized *=', line):
                        line = 'idealized = True'
                    cfg_out.write(line)


    def _run_modify_fsurdat(self):
        tool_path = os.path.join(self._ctsm_root,
                                 'tools/modify_input_files/fsurdat_modifier')

        self._case.load_env(reset=True)
        conda_env = ". "+self._get_caseroot()+"/.env_mach_specific.sh; "
        # Preprend the commands to get the conda environment for python first
        conda_env += self._get_conda_env()
        # Source the env
        try:
            subprocess.run( conda_env+"python3 "+tool_path+" "+self._cfg_file_path, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            print("ERROR while getting the conda environment and/or ")
            print("running the fsurdat_modifier tool: ")
            print("(1) If your ctsm_pylib environment is out of date or you ")
            print("have not created the ctsm_pylib environment, yet, you may ")
            print("get past this error by running ./py_env_create ")
            print("in your ctsm directory and trying this test again. ")
            print("(2) If conda is not available, install and load conda, ")
            print("run ./py_env_create, and then try this test again. ")
            print("(3) If (1) and (2) are not the issue, then you may be ")
            print("getting an error within the fsurdat_modifier tool itself. ")
            print("Default error message: ")
            print(error.output)
        except:
            print("ERROR trying to run fsurdat_modifier tool.")
            raise

    def _modify_user_nl(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "fsurdat = '{}'".format(self._fsurdat_out))

    def _get_conda_env(self):
        #
        # Add specific commands needed on different machines to get conda available
        # Use semicolon here since it's OK to fail
        #
        # Execute the module unload/load when "which conda" fails
        # eg on cheyenne
        try:
            subprocess.run( "which conda", shell=True, check=True)
            conda_env = " "
        except subprocess.CalledProcessError:
            # Remove python and add conda to environment for cheyennne
            conda_env = "module unload python; module load conda;"

        # Activate the python environment
        conda_env += " conda activate ctsm_pylib"
        # End above to get to actual command
        conda_env += " && "

        return( conda_env )
