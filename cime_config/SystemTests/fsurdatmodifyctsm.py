"""
CTSM-specific test to ...
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
            # If fsurdat_in does not exist, download it from the
            # server
            case.check_all_input_data()

            lnd_in_path = os.path.join(self._get_caseroot(), 'CaseDocs/lnd_in')
            with open (lnd_in_path,'r') as lnd_in:
                for line in lnd_in:
                    fsurdat_in = re.match(r" *fsurdat *= *'(.*)'", line)
                    if fsurdat_in:
                        self._fsurdat_in = fsurdat_in.group(1)
                        break

            self._fsurdat_out = os.path.join(
                self._get_caseroot(), 'fsurdat.nc')
            self._ctsm_root = self._case.get_value(
                'COMP_ROOT_DIR_LND')
            self._cfg_file_path = os.path.join(
                self._get_caseroot(), 'modify_fsurdat.cfg')

            self._create_config_file()
            self._run_modify_fsurdat()
            self._modify_user_nl()
            with open('done_FSURDATMODIFYCTSM_setup.txt', 'w') as fp:
                pass

    def _create_config_file(self):
        cfg_template_path = os.path.join(self._ctsm_root,
            'tools/modify_fsurdat/modify_template.cfg')

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
                                 'tools/modify_fsurdat/fsurdat_modifier')

        # Need to specify a specific python version that has the required dependencies
        python_path = _get_python_path()

        subprocess.check_call([python_path, tool_path, self._cfg_file_path])

    def _modify_user_nl(self):
        append_to_user_nl_files(caseroot = self._get_caseroot(),
                                component = "clm",
                                contents = "fsurdat = '{}'".format(self._fsurdat_out))

def _get_python_path():
    """Get path to ncar_pylib's python on cheyenne

    This is needed because we need a python environment that includes xarray and its
    dependencies. This is currently hard-coded for cheyenne until we come up with a robust
    way in CIME of ensuring the correct python environment is loaded.

    """
    out = subprocess.check_output(['/glade/u/apps/opt/ncar_pylib/ncar_pylib', '-l'],
                                  universal_newlines=True)

    # First look for a loaded ('L') python
    path = _find_path_from_pylib_output(out, 'L')
    # If no loaded python found, look for a default ('D') python
    if path is None:
        path = _find_path_from_pylib_output(out, 'D')

    if path is None:
        raise RuntimeError('No python found')

    return os.path.join(path, 'bin', 'python')

def _find_path_from_pylib_output(ncar_pylib_output, char):
    """Given line-by-line output from ncar_pylib, return the path to python if found

    Args:
    - ncar_pylib_output: line-by-line output from ncar_pylib
    - char: the character to look for in the leading parenthetical expression (typically 'L' or 'D')

    Returns a path to python, or None if not found
    """
    # The line of interest looks like the following (for char = 'L'):
    # (L) ... /path/to/python
    regex = r'\(' + char + r'\).* (/\S+)'
    for line in ncar_pylib_output.splitlines():
        match_line = re.match(regex, line)
        if match_line:
            return match_line.group(1)

    return None
