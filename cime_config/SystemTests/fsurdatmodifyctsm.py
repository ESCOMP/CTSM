"""
CTSM-specific test that first runs the fsurdat_modifier tool and then ensures
that the CTSM does not fail using the just-generated modified fsurdat file
"""

import os
import re
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

# For calling fsurdat_modifier
from argparse import Namespace

_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)

logger = logging.getLogger(__name__)


class FSURDATMODIFYCTSM(SystemTestsCommon):
    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

        if not os.path.exists(
            os.path.join(self._get_caseroot(), "done_FSURDATMODIFYCTSM_setup.txt")
        ):
            # Create out-of-the-box lnd_in to obtain fsurdat_in
            case.create_namelists(component="lnd")
            # If fsurdat_in does not exist, download it from the server
            case.check_all_input_data()

            lnd_in_path = os.path.join(self._get_caseroot(), "CaseDocs/lnd_in")
            with open(lnd_in_path, "r") as lnd_in:
                for line in lnd_in:
                    fsurdat_in = re.match(r" *fsurdat *= *'(.*)'", line)
                    if fsurdat_in:
                        self._fsurdat_in = fsurdat_in.group(1)
                        break

            self._fsurdat_out = os.path.join(self._get_caseroot(), "fsurdat.nc")
            self._ctsm_root = self._case.get_value("COMP_ROOT_DIR_LND")
            self._cfg_file_path = os.path.join(self._get_caseroot(), "modify_fsurdat.cfg")

            logger.info("  create config file to modify")
            self._create_config_file()
            logger.info("  run modify_fsurdat")
            self._run_modify_fsurdat()
            logger.info("  modify user_nl files")
            self._modify_user_nl()
            with open("done_FSURDATMODIFYCTSM_setup.txt", "w") as fp:
                pass

    def _create_config_file(self):
        cfg_template_path = os.path.join(
            self._ctsm_root, "tools/modify_input_files/modify_fsurdat_template.cfg"
        )

        with open(self._cfg_file_path, "w") as cfg_out:
            with open(cfg_template_path, "r") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *fsurdat_in *=", line):
                        line = "fsurdat_in = {}".format(self._fsurdat_in)
                    elif re.match(r" *fsurdat_out *=", line):
                        line = "fsurdat_out = {}".format(self._fsurdat_out)
                    elif re.match(r" *idealized *=", line):
                        line = "idealized = True"
                    cfg_out.write(line)

    def _run_modify_fsurdat(self):
        fsurdat_modifier_args = Namespace(
            cfg_path=self._cfg_file_path,
            debug=False,
            fsurdat_in="UNSET",
            fsurdat_out="UNSET",
            overwrite=False,
            silent=False,
            verbose=False,
        )
        from ctsm.modify_input_files.fsurdat_modifier import fsurdat_modifier

        fsurdat_modifier(fsurdat_modifier_args)

    def _modify_user_nl(self):
        append_to_user_nl_files(
            caseroot=self._get_caseroot(),
            component="clm",
            contents="fsurdat = '{}'".format(self._fsurdat_out),
        )
