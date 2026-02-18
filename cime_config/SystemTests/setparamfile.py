"""
CTSM-specific test that first runs the set_paramfile tool and then ensures that CTSM does not fail
using the just-generated parameter file
"""

import os
import sys
import logging
import re
from CIME.SystemTests.system_tests_common import SystemTestsCommon

# In case we need to import set_paramfile later
_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)

logger = logging.getLogger(__name__)


class SETPARAMFILE(SystemTestsCommon):
    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

        # Create out-of-the-box lnd_in to obtain paramfile
        case.create_namelists(component="lnd")

        # Find the paramfile to modify
        lnd_in_path = os.path.join(self._get_caseroot(), "CaseDocs", "lnd_in")
        self._paramfile_in = None
        with open(lnd_in_path, "r", encoding="utf-8") as lnd_in:
            for line in lnd_in:
                paramfile_in = re.match(r" *paramfile *= *'(.*)'", line)
                if paramfile_in:
                    self._paramfile_in = paramfile_in.group(1)
                    break
        if not self._paramfile_in:
            raise RuntimeError(f"paramfile not found in {lnd_in_path}")

        # Get the output file
        self.paramfile_out = os.path.join(self._get_caseroot(), "paramfile.nc")

        # Define set_paramfile command
        self.set_paramfile_cmd = [
            "set_paramfile",
            "-i",
            self._paramfile_in,
            "-o",
            self.paramfile_out,
            # Change two parameters for one PFT
            "-p",
            "needleleaf_deciduous_boreal_tree",
            "rswf_min=0.35",
            "rswf_max=0.7",
        ]

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Run set_paramfile and then build the model
        """

        # Run set_paramfile.
        # build_phase gets called twice:
        # - once with sharedlib_only = True and
        # - once with model_only = True
        # Because we only need set_paramfile run once, we only do it for the sharedlib_only call.
        # We could also check for the existence of the set_paramfile outputs, but that might lead to
        # a situation where the user expects set_paramfile to be called but it's not. Better to run
        # unnecessarily (e.g., if you fixed some FORTRAN code and just need to rebuild).
        if sharedlib_only:
            self._run_set_paramfile()

        # Do the build
        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def _run_set_paramfile(self):
        """
        Run set_paramfile
        """
        # Import set_paramfile. Do it here rather than at top because otherwise the import will
        # be attempted even during RUN phase.
        # pylint: disable=wrong-import-position,import-outside-toplevel
        from ctsm.param_utils.set_paramfile import main as set_paramfile

        # Run set_paramfile
        sys.argv = self.set_paramfile_cmd
        set_paramfile()

        # Append
        user_nl_clm_path = os.path.join(self._get_caseroot(), "user_nl_clm")
        with open(user_nl_clm_path, "a", encoding="utf-8") as user_nl_clm:
            user_nl_clm.write(f"paramfile = '{self.paramfile_out}'\n")
