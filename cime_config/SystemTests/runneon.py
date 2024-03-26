"""
CTSM-specific test that uses the run_neon.py methods for setting up, building, and running a NEON
site
"""

import os
import sys
from CIME.SystemTests.system_tests_common import SystemTestsCommon

_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)
import ctsm.site_and_regional.run_neon as rn  # python: disable=wrong-import-position


class RUNNEON(SystemTestsCommon):
    # pylint: disable=too-many-instance-attributes
    """
    CTSM-specific test that uses the run_neon.py methods for setting up, building, and running a
    NEON site
    """

    def __init__(self, case):
        # initialize an object interface to the SMS system test
        SystemTestsCommon.__init__(self, case)

        caseroot = self._case.get_value("CASEROOT")
        argv = [
            "run_neon",
            "--neon-sites",
            "ABBY",
            "--output-root",
            caseroot,
        ]

        (
            self._cesmroot,
            site_list,
            self._output_root,
            self._run_type,
            self._experiment,
            self._prism,
            self._overwrite,
            self._run_length,
            self._base_case_root,
            run_from_postad,
            self._setup_only,
            self._no_batch,
            self._rerun,
            self._user_version,
            available_list,
            self._res,
            self._compset,
        ) = rn.setup(description="", argv=argv)

        if not available_list:
            raise RuntimeError("available_list not defined")
        self._neon_site = None
        for neon_site in available_list:
            if neon_site.name in site_list:
                self._neon_site = neon_site
                break
        if self._neon_site is None:
            raise RuntimeError(
                f"available_list does not contain any member of site_list: {site_list}"
            )

        if run_from_postad:
            self._neon_site.finidat = None

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Builds the case
        """
        if not model_only:
            self._base_case_root = self._neon_site.build_base_case(
                self._cesmroot,
                self._output_root,
                self._res,
                self._compset,
                self._overwrite,
                self._setup_only,
            )

    def run_phase(self):
        """
        Runs the case
        """
        self._neon_site.run_case(
            self._base_case_root,
            self._run_type,
            self._prism,
            self._run_length,
            self._user_version,
            self._overwrite,
            self._setup_only,
            self._no_batch,
            self._rerun,
            self._experiment,
        )
