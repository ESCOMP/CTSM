"""
This module contains the NeonSite class and class functions which are used in run_neon.py
"""

# Import libraries
import glob
import logging
import os
import re
import shutil
import sys
import time

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# -- import local classes for this script
# pylint: disable=wrong-import-position
from ctsm.site_and_regional.tower_site import TowerSite

# pylint: disable=wrong-import-position, import-error, unused-import, wrong-import-order
from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

logger = logging.getLogger(__name__)


# pylint: disable=too-many-instance-attributes
class NeonSite(TowerSite):
    """
    A class for encapsulating neon sites.
    """

    def __init__(self, name, start_year, end_year, start_month, end_month, finidat):
        super().__init__(name, start_year, end_year, start_month, end_month, finidat)

    def build_base_case(
        self, cesmroot, output_root, res, compset, overwrite=False, setup_only=False
    ):
        user_mods_dirs = [
            os.path.join(self.cesmroot, "cime_config", "usermods_dirs", "NEON", self.name)
        ]
        case_path = super().build_base_case(cesmroot, output_root, res, compset, user_mods_dirs)

        return case_path

    def get_batch_query(self, case):
        return super().get_batch_query(case)

    # pylint: disable=too-many-statements
    def run_case(
        self,
        base_case_root,
        run_type,
        prism,
        run_length,
        user_version,
        overwrite=False,
        setup_only=False,
        no_batch=False,
        rerun=False,
        experiment=False,
    ):
        """
        Run case.

        Args:
        self
        base_case_root: str, opt
            file path of base case
        run_type: str, opt
            transient, post_ad, or ad case, default transient
        prism: bool, opt
            if True, use PRISM precipitation, default False
        run_length: str, opt
            length of run, default '4Y'
        user_version: str, opt
            default 'latest'
        overwrite: bool, opt
            default False
        setup_only: bool, opt
            default False; if True, set up but do not run case
        no_batch: bool, opt
            default False
        rerun: bool, opt
            default False
        experiment: str, opt
            name of experiment, default False
        """
        user_mods_dirs = [
            os.path.join(self.cesmroot, "cime_config", "usermods_dirs", "NEON", self.name)
        ]
        tower_type = "NEON"
        super().run_case(
            base_case_root, run_type, prism, run_length, tower_type, user_mods_dirs, user_version
        )

    def set_ref_case(self, case):
        super().set_ref_case(case)
        return True  ### Check if super returns false, if this will still return True?

    def modify_user_nl(self, case_root, run_type, rundir):
        # TODO: include neon-specific user namelist lines, using this as just an example currently
        site_lines = [
            """hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC',
                                 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO',"""
        ]
        super().modify_user_nl(case_root, run_type, rundir, site_lines)
