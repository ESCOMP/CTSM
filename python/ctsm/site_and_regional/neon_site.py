"""
This module contains the NeonSite class and class functions which extend the tower_site class for
things that are specific just for NEON sites.
"""

# Import libraries
import logging
import os
import sys

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# -- import local classes for this script
# pylint: disable=wrong-import-position
from ctsm.site_and_regional.tower_site import TowerSite

# pylint: disable=wrong-import-position, import-error, unused-import, wrong-import-order
from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root
from ctsm.utils import abort

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

logger = logging.getLogger(__name__)


# pylint: disable=too-many-instance-attributes
class NeonSite(TowerSite):
    """
    A class for encapsulating neon sites.
    """

    def __init__(self, *args, **kwargs):
        super().__init__("NEON", *args, **kwargs)

    def modify_user_nl(self, case_root, run_type, rundir, site_lines=None):
        # TODO: include neon-specific user namelist lines, using this as just an example currently
        if site_lines is None:
            site_lines = [
                """hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC',
                                 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO',"""
            ]
        super().modify_user_nl(case_root, run_type, rundir, site_lines)
