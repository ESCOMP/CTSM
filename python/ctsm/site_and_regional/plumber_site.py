"""
This module contains the Plumber2Site class and class functions that extend the tower_site class for
things that are specific just for PLUMBER2 sites.
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

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

logger = logging.getLogger(__name__)


# pylint: disable=too-many-instance-attributes
class Plumber2Site(TowerSite):
    """
    A class for encapsulating plumber sites.
    """

    def __init__(self, *args, **kwargs):
        super().__init__("PLUMBER2", *args, **kwargs)

    def set_ref_case(self, case):
        super().set_ref_case(case)
        return True  ### Check if super returns false, if this will still return True?
