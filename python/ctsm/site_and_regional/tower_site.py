"""
This module includes the definition for the TowerSite class,
which has NeonSite and Plumber2Site child classes. This class defines common
functionalities that are in both NeonSite and Plumber2Site classes.
"""
# -- Import libraries

# -- standard libraries
import os.path
import glob
import logging
import re
import shutil
import sys
import time

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# -- import local classes for this script
# pylint: disable=wrong-import-position
from ctsm.site_and_regional.base_site import BaseSite

# pylint: disable=wrong-import-position, import-error, unused-import, wrong-import-order
from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

logger = logging.getLogger(__name__)


class TowerSite(BaseSite):
    """
    Parent class to NeonSite and Plumber2Site classes.
    ...
    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self, name, start_year, end_year, start_month, end_month, finidat):
        """
        Initializes TowerSite with the given arguments.

        Parameters
        ----------
        """
        super().__init__(name, start_year, end_year, start_month, end_month, finidat)

    def __str__(self):
        """
        Converts ingredients of the TowerSite to string for printing.
        """
        return "{}\n{}".format(
            str(self.__class__),
            "\n".join(
                (
                    "{} = {}".format(str(key), str(self.__dict__[key]))
                    for key in sorted(self.__dict__)
                )
            ),
        )
