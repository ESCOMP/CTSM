"""
This module includes the definition for a parent class for the TowerSite class,
which has NeonSite and Plumber2Site child classes. This class defines common structures
that are in both NeonSite and Plumber2Site classes. The TowerSite class will define common
functionalities that are in both NeonSite and Plumber2Site classes.
"""
# -- Import libraries

# -- standard libraries
import os.path
import sys
import logging

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position, import-error, unused-import, wrong-import-order
from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root

logger = logging.getLogger(__name__)


class BaseSite:
    """
    Parent class to TowerSite
    ...
    Attributes
    ----------

    """

    def __init__(self, name, start_year, end_year, start_month, end_month, finidat):
        """
        Initializes BaseSite with the given arguments.

        Parameters
        ----------
        """
        self.name = name
        self.start_year = int(start_year)
        self.end_year = int(end_year)
        self.start_month = int(start_month)
        self.end_month = int(end_month)
        self.cesmroot = path_to_ctsm_root()
        self.finidat = finidat

    def __str__(self):
        """
        Converts ingredients of the BaseSite to string for printing.
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

    # pylint: disable=no-self-use
    def get_batch_query(self, case):
        """
        Function for querying the batch queue query command for a case, depending on the
        user's batch system.

        Args:
        case:
            case object
        """

        if case.get_value("BATCH_SYSTEM") == "none":
            return "none"
        return case.get_value("batch_query")
