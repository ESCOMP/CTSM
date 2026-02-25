"""Misc. shared utilities"""

import os
import sys

# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    INPUTDATA_PREFIX,
)


def convert_to_absolute_path(relative_path: str) -> str:
    """
    Convert a relative path to an absolute path.

    Args:
        relative_path: Relative path starting with OUR_PATH, or already absolute path

    Returns:
        Absolute path
    """
    # If the path is already absolute, return it as-is
    if os.path.isabs(relative_path):
        return relative_path

    # Otherwise, convert relative path to absolute
    return os.path.join(INPUTDATA_PREFIX, relative_path)
