"""
Module for handling NaN fill values in CTSM input files.
"""

import os
import sys

_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
print(_CTSM_PYTHON)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.no_nans_in_inputs.constants import ATTR, NEW_FILLVALUES_FILE, USER_REQ_DELETE, XML_FILE
