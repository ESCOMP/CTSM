"""
Tests for the no_nans_in_inputs module.
"""

import os
import sys

# Add the python directory to sys.path so we can import ctsm modules
_CTSM_PYTHON = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir)
)
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)
