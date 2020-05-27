"""Utility functions related to getting paths to various important places
"""

from __future__ import print_function

import os
import sys

# ========================================================================
# Constants that may need to be changed if directory structures change
# ========================================================================

# Path to the root directory of CTSM, based on the path of this file
#
# Note: It's important that this NOT end with a trailing slash;
# os.path.normpath guarantees this.
_CTSM_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                           os.pardir,
                                           os.pardir))

# Candidates for the last path components to the CTSM directory within a
# CESM checkout
_CESM_CTSM_PATHS = [os.path.join('components', 'ctsm'),
                    os.path.join('components', 'clm')]

# ========================================================================
# Public functions
# ========================================================================

def path_to_ctsm_root():
    """Returns the path to the root directory of CTSM"""
    return _CTSM_ROOT

def path_to_cime(standalone_only=False):
    """Returns the path to cime, if it can be found

    Raises a RuntimeError if it cannot be found

    We first check in the location where cime should be in a standalone
    checkout. If standalone_only is True, then we ONLY look for cime in
    that location. If standalone_only is False, then we fall back to
    checking where cime should be in a full CESM checkout.
    """
    cime_standalone_path = os.path.join(path_to_ctsm_root(), 'cime')
    if os.path.isdir(cime_standalone_path):
        return cime_standalone_path

    if standalone_only:
        raise RuntimeError("Cannot find cime within standalone CTSM checkout")

    cesm_path = _path_to_cesm_root()
    if cesm_path is None:
        raise RuntimeError("Cannot find cime within standalone CTSM checkout, "
                           "and we don't seem to be within a CESM checkout.")

    cime_in_cesm_path = os.path.join(cesm_path, 'cime')
    if os.path.isdir(cime_in_cesm_path):
        return cime_in_cesm_path

    raise RuntimeError("Cannot find cime within standalone CTSM checkout, "
                       "or within CESM checkout rooted at {}".format(cesm_path))

def prepend_to_python_path(path):
    """Adds the given path to python's sys.path if it isn't already in the path

    The path is added near the beginning, so that it takes precedence over existing
    entries in the path
    """
    if not path in sys.path:
        # Insert at location 1 rather than 0, because 0 is special
        sys.path.insert(1, path)

def add_cime_lib_to_path(standalone_only=False):
    """Adds the CIME python library to the python path, to allow importing
    modules from that library

    Returns the path to the top-level cime directory

    For documentation on standalone_only: See documentation in
    path_to_cime
    """
    cime_path = path_to_cime(standalone_only=standalone_only)
    cime_lib_path = os.path.join(cime_path,
                                 'scripts', 'lib')
    prepend_to_python_path(cime_lib_path)
    return cime_path

# ========================================================================
# Private functions
# ========================================================================

def _path_to_cesm_root():
    """Returns the path to the root directory of CESM, if we appear to
    be inside a CESM checkout. If we don't appear to be inside a CESM
    checkout, then returns None.
    """
    ctsm_root = path_to_ctsm_root()
    for candidate_path in _CESM_CTSM_PATHS:
        if ctsm_root.endswith(candidate_path):
            return os.path.normpath(ctsm_root[:-len(candidate_path)])

    return None
