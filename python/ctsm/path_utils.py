"""Utility functions related to getting paths to various important places"""

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
_CTSM_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
)

# Candidates for the last path components to the CTSM directory within a
# CESM checkout
_CESM_CTSM_PATHS = [
    os.path.join("components", "ctsm"),
    os.path.join("components", "clm"),
]

# ========================================================================
# Public functions
# ========================================================================


def path_to_ctsm_root():
    """Returns the path to the root directory of CTSM"""
    return _CTSM_ROOT


def path_to_top_root():
    """Returns the top-level repository root directory (CESM or CTSM)."""
    return _path_to_cesm_root() or path_to_ctsm_root()


def path_under_ctsm(submodule_name):
    """Returns the path to a submodule directory located inside the CTSM root directory.

    Raises a RuntimeError if it cannot be found within the standalone CTSM checkout.
    """
    submod_path = os.path.join(path_to_ctsm_root(), submodule_name)
    if os.path.isdir(submod_path):
        return submod_path

    raise RuntimeError(f"Cannot find {submodule_name} within standalone CTSM checkout")


def path_under_top_submodule(submodule_name):
    """Returns the path to a top-level submodule directory (e.g. 'cime', 'ccs_config').

    Checks <CTSM_ROOT>/<submodule> first (via path_under_ctsm), then falls back to checking
    <CESM_ROOT>/<submodule> if inside a CESM checkout.
    Raises a RuntimeError if it cannot be found.
    """
    try:
        return path_under_ctsm(submodule_name)
    except RuntimeError as exc:
        cesm_path = _path_to_cesm_root()
        if cesm_path is not None:
            cesm_submod = os.path.join(cesm_path, submodule_name)
            if os.path.isdir(cesm_submod):
                return cesm_submod

            raise RuntimeError(
                f"Cannot find {submodule_name} within standalone CTSM checkout, "
                f"or within CESM checkout rooted at {cesm_path}"
            ) from exc
        raise RuntimeError(
            f"Cannot find {submodule_name} within standalone CTSM checkout, "
            "and we don't seem to be within a CESM checkout."
        ) from exc


def path_to_cime(ctsm_only=False):
    """Returns the path to cime, if it can be found

    Raises a RuntimeError if it cannot be found

    We first check in the location where cime should be in a standalone
    checkout. If ctsm_only is True, then we ONLY look for cime in
    that location. If ctsm_only is False, then we fall back to
    checking where cime should be in a full CESM checkout.
    """
    if ctsm_only:
        return path_under_ctsm("cime")
    return path_under_top_submodule("cime")


def path_to_ccs_config(ctsm_only=False):
    """Returns the path to ccs_config, if it can be found

    Raises a RuntimeError if it cannot be found

    We first check in the location where ccs_config should be in a standalone
    checkout. If ctsm_only is True, then we ONLY look for ccs_config in
    that location. If ctsm_only is False, then we fall back to
    checking where ccs_config should be in a full CESM checkout.
    """
    if ctsm_only:
        return path_under_ctsm("ccs_config")
    return path_under_top_submodule("ccs_config")


def prepend_to_python_path(path):
    """Adds the given path to python's sys.path if it isn't already in the path

    The path is added near the beginning, so that it takes precedence over existing
    entries in the path
    """
    if not path in sys.path:
        # Insert at location 1 rather than 0, because 0 is special
        sys.path.insert(1, path)


def add_cime_lib_to_path(ctsm_only=False):
    """Adds the CIME python library to the python path, to allow importing
    modules from that library

    Returns the path to the top-level cime directory

    For documentation on ctsm_only: See documentation in
    path_to_cime
    """
    cime_path = path_to_cime(ctsm_only=ctsm_only)
    prepend_to_python_path(cime_path)
    cime_lib_path = os.path.join(cime_path, "CIME", "Tools")
    prepend_to_python_path(cime_lib_path)
    return cime_path


def add_ctsm_systests_to_path(ctsm_only=False):
    """Adds the CTSM python SystemTests to the python path, to allow importing
    modules from that library
    """
    cime_path = path_to_cime(ctsm_only=ctsm_only)
    ctsm_systest_dir = os.path.join(cime_path, os.pardir, "cime_config")
    prepend_to_python_path(ctsm_systest_dir)
    sys.path.insert(1, ctsm_systest_dir)


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
            return os.path.normpath(ctsm_root[: -len(candidate_path)])

    return None
