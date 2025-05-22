"""General-purpose git utility functions"""

import logging
import subprocess

from ctsm.path_utils import path_to_ctsm_root

logger = logging.getLogger(__name__)


def get_ctsm_git_short_hash():
    """
    Returns Git short SHA for the CTSM repository.

    Args:

    Raises:

    Returns:
        sha (str) : git short hash for ctsm repository
    """
    sha = (
        subprocess.check_output(["git", "-C", path_to_ctsm_root(), "rev-parse", "--short", "HEAD"])
        .strip()
        .decode()
    )
    return sha


def get_ctsm_git_long_hash():
    """
    Returns Git long SHA for the CTSM repository.

    Args:

    Raises:

    Returns:
        sha (str) : git long hash for ctsm repository
    """
    sha = (
        subprocess.check_output(["git", "-C", path_to_ctsm_root(), "rev-parse", "HEAD"])
        .strip()
        .decode()
    )
    return sha


def get_ctsm_git_describe():
    """
    Function for giving the recent tag of the CTSM repository

    Args:

    Raises:

    Returns:
        label (str) : ouput of running 'git describe' for the CTSM repository
    """
    label = subprocess.check_output(["git", "-C", path_to_ctsm_root(), "describe"]).strip().decode()
    return label
