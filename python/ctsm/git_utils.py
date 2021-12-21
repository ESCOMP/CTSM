"""General-purpose git utility functions"""

import logging
import subprocess
from ctsm.path_utils import path_to_ctsm_root

logger = logging.getLogger(__name__)


def get_git_short_hash():
    """
    Returns Git short SHA for the currect directory.
    """
    sha = subprocess.check_output(['git', '-C', path_to_ctsm_root(),
                'rev-parse', '--short', 'HEAD']).strip().decode()
    return sha


def get_git_long_hash():
    """
    Returns Git long SHA for the currect directory.
    """
    sha = subprocess.check_output(["git", '-C', path_to_ctsm_root(), "rev-parse", "HEAD"]).strip().decode()
    return sha

def get_git_describe():
    """
    Returns git describe output
    """
    label = subprocess.check_output(["git", "describe"]).strip().decode()
    return label
