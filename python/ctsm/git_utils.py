"""General-purpose git utility functions"""

import logging
import subprocess

logger = logging.getLogger(__name__)


def get_git_short_hash():
    """
    Returns Git short SHA for the currect directory.
    """
    try:

        # os.abspath(__file__)
        sha = (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
            .strip()
            .decode()
        )
    except subprocess.CalledProcessError:
        sha = "NOT-A-GIT-REPOSITORY"
    return sha


def get_git_long_hash():
    """
    Returns Git long SHA for the currect directory.
    """
    try:

        # os.abspath(__file__)
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode()
    except subprocess.CalledProcessError:
        sha = "NOT-A-GIT-REPOSITORY"
    return sha

def get_git_describe():
    """
    Returns git describe output
    """
    label = subprocess.check_output(["git", "describe"]).strip()
    return label.decode()
