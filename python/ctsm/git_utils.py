"""General-purpose git utility functions"""

import logging
import subprocess

logger = logging.getLogger(__name__)


def tag_describe():
    """
    Function for giving the recent tag of the git repo

    Args:

    Raises:

    Returns:
        label.decode (str) : ouput of running 'git describe' in shell
    """
    label = subprocess.check_output(["git", "describe"]).strip()
    return label.decode()
