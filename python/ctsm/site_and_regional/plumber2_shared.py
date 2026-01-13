"""
Things shared between plumber2 scripts
"""

import os
import pandas as pd
from ctsm.path_utils import path_to_ctsm_root

PLUMBER2_SITES_CSV = os.path.join(
    path_to_ctsm_root(),
    "tools",
    "site_and_regional",
    "PLUMBER2_sites.csv",
)


def read_plumber2_sites_csv(file=PLUMBER2_SITES_CSV):
    """
    Read PLUMBER2_sites.csv using pandas
    """
    return pd.read_csv(file, skiprows=5)
