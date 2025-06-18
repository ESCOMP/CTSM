"""
Things shared between plumber2 scripts
"""

import os
import pandas as pd

PLUMBER2_SITES_CSV = os.path.realpath(
    os.path.join(
        os.path.dirname(__file__),
        os.pardir,
        os.pardir,
        os.pardir,
        "tools",
        "site_and_regional",
        "PLUMBER2_sites.csv",
    )
)


def read_plumber2_sites_csv():
    """
    Read PLUMBER2_sites.csv using pandas
    """
    return pd.read_csv(PLUMBER2_SITES_CSV, skiprows=5)
