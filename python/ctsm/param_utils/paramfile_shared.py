"""
Functions etc. shared among parameter file utilities
"""

import argparse
import xarray as xr


def open_paramfile(file_in):
    return xr.open_dataset(file_in, decode_timedelta=False)

def paramfile_parser_setup(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input", required=True, help="Input netCDF file")

    # Flags that can be used for the PFT argument
    pft_flags = ["-p", "--pft"]

    return parser, pft_flags
