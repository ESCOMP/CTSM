"""
Functions etc. shared among parameter file utilities
"""

import argparse
import xarray as xr

PFTNAME_VAR = "pftname"


def check_pfts_in_paramfile(selected_pfts, ds):
    """
    Check that given PFTs are in parameter file
    """
    if PFTNAME_VAR not in ds:
        raise KeyError(f"paramfile missing variable: {PFTNAME_VAR}")
    pft_names = get_pft_names(ds)
    pfts_not_in_file = []
    for pft in selected_pfts:
        if pft not in pft_names:
            pfts_not_in_file += [pft]
    if pfts_not_in_file:
        raise KeyError(f"PFT(s) not found in parameter file: {', '.join(pfts_not_in_file)}")

    return pft_names


def get_pft_names(ds):
    pft_names = [pft.decode().strip() for pft in ds[PFTNAME_VAR].values]
    return pft_names


def get_selected_pft_indices(selected_pfts, pft_names):
    indices = [i for i, name in enumerate(pft_names) if name in selected_pfts]
    return indices


def open_paramfile(file_in, mask_and_scale=False):
    return xr.open_dataset(file_in, decode_timedelta=False, mask_and_scale=mask_and_scale)


def paramfile_parser_setup(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input", required=True, help="Input netCDF file")

    # Flags that can be used for the PFT argument
    pft_flags = ["-p", "--pft"]

    return parser, pft_flags
