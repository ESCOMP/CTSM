"""
Functions etc. shared among parameter file utilities
"""

import argparse
import xarray as xr

from ctsm.netcdf_utils import are_xr_dataarrays_identical

PFTNAME_VAR = "pftname"


def are_paramfile_dataarrays_identical(da0: xr.DataArray, da1: xr.DataArray):
    """
    Check whether parameter DataArrays are identical enough, ignoring some metadata
    """
    return are_xr_dataarrays_identical(da0, da1, keys_to_ignore=["source", "original_shape"])


def check_pfts_in_paramfile(selected_pfts, ds):
    """
    Check that the given PFTs are present in the parameter file.

    Parameters
    ----------
    selected_pfts : list of str
        List of PFT names to check.
    ds : xarray.Dataset
        The parameter file dataset.

    Returns
    -------
    list of str
        List of all PFT names in the file.

    Raises
    ------
    KeyError
        If any selected PFT is not found in the file, or if PFTNAME_VAR is missing.
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
    """
    Get the list of PFT names from the parameter file dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        The parameter file dataset.

    Returns
    -------
    list of str
        List of PFT names.
    """
    pft_names = [pft.decode().strip() for pft in ds[PFTNAME_VAR].values]
    return pft_names


def get_selected_pft_indices(selected_pfts, pft_names):
    """
    Get indices of selected PFTs in the list of all PFT names.

    Parameters
    ----------
    selected_pfts : list of str
        List of PFT names to select.
    pft_names : list of str
        List of all PFT names.

    Returns
    -------
    list of int
        Indices of selected PFTs.
    """
    indices = [i for i, name in enumerate(pft_names) if name in selected_pfts]
    return indices


def open_paramfile(file_in, mask_and_scale=False):
    """
    Open a parameter file as an xarray.Dataset.

    Parameters
    ----------
    file_in : str
        Path to the input netCDF file.
    mask_and_scale : bool, optional
        Whether to apply mask and scale (default: False).

    Returns
    -------
    xarray.Dataset
        The opened dataset.
    """
    return xr.open_dataset(file_in, decode_timedelta=False, mask_and_scale=mask_and_scale)


def paramfile_parser_setup(description):
    """
    Set up an argument parser for parameter file utilities.

    Parameters
    ----------
    description : str
        Description for the argument parser.

    Returns
    -------
    tuple
        (parser, pft_flags) where parser is an ArgumentParser and pft_flags is a list of flags for
        PFT argument.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, help="Input netCDF file")

    # Flags that can be used for the PFT argument
    pft_flags = ["-p", "--pft"]

    return parser, pft_flags
