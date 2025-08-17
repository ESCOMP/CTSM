"""
Helper functions for working with netCDF files
"""

import numpy as np
import xarray as xr
from netCDF4 import Dataset  # pylint: disable=no-name-in-module


def _is_dtype_nan_capable(ndarray: np.ndarray):
    """
    Given a numpy array, return True if it's capable of taking a NaN
    """
    try:
        np.isnan(ndarray)
        return True
    except TypeError:
        return False


def _are_dicts_identical_nansequal(dict0: dict, dict1: dict, keys_to_ignore=None):
    """
    Compare two dictionaries, considering NaNs to be equal. Don't be strict here about types; if
    they can be coerced to comparable types and then they match, return True.
    """
    # pylint: disable=too-many-return-statements

    if keys_to_ignore is None:
        keys_to_ignore = []
    keys_to_ignore = np.array(keys_to_ignore)

    if len(dict0) != len(dict1):
        return False
    for key, value0 in dict0.items():
        if key in keys_to_ignore:
            continue
        if key not in dict1:
            return False
        value1 = dict1[key]

        # Coerce to numpy arrays to simplify comparison code
        value0 = np.array(value0)
        value1 = np.array(value1)

        # Compare, only asking to check equal NaNs if both are capable of taking NaN values
        both_are_nan_capable = _is_dtype_nan_capable(value0) and _is_dtype_nan_capable(value1)
        if not np.array_equal(value0, value1, equal_nan=both_are_nan_capable):
            return False

    return True


def get_netcdf_format(file_path):
    """
    Get format of netCDF file
    """
    with Dataset(file_path, "r") as netcdf_file:
        netcdf_format = netcdf_file.data_model
    return netcdf_format


def _is_dataarray_metadata_identical(da0: xr.DataArray, da1: xr.DataArray, keys_to_ignore=None):
    """
    Check whether two DataArrays have identical-enough metadata
    """

    # Check data type
    if da0.dtype != da1.dtype:
        return False

    # Check encoding
    if not _are_dicts_identical_nansequal(
        da0.encoding, da1.encoding, keys_to_ignore=keys_to_ignore
    ):
        return False

    # Check attributes
    if not _are_dicts_identical_nansequal(da0.attrs, da1.attrs):
        return False

    # Check name
    if da0.name != da1.name:
        return False

    # Check dims
    if da0.dims != da1.dims:
        return False

    return True


def _is_dataarray_data_identical(da0: xr.DataArray, da1: xr.DataArray):
    """
    Check whether two DataArrays have identical data
    """
    # pylint: disable=too-many-return-statements

    # Check sizes
    if da0.sizes != da1.sizes:
        return False

    # Check coordinates
    if bool(da0.coords) or bool(da1.coords):
        if not bool(da0.coords) or not bool(da1.coords):
            return False
        if not da0.coords.equals(da1.coords):
            return False

    # Check values ("The array's data converted to numpy.ndarray")
    if not np.array_equal(da0.values, da1.values):
        # Try-except to avoid TypeError from putting NaN-incapable dtypes through
        # np.array_equal(..., equal_nan=True)
        try:
            if not np.array_equal(da0.values, da1.values, equal_nan=True):
                return False
        except TypeError:
            return False

    # Check data ("The DataArray's data as an array. The underlying array type (e.g. dask, sparse,
    # pint) is preserved.")
    da0_data_type = type(da0.data)
    if not isinstance(da1.data, da0_data_type):
        return False
    if not isinstance(da0.data, np.ndarray):
        raise NotImplementedError(f"Add support for comparing two objects of type {da0_data_type}")

    return True


def are_xr_dataarrays_identical(da0: xr.DataArray, da1: xr.DataArray, keys_to_ignore=None):
    """
    Comprehensively check whether two DataArrays are identical
    """
    if not _is_dataarray_metadata_identical(da0, da1, keys_to_ignore=keys_to_ignore):
        return False

    if not _is_dataarray_data_identical(da0, da1):
        return False

    # Fallback to however xarray defines equality, in case we missed something above
    return da0.equals(da1)
