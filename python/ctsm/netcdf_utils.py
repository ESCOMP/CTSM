"""
Helper functions for working with netCDF files
"""

import numpy as np
import xarray as xr
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm.utils import are_dicts_identical_nansequal


def get_netcdf_format(file_path):
    """
    Get format of netCDF file
    """
    with Dataset(file_path, "r") as netcdf_file:
        netcdf_format = netcdf_file.data_model
    return netcdf_format


def _is_dataarray_metadata_identical(da0: xr.DataArray, da1: xr.DataArray):
    """
    Check whether two DataArrays have identical-enough metadata
    """

    # Check data type
    if da0.dtype != da1.dtype:
        return False

    # Check encoding
    if not are_dicts_identical_nansequal(
        da0.encoding, da1.encoding, keys_to_ignore=["source", "original_shape"]
    ):
        return False

    # Check attributes
    if not are_dicts_identical_nansequal(da0.attrs, da1.attrs):
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


def are_xr_dataarrays_identical(da0: xr.DataArray, da1: xr.DataArray):
    """
    Comprehensively check whether two DataArrays are identical
    """
    if not _is_dataarray_metadata_identical(da0, da1):
        return False

    if not _is_dataarray_data_identical(da0, da1):
        return False

    # Fallback to however xarray defines equality, in case we missed something above
    return da0.equals(da1)
