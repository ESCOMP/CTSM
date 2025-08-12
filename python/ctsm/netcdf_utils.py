"""
Helper functions for working with netCDF files
"""

from netCDF4 import Dataset  # pylint: disable=no-name-in-module


def get_netcdf_format(file_path):
    """
    Get format of netCDF file
    """
    with Dataset(file_path, "r") as netcdf_file:
        netcdf_format = netcdf_file.data_model
    return netcdf_format
