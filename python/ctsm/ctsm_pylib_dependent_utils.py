"""
Utilities that are dependent on non-standard modules (i.e., require ctsm_pylib).
"""

import numpy as np
from ctsm.utils import abort


def import_coord_1d(data_set, coord_name):
    """Import 1-d coordinate variable

    Args:
        data_set (xarray Dataset): Dataset whose coordinate you want to import.
        coord_name (str): Name of coordinate to import

    Returns:
        xarray DataArray: DataArray corresponding to the requested coordinate.
    """
    data_array = data_set[coord_name]
    if len(data_array.dims) != 1:
        abort(
            f"Expected 1 dimension for {coord_name}; "
            + f"found {len(data_array.dims)}: {data_array.dims}"
        )
    return data_array, len(data_array)


def import_coord_2d(data_set, coord_name, var_name):
    """
    Import 2-d latitude or longitude variable from a CESM history file (e.g., name LATIXY
    or LONGXY and return it as a 1-d DataArray that can be used as a coordinate for writing
    CESM input files

    Args:
        data_set (xarray Dataset): Dataset whose coordinate you want to import.
        coord_name (str): Name of coordinate to import
        var_name (str): Name of variable with dimension coord_name

    Returns:
        xarray DataArray: 1-d variable that can be used as a coordinate for writing CESM input files
        int: Length of that variable
    """
    data_array = data_set[var_name]
    this_dim = [x for x in data_array.dims if coord_name in x]
    if len(this_dim) != 1:
        abort(
            f"Expected 1 dimension name containing {coord_name}; "
            + f"found {len(this_dim)}: {this_dim}"
        )
    this_dim = this_dim[0]
    other_dim = [x for x in data_array.dims if coord_name not in x]
    if len(other_dim) != 1:
        abort(
            f"Expected 1 dimension name not containing {coord_name}; "
            + f"found {len(other_dim)}: {other_dim}"
        )
    other_dim = other_dim[0]
    data_array = data_array.astype(np.float32)
    data_array = data_array.isel({other_dim: [0]}).squeeze()
    data_array = data_array.rename({this_dim: coord_name}).rename(coord_name)
    data_array = data_array.assign_coords({coord_name: data_array.values})
    data_array.attrs["long_name"] = "coordinate " + data_array.attrs["long_name"]
    data_array.attrs["units"] = data_array.attrs["units"].replace(" ", "_")
    return data_array, len(data_array)
