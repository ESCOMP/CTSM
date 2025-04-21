"""
Make a geographically gridded DataArray (with dimensions time, vegetation type [as string], lat,
lon) of one variable within a Dataset.

- Optional keyword arguments will be passed to xr_flexsel() to select single steps or slices
    along the specified ax(ie)s.
- fill_value: Default None means grid will be filled with NaN, unless the variable in question
    already has a _FillValue, in which case that will be used.
"""

import numpy as np
import xarray as xr
from ctsm.crop_calendars.xr_flexsel import xr_flexsel

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments


def get_thisvar_da(var, this_ds):
    """
    Return a DataArray, with defined coordinates, for a given variable in a dataset.
    """
    # Make DataArray for this variable
    thisvar_da = np.array(this_ds.variables[var])
    these_dims = this_ds.variables[var].dims
    thisvar_da = xr.DataArray(thisvar_da, dims=these_dims)

    # Define coordinates of this variable's DataArray
    dims_dict = {}
    for dim in these_dims:
        dims_dict[dim] = this_ds[dim]
    thisvar_da = thisvar_da.assign_coords(dims_dict)
    thisvar_da.attrs = this_ds[var].attrs

    return thisvar_da


def convert_to_da(this_ds, var, fill_value, thisvar_da, new_dims, thisvar_gridded):
    """
    Convert Numpy array to DataArray with coordinates, attributes and name
    """
    thisvar_gridded = xr.DataArray(thisvar_gridded, dims=tuple(new_dims), attrs=thisvar_da.attrs)
    for dim in new_dims:
        if dim == "ivt_str":
            values = this_ds.vegtype_str.values
        elif dim in thisvar_da.coords:
            values = thisvar_da[dim]
        else:
            values = this_ds[dim].values
        thisvar_gridded = thisvar_gridded.assign_coords({dim: values})
    thisvar_gridded.name = var

    # Add FillValue attribute
    if fill_value:
        thisvar_gridded.attrs["_FillValue"] = fill_value
    return thisvar_gridded


def grid_the_data(thisvar_da, vt_da, ixy_da, jxy_da, new_dims, thisvar_gridded):
    """
    Fill lat-lon array with previously-ungridded data
    """
    fill_indices = []
    for dim in new_dims:
        if dim == "lat":
            fill_indices.append(jxy_da.values.astype(int) - 1)
        elif dim == "lon":
            fill_indices.append(ixy_da.values.astype(int) - 1)
        elif dim == "ivt_str":
            fill_indices.append(vt_da)
        elif not fill_indices:
            # I.e., if fill_indices is empty. Could also do "elif len(fill_indices)==0".
            fill_indices.append(Ellipsis)
    try:
        thisvar_gridded[tuple(fill_indices[: len(fill_indices)])] = thisvar_da.values
    except:  # pylint: disable=bare-except
        thisvar_gridded[tuple(fill_indices[: len(fill_indices)])] = thisvar_da.values.transpose()
    if not np.any(np.bitwise_not(np.isnan(thisvar_gridded))):
        if np.all(np.isnan(thisvar_da.values)):
            print("Warning: This DataArray (and thus map) is all NaN")
        else:
            raise RuntimeError("thisvar_gridded was not filled!")


def create_filled_array(this_ds, fill_value, thisvar_da, new_dims):
    """
    Create a Numpy array to be filled with gridded data
    """
    dim_size_list = []
    for dim in new_dims:
        if dim == "ivt_str":
            dim_size = this_ds.sizes["ivt"]
        elif dim in thisvar_da.coords:
            dim_size = thisvar_da.sizes[dim]
        else:
            dim_size = this_ds.sizes[dim]
        dim_size_list = dim_size_list + [dim_size]
    thisvar_gridded = np.empty(dim_size_list)
    if fill_value:
        thisvar_gridded[:] = fill_value
    else:
        thisvar_gridded[:] = np.NaN
    return thisvar_gridded


def get_ixy_jxy_das(this_ds, var):
    """
    Get DataArrays needed for gridding
    """
    thisvar_da = get_thisvar_da(var, this_ds)
    vt_da = None
    if "patch" in thisvar_da.dims:
        spatial_unit = "patch"
        xy_1d_prefix = "patches"
        if "patches1d_itype_veg" in this_ds:
            vt_da = get_thisvar_da("patches1d_itype_veg", this_ds)
    elif "gridcell" in thisvar_da.dims:
        spatial_unit = "gridcell"
        xy_1d_prefix = "grid"
    else:
        raise RuntimeError(
            f"What variables to use for _ixy and _jxy of variable with dims {thisvar_da.dims}?"
        )
    ixy_da = get_thisvar_da(xy_1d_prefix + "1d_ixy", this_ds)
    jxy_da = get_thisvar_da(xy_1d_prefix + "1d_jxy", this_ds)
    return thisvar_da, vt_da, spatial_unit, ixy_da, jxy_da


def get_new_dim_list(this_ds, thisvar_da, spatial_unit):
    """
    Get new dimension list
    """
    new_dims = list(thisvar_da.dims)
    ### Remove "[spatial_unit]".
    if spatial_unit in new_dims:
        new_dims.remove(spatial_unit)
    #  Add "ivt_str" (vegetation type, as string). This needs to go at the end, to avoid a possible
    # situation where you wind up with multiple Ellipsis members of fill_indices.
    if "ivt" in this_ds and spatial_unit == "patch":
        new_dims.append("ivt_str")
    ### Add lat and lon to end of list
    new_dims = new_dims + ["lat", "lon"]
    return new_dims


def grid_one_variable(this_ds, var, fill_value=None, **kwargs):
    """
    Make a geographically gridded DataArray (with dimensions time, vegetation type [as string], lat,
    lon) of one variable within a Dataset.

    - Optional keyword arguments will be passed to xr_flexsel() to select single steps or slices
      along the specified ax(ie)s.
    - fill_value: Default None means grid will be filled with NaN, unless the variable in question
      already has a _FillValue, in which case that will be used.
    """
    # Get this Dataset's values for selection(s), if provided
    this_ds = xr_flexsel(this_ds, **kwargs)

    # Get DataArrays needed for gridding
    thisvar_da, vt_da, spatial_unit, ixy_da, jxy_da = get_ixy_jxy_das(this_ds, var)

    if not fill_value and "_FillValue" in thisvar_da.attrs:
        fill_value = thisvar_da.attrs["_FillValue"]

    # Renumber vt_da to work as indices on new ivt dimension, if needed.
    ### Ensures that the unique set of vt_da values begins with 1 and
    ### contains no missing steps.
    if "ivt" in this_ds and vt_da is not None:
        vt_da.values = np.array([np.where(this_ds.ivt.values == x)[0][0] for x in vt_da.values])

    # Get new dimension list
    new_dims = get_new_dim_list(this_ds, thisvar_da, spatial_unit)

    # Create a Numpy array to be filled with gridded data
    thisvar_gridded = create_filled_array(this_ds, fill_value, thisvar_da, new_dims)

    # Fill lat-lon array with previously-ungridded data
    grid_the_data(thisvar_da, vt_da, ixy_da, jxy_da, new_dims, thisvar_gridded)

    # Convert Numpy array to DataArray with coordinates, attributes and name
    thisvar_gridded = convert_to_da(this_ds, var, fill_value, thisvar_da, new_dims, thisvar_gridded)

    return thisvar_gridded
