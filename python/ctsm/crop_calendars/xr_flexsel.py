"""
Flexibly subset time(s) and/or vegetation type(s) from an xarray Dataset or DataArray.
"""

import re
import numpy as np
import xarray as xr

from ctsm.crop_calendars.cropcal_utils import vegtype_str2int, is_each_vegtype

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments


def xr_flexsel(xr_object, patches1d_itype_veg=None, warn_about_seltype_interp=True, **kwargs):
    """
    Flexibly subset time(s) and/or vegetation type(s) from an xarray Dataset or DataArray.

    - Keyword arguments like dimension=selection.
    - Selections can be individual values or slice()s.
    - Optimize memory usage by beginning keyword argument list with the selections that will result
      in the largest reduction of object size.
    - Use dimension "vegtype" to extract patches of designated vegetation type (can be string or
      integer).
    - Can also do dimension=function---e.g., time=np.mean will take the mean over the time
      dimension.
    """
    # Setup
    havewarned = False
    delimiter = "__"

    for key, selection in kwargs.items():
        if callable(selection):
            xr_object = handle_callable(xr_object, key, selection)

        elif key == "vegtype":
            xr_object = handle_vegtype(xr_object, patches1d_itype_veg, selection)

        else:
            # Parse selection type, if provided
            if delimiter in key:
                key, selection_type = key.split(delimiter)

            # Check type of selection
            else:
                is_inefficient = False
                if isinstance(selection, slice):
                    this_type = set_type_from_slice(selection)
                elif isinstance(selection, np.ndarray):
                    selection, is_inefficient, this_type = set_type_from_ndarray(selection)
                else:
                    this_type = type(selection)

                warn_about_this_seltype_interp = warn_about_seltype_interp
                if this_type == list and isinstance(selection[0], str):
                    selection_type = "values"
                    warn_about_this_seltype_interp = False
                elif this_type == int:
                    selection_type = "indices"
                else:
                    selection_type = "values"

                if warn_about_this_seltype_interp:
                    do_warn_about_seltype_interp(
                        havewarned, delimiter, key, selection_type, is_inefficient, this_type
                    )

            # Trim along relevant 1d axes
            if isinstance(xr_object, xr.Dataset) and key in ["lat", "lon"]:
                xr_object = trim_along_relevant_1d_axes(xr_object, selection, selection_type, key)

            # Perform selection
            xr_object = perform_selection(xr_object, key, selection, selection_type)

    return xr_object


def perform_selection(xr_object, key, selection, selection_type):
    """
    Perform selection
    """
    if selection_type == "indices":
        # Have to select like this instead of with index directly because otherwise assign_coords()
        # will throw an error. Not sure why.
        if isinstance(selection, int):
            # Single integer? Turn it into a slice.
            selection = slice(selection, selection + 1)
        elif (
            isinstance(selection, np.ndarray)
            and not selection.dtype.kind in np.typecodes["AllInteger"]
        ):
            selection = selection.astype(int)
        xr_object = xr_object.isel({key: selection})
    elif selection_type == "values":
        xr_object = xr_object.sel({key: selection})
    else:
        raise TypeError(f"selection_type {selection_type} not recognized")
    return xr_object


def trim_along_relevant_1d_axes(xr_object, selection, selection_type, key):
    """
    Trim along relevant 1d axes
    """
    if selection_type == "indices":
        incl_coords = xr_object[key].values[selection]
    elif selection_type == "values":
        if isinstance(selection, slice):
            incl_coords = xr_object.sel({key: selection}, drop=False)[key].values
        else:
            incl_coords = selection
    else:
        raise TypeError(f"selection_type {selection_type} not recognized")
    if key == "lat":
        this_xy = "jxy"
    elif key == "lon":
        this_xy = "ixy"
    else:
        raise KeyError(
            f"Key '{key}' not recognized: What 1d_ suffix should I use for variable name?"
        )
    pattern = re.compile(f"1d_{this_xy}")
    matches = [x for x in list(xr_object.keys()) if pattern.search(x) is not None]
    for var in matches:
        if len(xr_object[var].dims) != 1:
            raise RuntimeError(
                f"Expected {var} to have 1 dimension, but it has"
                f" {len(xr_object[var].dims)}: {xr_object[var].dims}"
            )
        dim = xr_object[var].dims[0]
        # print(f"Variable {var} has dimension {dim}")
        coords = xr_object[key].values[xr_object[var].values.astype(int) - 1]
        # print(f"{dim} size before: {xr_object.sizes[dim]}")
        ok_ind = []
        new_1d_this_xy = []
        for i, member in enumerate(coords):
            if member in incl_coords:
                ok_ind = ok_ind + [i]
                new_1d_this_xy = new_1d_this_xy + [(incl_coords == member).nonzero()[0] + 1]
        xr_object = xr_object.isel({dim: ok_ind})
        new_1d_this_xy = np.array(new_1d_this_xy).squeeze()
        xr_object[var].values = new_1d_this_xy
        # print(f"{dim} size after: {xr_object.sizes[dim]}")
    return xr_object


def do_warn_about_seltype_interp(
    havewarned, delimiter, key, selection_type, is_inefficient, this_type
):
    """
    Suggest suppressing selection type interpretation warnings
    """
    if not havewarned:
        print(
            "xr_flexsel(): Suppress all 'selection type interpretation' messages by specifying"
            + "warn_about_seltype_interp=False"
        )
        havewarned = True
    if is_inefficient:
        extra = " This will also improve efficiency for large selections."
    else:
        extra = ""
    print(
        f"xr_flexsel(): Selecting {key} as {selection_type} because selection was"
        f" interpreted as {this_type}. If not correct, specify selection type"
        " ('indices' or 'values') in keyword like"
        f" '{key}{delimiter}SELECTIONTYPE=...' instead of '{key}=...'.{extra}"
    )


def set_type_from_ndarray(selection):
    """
    Sets selection type if given a Numpy array
    """
    if selection.dtype.kind in np.typecodes["AllInteger"]:
        this_type = int
    else:
        is_inefficient = True
        this_type = None
        for member in selection:
            if member < 0 or member % 1 > 0:
                if isinstance(member, int):
                    this_type = "values"
                else:
                    this_type = type(member)
                break
        if this_type is None:
            this_type = int
            selection = selection.astype(int)
    return selection, is_inefficient, this_type


def set_type_from_slice(selection):
    """
    Sets selection type if given a slice
    """
    slice_members = []
    if selection == slice(0):
        raise ValueError("slice(0) will be empty")
    if selection.start is not None:
        slice_members = slice_members + [selection.start]
    if selection.stop is not None:
        slice_members = slice_members + [selection.stop]
    if selection.step is not None:
        slice_members = slice_members + [selection.step]
    if not slice_members:
        raise TypeError("slice is all None?")
    this_type = int
    for member in slice_members:
        if member < 0 or not isinstance(member, int):
            this_type = "values"
            break
    return this_type


def handle_vegtype(xr_object, patches1d_itype_veg, selection):
    """
    Handle selection "vegtype
    """
    # Convert to list, if needed
    if not isinstance(selection, list):
        selection = [selection]

    # Convert to indices, if needed
    if isinstance(selection[0], str):
        selection = vegtype_str2int(selection)

    # Get list of boolean(s)
    if isinstance(selection[0], int):
        if isinstance(patches1d_itype_veg, type(None)):
            patches1d_itype_veg = xr_object.patches1d_itype_veg.values
        elif isinstance(patches1d_itype_veg, xr.core.dataarray.DataArray):
            patches1d_itype_veg = patches1d_itype_veg.values
        is_vegtype = is_each_vegtype(patches1d_itype_veg, selection, "ok_exact")
    elif isinstance(selection[0], bool):
        if len(selection) != len(xr_object.patch):
            raise ValueError(
                "If providing boolean 'vegtype' argument to xr_flexsel(), it must be the"
                f" same length as xr_object.patch ({len(selection)} vs."
                f" {len(xr_object.patch)})"
            )
        is_vegtype = selection
    else:
        raise TypeError(f"Not sure how to handle 'vegtype' of type {type(selection[0])}")
    xr_object = xr_object.isel(patch=[i for i, x in enumerate(is_vegtype) if x])
    if "ivt" in xr_object:
        xr_object = xr_object.isel(ivt=is_each_vegtype(xr_object.ivt.values, selection, "ok_exact"))

    return xr_object


def handle_callable(xr_object, key, selection):
    """
    Handle selection that's a callable
    """
    # It would have been really nice to do selection(xr_object, axis=key), but numpy methods and
    # xarray methods disagree on "axis" vs. "dimension." So instead, just do this manually.
    if selection == np.mean:  # pylint: disable=comparison-with-callable
        try:
            xr_object = xr_object.mean(dim=key)
        except:  # pylint: disable=raise-missing-from
            raise ValueError(
                f"Failed to take mean of dimension {key}. Try doing so outside of xr_flexsel()."
            )
    else:
        raise ValueError(f"xr_flexsel() doesn't recognize function {selection}")
    return xr_object
