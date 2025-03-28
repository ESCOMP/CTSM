"""
Import a dataset that can be spread over multiple files, only including specified variables
and/or vegetation types and/or timesteps, concatenating by time.

- DOES actually read the dataset into memory, but only AFTER dropping unwanted variables and/or
    vegetation types.
"""

import re
import warnings
from importlib.util import find_spec
import numpy as np
import xarray as xr
from ctsm.utils import is_instantaneous
import ctsm.crop_calendars.cropcal_utils as utils
from ctsm.crop_calendars.xr_flexsel import xr_flexsel


def compute_derived_vars(ds_in, var):
    """
    Compute derived variables
    """
    if (
        var == "HYEARS"
        and "HDATES" in ds_in
        and ds_in.HDATES.dims == ("time", "mxharvests", "patch")
    ):
        year_adj = 1 if is_instantaneous(ds_in["time"]) else 0
        year_list = np.array([np.float32(x.year - year_adj) for x in ds_in.time.values])
        hyears = ds_in["HDATES"].copy()
        hyears.values = np.tile(
            np.expand_dims(year_list, (1, 2)),
            (1, ds_in.dims["mxharvests"], ds_in.dims["patch"]),
        )
        with np.errstate(invalid="ignore"):
            is_le_zero = ~np.isnan(ds_in.HDATES.values) & (ds_in.HDATES.values <= 0)
        hyears.values[is_le_zero] = ds_in.HDATES.values[is_le_zero]
        hyears.values[np.isnan(ds_in.HDATES.values)] = np.nan
        hyears.attrs["long_name"] = "DERIVED: actual crop harvest years"
        hyears.attrs["units"] = "year"
        ds_in["HYEARS"] = hyears
    else:
        raise RuntimeError(f"Unable to compute derived variable {var}")
    return ds_in


def manual_mfdataset(filelist, my_vars, my_vegtypes, time_slice):
    """
    Opening a list of files with Xarray's open_mfdataset requires dask. This function is a
    workaround for Python environments that don't have dask.
    """
    ds_out = None
    for filename in filelist:
        ds_in = xr.open_dataset(filename)
        ds_in = mfdataset_preproc(ds_in, my_vars, my_vegtypes, time_slice)
        if ds_out is None:
            ds_out = ds_in
        else:
            ds_out = xr.concat(
                [ds_out, ds_in],
                data_vars="minimal",
                compat="override",
                coords="all",
                dim="time",
            )
    return ds_out


def mfdataset_preproc(ds_in, vars_to_import, vegtypes_to_import, time_slice):
    """
    Function to drop unwanted variables in preprocessing of open_mfdataset().

    - Makes sure to NOT drop any unspecified variables that will be useful in gridding.
    - Also adds vegetation type info in the form of a DataArray of strings.
    - Also renames "pft" dimension (and all like-named variables, e.g., pft1d_itype_veg_str) to be
      named like "patch". This can later be reversed, for compatibility with other code, using
      patch2pft().
    """
    # Rename "pft" dimension and variables to "patch", if needed
    if "pft" in ds_in.dims:
        pattern = re.compile("pft.*1d")
        matches = [x for x in list(ds_in.keys()) if pattern.search(x) is not None]
        pft2patch_dict = {"pft": "patch"}
        for match in matches:
            pft2patch_dict[match] = match.replace("pft", "patch").replace("patchs", "patches")
        ds_in = ds_in.rename(pft2patch_dict)

    derived_vars = []
    if vars_to_import is not None:
        # Split vars_to_import into variables that are vs. aren't already in ds
        derived_vars = [v for v in vars_to_import if v not in ds_in]
        present_vars = [v for v in vars_to_import if v in ds_in]
        vars_to_import = present_vars

        # Get list of dimensions present in variables in vars_to_import.
        dim_list = []
        for var in vars_to_import:
            # list(set(x)) returns a list of the unique items in x
            dim_list = list(set(dim_list + list(ds_in.variables[var].dims)))

        # Get any _1d variables that are associated with those dimensions. These will be useful in
        # gridding. Also, if any dimension is "pft", set up to rename it and all like-named
        # variables to "patch"
        oned_vars = []
        for dim in dim_list:
            pattern = re.compile(f"{dim}.*1d")
            matches = [x for x in list(ds_in.keys()) if pattern.search(x) is not None]
            oned_vars = list(set(oned_vars + matches))

        # Add dimensions and _1d variables to vars_to_import
        vars_to_import = list(set(vars_to_import + list(ds_in.dims) + oned_vars))

        # Add any _bounds variables
        bounds_vars = []
        for var in vars_to_import:
            bounds_var = var + "_bounds"
            if bounds_var in ds_in:
                bounds_vars = bounds_vars + [bounds_var]
        vars_to_import = vars_to_import + bounds_vars

        # Get list of variables to drop
        varlist = list(ds_in.variables)
        vars_to_drop = list(np.setdiff1d(varlist, vars_to_import))

        # Drop them
        ds_in = ds_in.drop_vars(vars_to_drop)

    # Add vegetation type info
    if "patches1d_itype_veg" in list(ds_in):
        this_pftlist = utils.define_pftlist()
        utils.get_patch_ivts(
            ds_in, this_pftlist
        )  # Includes check of whether vegtype changes over time anywhere
        vegtype_da = utils.get_vegtype_str_da(this_pftlist)
        patches1d_itype_veg_str = vegtype_da.values[
            ds_in.isel(time=0).patches1d_itype_veg.values.astype(int)
        ]
        npatch = len(patches1d_itype_veg_str)
        patches1d_itype_veg_str = xr.DataArray(
            patches1d_itype_veg_str,
            coords={"patch": np.arange(0, npatch)},
            dims=["patch"],
            name="patches1d_itype_veg_str",
        )
        ds_in = xr.merge([ds_in, vegtype_da, patches1d_itype_veg_str])

    # Restrict to veg. types of interest, if any
    if vegtypes_to_import is not None:
        ds_in = xr_flexsel(ds_in, vegtype=vegtypes_to_import)

    # Restrict to time slice, if any
    if time_slice:
        ds_in = utils.safer_timeslice(ds_in, time_slice)

    # Finish import
    ds_in = xr.decode_cf(ds_in, decode_times=True)

    # Compute derived variables
    for var in derived_vars:
        ds_in = compute_derived_vars(ds_in, var)

    return ds_in


def process_inputs(filelist, my_vars, my_vegtypes, my_vars_missing_ok):
    """
    Process inputs to import_ds()
    """
    if my_vars_missing_ok is None:
        my_vars_missing_ok = []
    # Convert my_vegtypes here, if needed, to avoid repeating the process each time you read a file
    # in xr.open_mfdataset().
    if my_vegtypes is not None:
        if not isinstance(my_vegtypes, list):
            my_vegtypes = [my_vegtypes]
        if isinstance(my_vegtypes[0], str):
            my_vegtypes = utils.vegtype_str2int(my_vegtypes)

    # Same for these variables.
    if my_vars is not None:
        if not isinstance(my_vars, list):
            my_vars = [my_vars]
    if my_vars_missing_ok:
        if not isinstance(my_vars_missing_ok, list):
            my_vars_missing_ok = [my_vars_missing_ok]

    # Make sure lists are actually lists
    if not isinstance(filelist, list):
        filelist = [filelist]
    if not isinstance(my_vars_missing_ok, list):
        my_vars_missing_ok = [my_vars_missing_ok]
    return filelist, my_vars, my_vegtypes, my_vars_missing_ok


def import_ds(
    filelist,
    *,
    my_vars=None,
    my_vegtypes=None,
    time_slice=None,
    my_vars_missing_ok=None,
    rename_lsmlatlon=False,
    chunks=None,
):
    """
    Import a dataset that can be spread over multiple files, only including specified variables
    and/or vegetation types and/or timesteps, concatenating by time.

    - DOES actually read the dataset into memory, but only AFTER dropping unwanted variables and/or
      vegetation types.
    """
    filelist, my_vars, my_vegtypes, my_vars_missing_ok = process_inputs(
        filelist, my_vars, my_vegtypes, my_vars_missing_ok
    )

    # Remove files from list if they don't contain requested timesteps.
    # time_slice should be in the format slice(start,end[,step]). start or end can be None to be
    # unbounded on one side. Note that the standard slice() documentation suggests that only
    # elements through end-1 will be selected, but that seems not to be the case in the xarray
    # implementation.
    if time_slice:
        new_filelist = []
        for file in sorted(filelist):
            filetime = xr.open_dataset(file).time
            filetime_sel = utils.safer_timeslice(filetime, time_slice)
            include_this_file = filetime_sel.size
            if include_this_file:
                new_filelist.append(file)

            # If you found some matching files, but then you find one that doesn't, stop going
            # through the list.
            elif new_filelist:
                break
        if not new_filelist:
            raise RuntimeError(f"No files found in time_slice {time_slice}")
        filelist = new_filelist

    # The xarray open_mfdataset() "preprocess" argument requires a function that takes exactly one
    # variable (an xarray.Dataset object). Wrapping mfdataset_preproc() in this lambda function
    # allows this. Could also just allow mfdataset_preproc() to access my_vars and my_vegtypes
    # directly, but that's bad practice as it could lead to scoping issues.
    # pylint: disable=unnecessary-lambda-assignment
    mfdataset_preproc_closure = lambda ds: mfdataset_preproc(ds, my_vars, my_vegtypes, time_slice)

    # Import
    if isinstance(filelist, list) and len(filelist) == 1:
        filelist = filelist[0]
    if isinstance(filelist, list):
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=DeprecationWarning)
            dask_unavailable = find_spec("dask") is None
        if dask_unavailable:
            this_ds = manual_mfdataset(filelist, my_vars, my_vegtypes, time_slice)
        else:
            this_ds = xr.open_mfdataset(
                sorted(filelist),
                data_vars="minimal",
                preprocess=mfdataset_preproc_closure,
                compat="override",
                coords="all",
                concat_dim="time",
                combine="nested",
                chunks=chunks,
            )
    elif isinstance(filelist, str):
        this_ds = xr.open_dataset(filelist, chunks=chunks)
        this_ds = mfdataset_preproc(this_ds, my_vars, my_vegtypes, time_slice)
        this_ds = this_ds.compute()

    # Warn and/or error about variables that couldn't be imported or derived
    if my_vars:
        missing_vars = [v for v in my_vars if v not in this_ds]
        ok_missing_vars = [v for v in missing_vars if v in my_vars_missing_ok]
        bad_missing_vars = [v for v in missing_vars if v not in my_vars_missing_ok]
        if ok_missing_vars:
            print(
                "Could not import some variables; either not present or not deriveable:"
                f" {ok_missing_vars}"
            )
        if bad_missing_vars:
            raise RuntimeError(
                "Could not import some variables; either not present or not deriveable:"
                f" {bad_missing_vars}"
            )

    if rename_lsmlatlon:
        if "lsmlat" in this_ds.dims:
            this_ds = this_ds.rename({"lsmlat": "lat"})
        if "lsmlon" in this_ds.dims:
            this_ds = this_ds.rename({"lsmlon": "lon"})

    return this_ds
