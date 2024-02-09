"""
utility functions
copied from klindsay, https://github.com/klindsay28/CESM2_coup_carb_cycle_JAMES/blob/master/utils.py
"""

import re
import warnings
import importlib

import numpy as np
import xarray as xr
from ctsm.crop_calendars.xr_flexsel import xr_flexsel


def define_pftlist():
    """
    Return list of PFTs used in CLM
    """
    pftlist = [
        "not_vegetated",
        "needleleaf_evergreen_temperate_tree",
        "needleleaf_evergreen_boreal_tree",
        "needleleaf_deciduous_boreal_tree",
        "broadleaf_evergreen_tropical_tree",
        "broadleaf_evergreen_temperate_tree",
        "broadleaf_deciduous_tropical_tree",
        "broadleaf_deciduous_temperate_tree",
        "broadleaf_deciduous_boreal_tree",
        "broadleaf_evergreen_shrub",
        "broadleaf_deciduous_temperate_shrub",
        "broadleaf_deciduous_boreal_shrub",
        "c3_arctic_grass",
        "c3_non-arctic_grass",
        "c4_grass",
        "unmanaged_c3_crop",
        "unmanaged_c3_irrigated",
        "temperate_corn",
        "irrigated_temperate_corn",
        "spring_wheat",
        "irrigated_spring_wheat",
        "winter_wheat",
        "irrigated_winter_wheat",
        "soybean",
        "irrigated_soybean",
        "barley",
        "irrigated_barley",
        "winter_barley",
        "irrigated_winter_barley",
        "rye",
        "irrigated_rye",
        "winter_rye",
        "irrigated_winter_rye",
        "cassava",
        "irrigated_cassava",
        "citrus",
        "irrigated_citrus",
        "cocoa",
        "irrigated_cocoa",
        "coffee",
        "irrigated_coffee",
        "cotton",
        "irrigated_cotton",
        "datepalm",
        "irrigated_datepalm",
        "foddergrass",
        "irrigated_foddergrass",
        "grapes",
        "irrigated_grapes",
        "groundnuts",
        "irrigated_groundnuts",
        "millet",
        "irrigated_millet",
        "oilpalm",
        "irrigated_oilpalm",
        "potatoes",
        "irrigated_potatoes",
        "pulses",
        "irrigated_pulses",
        "rapeseed",
        "irrigated_rapeseed",
        "rice",
        "irrigated_rice",
        "sorghum",
        "irrigated_sorghum",
        "sugarbeet",
        "irrigated_sugarbeet",
        "sugarcane",
        "irrigated_sugarcane",
        "sunflower",
        "irrigated_sunflower",
        "miscanthus",
        "irrigated_miscanthus",
        "switchgrass",
        "irrigated_switchgrass",
        "tropical_corn",
        "irrigated_tropical_corn",
        "tropical_soybean",
        "irrigated_tropical_soybean",
    ]
    return pftlist


def ivt_str2int(ivt_str):
    """
    Get CLM ivt number corresponding to a given name
    """
    pftlist = define_pftlist()
    if isinstance(ivt_str, str):
        ivt_int = pftlist.index(ivt_str)
    elif isinstance(ivt_str, list) or isinstance(ivt_str, np.ndarray):
        ivt_int = [ivt_str2int(x) for x in ivt_str]
        if isinstance(ivt_str, np.ndarray):
            ivt_int = np.array(ivt_int)
    else:
        raise RuntimeError(
            f"Update ivt_str_to_int() to handle input of type {type(ivt_str)} (if possible)"
        )

    return ivt_int


def ivt_int2str(ivt_int):
    """
    Get CLM ivt name corresponding to a given number
    """
    pftlist = define_pftlist()
    if np.issubdtype(type(ivt_int), np.integer) or int(ivt_int) == ivt_int:
        ivt_str = pftlist[int(ivt_int)]
    elif isinstance(ivt_int, list) or isinstance(ivt_int, np.ndarray):
        ivt_str = [ivt_int2str(x) for x in ivt_int]
        if isinstance(ivt_int, np.ndarray):
            ivt_str = np.array(ivt_str)
    elif isinstance(ivt_int, float):
        raise RuntimeError("List indices must be integers")
    else:
        raise RuntimeError(
            f"Update ivt_str_to_int() to handle input of type {type(ivt_int)} (if possible)"
        )

    return ivt_str


def is_this_vegtype(this_vegtype, this_filter, this_method):
    """
    Does this vegetation type's name match (for a given comparison method) any member of a filtering
    list?

    Methods:
    ok_contains:    True if any member of this_filter is found in this_vegtype.
    notok_contains: True of no member of this_filter is found in this_vegtype.
    ok_exact:       True if this_vegtype matches any member of this_filter
                    exactly.
    notok_exact:    True if this_vegtype does not match any member of
                    this_filter exactly.
    """
    # Make sure data type of this_vegtype is acceptable
    if isinstance(this_vegtype, float) and int(this_vegtype) == this_vegtype:
        this_vegtype = int(this_vegtype)
    data_type_ok = lambda x: isinstance(x, str) or isinstance(x, int) or isinstance(x, np.int64)
    ok_input = True
    if not data_type_ok(this_vegtype):
        if isinstance(this_vegtype, xr.core.dataarray.DataArray):
            this_vegtype = this_vegtype.values
        if isinstance(this_vegtype, (list, np.ndarray)):
            if len(this_vegtype) == 1 and data_type_ok(this_vegtype[0]):
                this_vegtype = this_vegtype[0]
            elif data_type_ok(this_vegtype[0]):
                raise TypeError(
                    "is_this_vegtype(): this_vegtype must be a single string or integer, not a list"
                    " of them. Did you mean to call is_each_vegtype() instead?"
                )
            else:
                ok_input = False
        else:
            ok_input = False
    if not ok_input:
        raise TypeError(
            "is_this_vegtype(): First argument (this_vegtype) must be a string or integer, not"
            f" {type(this_vegtype)}"
        )

    # Make sure data type of this_filter is acceptable
    if not np.iterable(this_filter):
        raise TypeError(
            "is_this_vegtype(): Second argument (this_filter) must be iterable (e.g., a list), not"
            f" {type(this_filter)}"
        )

    # Perform the comparison
    if this_method == "ok_contains":
        return any(n in this_vegtype for n in this_filter)
    if this_method == "notok_contains":
        return not any(n in this_vegtype for n in this_filter)
    if this_method == "ok_exact":
        return any(n == this_vegtype for n in this_filter)
    if this_method == "notok_exact":
        return not any(n == this_vegtype for n in this_filter)
    raise ValueError(f"Unknown comparison method: '{this_method}'")


def is_each_vegtype(this_vegtypelist, this_filter, this_method):
    """
    Get boolean list of whether each vegetation type in list is a managed crop

    this_vegtypelist: The list of vegetation types whose members you want to test.
    this_filter:      The list of strings against which you want to compare each member of
                      this_vegtypelist.
    this_method:      How you want to do the comparison. See is_this_vegtype().
    """
    if isinstance(this_vegtypelist, xr.DataArray):
        this_vegtypelist = this_vegtypelist.values

    return [is_this_vegtype(x, this_filter, this_method) for x in this_vegtypelist]


def define_mgdcrop_list():
    """
    List (strings) of managed crops in CLM.
    """
    notcrop_list = ["tree", "grass", "shrub", "unmanaged", "not_vegetated"]
    defined_pftlist = define_pftlist()
    is_crop = is_each_vegtype(defined_pftlist, notcrop_list, "notok_contains")
    return [defined_pftlist[i] for i, x in enumerate(is_crop) if x]


def vegtype_str2int(vegtype_str, vegtype_mainlist=None):
    """
    Convert list of vegtype strings to integer index equivalents.
    """
    convert_to_ndarray = not isinstance(vegtype_str, np.ndarray)
    if convert_to_ndarray:
        vegtype_str = np.array(vegtype_str)

    if isinstance(vegtype_mainlist, xr.Dataset):
        vegtype_mainlist = vegtype_mainlist.vegtype_str.values
    elif isinstance(vegtype_mainlist, xr.DataArray):
        vegtype_mainlist = vegtype_mainlist.values
    elif vegtype_mainlist is None:
        vegtype_mainlist = define_pftlist()
    if not isinstance(vegtype_mainlist, list) and isinstance(vegtype_mainlist[0], str):
        if isinstance(vegtype_mainlist, list):
            raise TypeError(
                f"Not sure how to handle vegtype_mainlist as list of {type(vegtype_mainlist[0])}"
            )
        raise TypeError(
            f"Not sure how to handle vegtype_mainlist as type {type(vegtype_mainlist[0])}"
        )

    if vegtype_str.shape == ():
        indices = np.array([-1])
    else:
        indices = np.full(len(vegtype_str), -1)
    for vegtype_str_2 in np.unique(vegtype_str):
        indices[np.where(vegtype_str == vegtype_str_2)] = vegtype_mainlist.index(vegtype_str_2)
    if convert_to_ndarray:
        indices = [int(x) for x in indices]
    return indices


def get_patch_ivts(this_ds, this_pftlist):
    """
    Get PFT of each patch, in both integer and string forms.
    """
    # First, get all the integer values; should be time*pft or pft*time. We will eventually just take the first timestep.
    vegtype_int = this_ds.patches1d_itype_veg
    vegtype_int.values = vegtype_int.values.astype(int)

    # Convert to strings.
    vegtype_str = list(np.array(this_pftlist)[vegtype_int.values])

    # Return a dictionary with both results
    return {"int": vegtype_int, "str": vegtype_str, "all_str": this_pftlist}


def get_vegtype_str_da(vegtype_str):
    """
    Convert a list of strings with vegetation type names into a DataArray. Used to add vegetation type info in import_ds().
    """
    nvt = len(vegtype_str)
    vegtype_str_da = xr.DataArray(
        vegtype_str, coords={"ivt": np.arange(0, nvt)}, dims=["ivt"], name="vegtype_str"
    )
    return vegtype_str_da


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

        # Get any _1d variables that are associated with those dimensions. These will be useful in gridding. Also, if any dimension is "pft", set up to rename it and all like-named variables to "patch"
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
        this_pftlist = define_pftlist()
        get_patch_ivts(
            ds_in, this_pftlist
        )  # Includes check of whether vegtype changes over time anywhere
        vegtype_da = get_vegtype_str_da(this_pftlist)
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
        ds_in = safer_timeslice(ds_in, time_slice)

    # Finish import
    ds_in = xr.decode_cf(ds_in, decode_times=True)

    # Compute derived variables
    for var in derived_vars:
        if (
            var == "HYEARS"
            and "HDATES" in ds_in
            and ds_in.HDATES.dims == ("time", "mxharvests", "patch")
        ):
            year_list = np.array([np.float32(x.year - 1) for x in ds_in.time.values])
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

    return ds_in


def import_ds(
    filelist,
    my_vars=None,
    my_vegtypes=None,
    time_slice=None,
    my_vars_missing_ok=[],
    only_active_patches=False,
    rename_lsmlatlon=False,
    chunks=None,
):
    """
    Import a dataset that can be spread over multiple files, only including specified variables
    and/or vegetation types and/or timesteps, concatenating by time.

    - DOES actually read the dataset into memory, but only AFTER dropping unwanted variables and/or
      vegetation types.
    """
    # Convert my_vegtypes here, if needed, to avoid repeating the process each time you read a file in xr.open_mfdataset().
    if my_vegtypes is not None:
        if not isinstance(my_vegtypes, list):
            my_vegtypes = [my_vegtypes]
        if isinstance(my_vegtypes[0], str):
            my_vegtypes = vegtype_str2int(my_vegtypes)

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

    # Remove files from list if they don't contain requested timesteps.
    # time_slice should be in the format slice(start,end[,step]). start or end can be None to be unbounded on one side. Note that the standard slice() documentation suggests that only elements through end-1 will be selected, but that seems not to be the case in the xarray implementation.
    if time_slice:
        new_filelist = []
        for file in sorted(filelist):
            filetime = xr.open_dataset(file).time
            filetime_sel = safer_timeslice(filetime, time_slice)
            include_this_file = filetime_sel.size
            if include_this_file:
                new_filelist.append(file)

            # If you found some matching files, but then you find one that doesn't, stop going through the list.
            elif new_filelist:
                break
        if not new_filelist:
            raise RuntimeError(f"No files found in time_slice {time_slice}")
        filelist = new_filelist

    # The xarray open_mfdataset() "preprocess" argument requires a function that takes exactly one variable (an xarray.Dataset object). Wrapping mfdataset_preproc() in this lambda function allows this. Could also just allow mfdataset_preproc() to access my_vars and my_vegtypes directly, but that's bad practice as it could lead to scoping issues.
    mfdataset_preproc_closure = lambda ds: mfdataset_preproc(ds, my_vars, my_vegtypes, time_slice)

    # Import
    if isinstance(filelist, list) and len(filelist) == 1:
        filelist = filelist[0]
    if isinstance(filelist, list):
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=DeprecationWarning)
            if importlib.find_loader("dask") is None:
                raise ModuleNotFoundError(
                    "You have asked xarray to import a list of files as a single Dataset using"
                    " open_mfdataset(), but this requires dask, which is not available.\nFile"
                    f" list: {filelist}"
                )
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

    # Include only active patches (or whatever)
    if only_active_patches:
        is_active = this_ds.patches1d_active.values
        p_active = np.where(is_active)[0]
        this_ds_active = this_ds.isel(patch=p_active)

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


def get_thisvar_da(var, this_ds):
    """
    Return a DataArray, with defined coordinates, for a given variable in a dataset.
    """
    # Make DataArray for this variable
    thisvar_da = np.array(this_ds.variables[var])
    these_dims = this_ds.variables[var].dims
    thisvar_da = xr.DataArray(thisvar_da, dims=these_dims)

    # Define coordinates of this variable's DataArray
    dims_dict = dict()
    for dim in these_dims:
        dims_dict[dim] = this_ds[dim]
    thisvar_da = thisvar_da.assign_coords(dims_dict)
    thisvar_da.attrs = this_ds[var].attrs

    return thisvar_da


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

    if not fill_value and "_FillValue" in thisvar_da.attrs:
        fill_value = thisvar_da.attrs["_FillValue"]

    # Renumber vt_da to work as indices on new ivt dimension, if needed.
    ### Ensures that the unique set of vt_da values begins with 1 and
    ### contains no missing steps.
    if "ivt" in this_ds and vt_da is not None:
        vt_da.values = np.array([np.where(this_ds.ivt.values == x)[0][0] for x in vt_da.values])

    # Get new dimension list
    new_dims = list(thisvar_da.dims)
    ### Remove "[spatial_unit]".
    if spatial_unit in new_dims:
        new_dims.remove(spatial_unit)
    #  Add "ivt_str" (vegetation type, as string). This needs to go at the end, to avoid a possible situation where you wind up with multiple Ellipsis members of fill_indices.
    if "ivt" in this_ds and spatial_unit == "patch":
        new_dims.append("ivt_str")
    ### Add lat and lon to end of list
    new_dims = new_dims + ["lat", "lon"]

    # Set up empty array
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

    # Fill with this variable
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

    # Assign coordinates, attributes and name
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


def safer_timeslice(ds_in, time_slice, time_var="time"):
    """
    ctsm_pylib can't handle time slicing like Dataset.sel(time=slice("1998-01-01", "2005-12-31"))
    for some reason. This function tries to fall back to slicing by integers. It should work with
    both Datasets and DataArrays.
    """
    try:
        ds_in = ds_in.sel({time_var: time_slice})
    except:  # pylint: disable=bare-except
        # If the issue might have been slicing using strings, try to fall back to integer slicing
        if (
            isinstance(time_slice.start, str)
            and isinstance(time_slice.stop, str)
            and len(time_slice.start.split("-")) == 3
            and time_slice.start.split("-")[1:] == ["01", "01"]
            and len(time_slice.stop.split("-")) == 3
            and (
                time_slice.stop.split("-")[1:] == ["12", "31"]
                or time_slice.stop.split("-")[1:] == ["01", "01"]
            )
        ):
            fileyears = np.array([x.year for x in ds_in.time.values])
            if len(np.unique(fileyears)) != len(fileyears):
                print("Could not fall back to integer slicing of years: Time axis not annual")
                raise
            y_start = int(time_slice.start.split("-")[0])
            y_stop = int(time_slice.stop.split("-")[0])
            where_in_timeslice = np.where((fileyears >= y_start) & (fileyears <= y_stop))[0]
            ds_in = ds_in.isel({time_var: where_in_timeslice})
        else:
            print(f"Could not fall back to integer slicing for time_slice {time_slice}")
            raise

    return ds_in


def lon_idl2pm(lons_in, fail_silently=False):
    """
    Convert a longitude axis that's -180 to 180 around the international date line to one that's 0
    to 360 around the prime meridian.

    - If you pass in a Dataset or DataArray, the "lon" coordinates will be changed. Otherwise, it
      assumes you're passing in numeric data.
    """

    def check_ok(tmp, fail_silently):
        msg = ""

        if np.any(tmp > 180):
            msg = f"Maximum longitude is already > 180 ({np.max(tmp)})"
        elif np.any(tmp < -180):
            msg = f"Minimum longitude is < -180 ({np.min(tmp)})"

        if msg == "":
            return True
        if fail_silently:
            return False
        raise ValueError(msg)

    def do_it(tmp):
        tmp = tmp + 360
        tmp = np.mod(tmp, 360)
        return tmp

    if isinstance(lons_in, (xr.DataArray, xr.Dataset)):
        if not check_ok(lons_in.lon.values, fail_silently):
            return lons_in
        lons_out = lons_in
        lons_out = lons_out.assign_coords(lon=do_it(lons_in.lon.values))
        lons_out = make_lon_increasing(lons_out)
    else:
        if not check_ok(lons_in, fail_silently):
            return lons_in
        lons_out = do_it(lons_in)
        if not is_strictly_increasing(lons_out):
            print(
                "WARNING: You passed in numeric longitudes to lon_idl2pm() and these have been"
                " converted, but they're not strictly increasing."
            )
        print(
            "To assign the new longitude coordinates to an Xarray object, use"
            " xarrayobject.assign_coordinates()! (Pass the object directly in to lon_idl2pm() in"
            " order to suppress this message.)"
        )

    return lons_out


def is_strictly_increasing(this_list):
    """
    Helper function to check that a list is strictly increasing

    https://stackoverflow.com/a/4983359/2965321
    """
    return all(x < y for x, y in zip(this_list, this_list[1:]))


def make_lon_increasing(xr_obj):
    """
    Ensure that longitude axis coordinates are monotonically increasing
    """
    if not "lon" in xr_obj.dims:
        return xr_obj

    lons = xr_obj.lon.values
    if is_strictly_increasing(lons):
        return xr_obj

    shift = 0
    while not is_strictly_increasing(lons) and shift < lons.size:
        shift = shift + 1
        lons = np.roll(lons, 1, axis=0)
    if not is_strictly_increasing(lons):
        raise RuntimeError("Unable to rearrange longitude axis so it's monotonically increasing")

    return xr_obj.roll(lon=shift, roll_coords=True)
