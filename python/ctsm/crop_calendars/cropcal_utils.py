"""
utility functions
copied from klindsay, https://github.com/klindsay28/CESM2_coup_carb_cycle_JAMES/blob/master/utils.py
"""

import numpy as np
import xarray as xr

from ctsm.utils import is_instantaneous


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
    elif isinstance(ivt_str, (list, np.ndarray)):
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
    elif isinstance(ivt_int, (list, np.ndarray)):
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
    def data_type_ok(x):
        return isinstance(x, (int, np.int64, str))

    if isinstance(this_vegtype, float) and int(this_vegtype) == this_vegtype:
        this_vegtype = int(this_vegtype)
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


def define_crop_list():
    """
    List (strings) of managed crops in CLM.
    """
    notcrop_list = [
        "tree",
        "c3_arctic_grass",
        "c3_non-arctic_grass",
        "c4_grass",
        "shrub",
        "not_vegetated",
    ]
    defined_pftlist = define_pftlist()
    is_crop = is_each_vegtype(defined_pftlist, notcrop_list, "notok_contains")
    return [defined_pftlist[i] for i, x in enumerate(is_crop) if x]


def define_mgdcrop_list_nograsses():
    """
    List (strings) of managed crops in CLM.
    """
    notcrop_list = ["tree", "grass", "shrub", "unmanaged", "not_vegetated"]
    defined_pftlist = define_pftlist()
    is_crop = is_each_vegtype(defined_pftlist, notcrop_list, "notok_contains")
    return [defined_pftlist[i] for i, x in enumerate(is_crop) if x]


def define_mgdcrop_list_withgrasses():
    """
    List (strings) of managed crops in CLM.
    """
    notcrop_list = [
        "tree",
        "c3_arctic_grass",
        "c3_non-arctic_grass",
        "c4_grass",
        "shrub",
        "unmanaged",
        "not_vegetated",
    ]
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
    # First, get all the integer values; should be time*pft or pft*time. We will eventually just
    # take the first timestep.
    vegtype_int = this_ds.patches1d_itype_veg
    vegtype_int.values = vegtype_int.values.astype(int)

    # Convert to strings.
    vegtype_str = list(np.array(this_pftlist)[vegtype_int.values])

    # Return a dictionary with both results
    return {"int": vegtype_int, "str": vegtype_str, "all_str": this_pftlist}


def get_vegtype_str_da(vegtype_str):
    """
    Convert a list of strings with vegetation type names into a DataArray.
    """
    nvt = len(vegtype_str)
    vegtype_str_da = xr.DataArray(
        vegtype_str, coords={"ivt": np.arange(0, nvt)}, dims=["ivt"], name="vegtype_str"
    )
    return vegtype_str_da


def safer_timeslice(ds_in, time_slice, time_var="time"):
    """
    ctsm_pylib can't handle time slicing like Dataset.sel(time=slice("1998-01-01", "2005-12-31"))
    for some reason. This function tries to fall back to slicing by integers. It should work with
    both Datasets and DataArrays. NOTE: This isn't a problem for more modern Python environments.
    Even npl-2022b can use the straightforward slicing in the "try" block.
    """
    try:
        ds_in = ds_in.sel({time_var: time_slice})
    except Exception as this_exception:  # pylint: disable=broad-except
        # If the issue might have been slicing using strings, try to fall back to integer slicing
        can_try_integer_slicing = (
            isinstance(time_slice.start, str)
            and isinstance(time_slice.stop, str)
            and len(time_slice.start.split("-")) == 3
            and time_slice.start.split("-")[1:] == ["01", "01"]
            and len(time_slice.stop.split("-")) == 3
            and (
                time_slice.stop.split("-")[1:] == ["12", "31"]
                or time_slice.stop.split("-")[1:] == ["01", "01"]
            )
        )
        if can_try_integer_slicing:
            fileyears = np.array([x.year for x in ds_in.time.values])
            if len(np.unique(fileyears)) != len(fileyears):
                msg = "Could not fall back to integer slicing of years: Time axis not annual"
                raise RuntimeError(msg) from this_exception
            y_start = int(time_slice.start.split("-")[0])
            y_stop = int(time_slice.stop.split("-")[0])
            where_in_timeslice = np.where((fileyears >= y_start) & (fileyears <= y_stop))[0]
            ds_in = ds_in.isel({time_var: where_in_timeslice})
        else:
            msg = f"Could not fall back to integer slicing for time_slice {time_slice}"
            raise RuntimeError(msg) from this_exception

    return ds_in


def lon_axis_type180_to_type360(lons_in, fail_silently=False):
    """
    Convert a longitude axis that's -180 to 180 to one that's 0 to 360

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
                "WARNING: You passed in numeric longitudes to lon_axis_type180_to_type360() and"
                " these have been converted, but they're not strictly increasing."
            )
        print(
            "To assign the new longitude coordinates to an Xarray object, use"
            " xarrayobject.assign_coordinates()! (Pass the object directly in to"
            " lon_axis_type180_to_type360() in order to suppress this message.)"
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


def get_beg_inst_timestep_year(timestep):
    """
    Get year associated with the BEGINNING of a timestep in an
    instantaneous file
    """
    year = timestep.year

    is_jan1 = timestep.dayofyr == 1
    is_midnight = timestep.hour == timestep.minute == timestep.second == 0
    if is_jan1 and is_midnight:
        year -= 1

    return year


def get_timestep_year(dsa, timestep):
    """
    Get the year associated with a timestep, with different handling
    depending on whether the file is instantaneous
    """
    if is_instantaneous(dsa["time"]):
        year = get_beg_inst_timestep_year(timestep)
    else:
        year = timestep.year
    return year


def get_integer_years(dsa):
    """
    Convert time axis to numpy array of integer years
    """
    out_array = [get_timestep_year(dsa, t) for t in dsa["time"].values]
    return out_array
