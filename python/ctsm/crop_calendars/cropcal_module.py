"""
Helper functions for various crop calendar stuff
"""

import os
import glob
import numpy as np
import xarray as xr

import ctsm.crop_calendars.cropcal_utils as utils
from ctsm.crop_calendars.convert_axis_time2gs import convert_axis_time2gs
from ctsm.crop_calendars.check_rx_obeyed import check_rx_obeyed
from ctsm.crop_calendars.cropcal_constants import DEFAULT_GDD_MIN
from ctsm.crop_calendars.import_ds import import_ds
from ctsm.utils import is_instantaneous

MISSING_RX_GDD_VAL = -1


def check_and_trim_years(year_1, year_n, ds_in):
    """
    After importing a file, restrict it to years of interest.
    """

    # Check that all desired years are included
    year = utils.get_timestep_year(ds_in, ds_in.time.values[0])
    if year > year_1:
        raise RuntimeError(f"Requested year_1 is {year_1} but first year in outputs is {year}")
    year = utils.get_timestep_year(ds_in, ds_in.time.values[-1])
    if year < year_1:
        raise RuntimeError(f"Requested year_n is {year_n} but last year in outputs is {year}")

    # Remove years outside range of interest
    ### Include an extra year at the end to finish out final seasons.
    slice_yr_1 = year_1
    slice_yr_n = year_n + 1
    if is_instantaneous(ds_in["time"]):
        slice_yr_1 += 1
        slice_yr_n += 1
    ds_in = utils.safer_timeslice(ds_in, slice(f"{slice_yr_1}-01-01", f"{slice_yr_n}-12-31"))

    # Make sure you have the expected number of timesteps (including extra year)
    n_years_expected = year_n - year_1 + 2
    if ds_in.dims["time"] != n_years_expected:
        raise RuntimeError(
            f"Expected {n_years_expected} timesteps in output but got {ds_in.dims['time']}"
        )

    return ds_in


def open_lu_ds(filename, year_1, year_n, existing_ds, ungrid=True):
    """
    Open land-use dataset
    """
    # Open and trim to years of interest
    this_ds_gridded = xr.open_dataset(filename).sel(time=slice(year_1, year_n))

    # Assign actual lon/lat coordinates
    this_ds_gridded = this_ds_gridded.assign_coords(
        lon=("lsmlon", existing_ds.lon.values), lat=("lsmlat", existing_ds.lat.values)
    )
    this_ds_gridded = this_ds_gridded.swap_dims({"lsmlon": "lon", "lsmlat": "lat"})

    if "AREA" in this_ds_gridded:
        this_ds_gridded["AREA_CFT"] = (
            this_ds_gridded.AREA
            * 1e6
            * this_ds_gridded.LANDFRAC_PFT
            * this_ds_gridded.PCT_CROP
            / 100
            * this_ds_gridded.PCT_CFT
            / 100
        )
        this_ds_gridded["AREA_CFT"].attrs = {"units": "m2"}
        this_ds_gridded["AREA_CFT"].load()
    else:
        print("Warning: AREA missing from Dataset, so AREA_CFT will not be created")

    if not ungrid:
        return this_ds_gridded

    # Un-grid
    query_ilons = [int(x) - 1 for x in existing_ds["patches1d_ixy"].values]
    query_ilats = [int(x) - 1 for x in existing_ds["patches1d_jxy"].values]
    query_ivts = [
        list(this_ds_gridded.cft.values).index(x) for x in existing_ds["patches1d_itype_veg"].values
    ]

    this_ds = xr.Dataset(attrs=this_ds_gridded.attrs)
    for var in ["AREA", "LANDFRAC_PFT", "PCT_CFT", "PCT_CROP", "AREA_CFT"]:
        if var not in this_ds_gridded:
            continue
        if "time" in this_ds_gridded[var].dims:
            new_coords = existing_ds["GRAINC_TO_FOOD_ANN"].coords
        else:
            new_coords = existing_ds["patches1d_lon"].coords
        if "cft" in this_ds_gridded[var].dims:
            this_ds[var] = (
                this_ds_gridded[var]
                .isel(
                    lon=xr.DataArray(query_ilons, dims="patch"),
                    lat=xr.DataArray(query_ilats, dims="patch"),
                    cft=xr.DataArray(query_ivts, dims="patch"),
                    drop=True,
                )
                .assign_coords(new_coords)
            )
        else:
            this_ds[var] = (
                this_ds_gridded[var]
                .isel(
                    lon=xr.DataArray(query_ilons, dims="patch"),
                    lat=xr.DataArray(query_ilats, dims="patch"),
                    drop=True,
                )
                .assign_coords(new_coords)
            )
    for var in existing_ds:
        if "patches1d_" in var or "grid1d_" in var:
            this_ds[var] = existing_ds[var]
    this_ds["lon"] = this_ds_gridded["lon"]
    this_ds["lat"] = this_ds_gridded["lat"]

    # Which crops are irrigated?
    is_irrigated = np.full_like(this_ds["patches1d_itype_veg"], False)
    for vegtype_str in np.unique(this_ds["patches1d_itype_veg_str"].values):
        if "irrigated" not in vegtype_str:
            continue
        vegtype_int = utils.ivt_str2int(vegtype_str)
        is_this_vegtype = np.where(this_ds["patches1d_itype_veg"].values == vegtype_int)[0]
        is_irrigated[is_this_vegtype] = True
    this_ds["IRRIGATED"] = xr.DataArray(
        data=is_irrigated,
        coords=this_ds["patches1d_itype_veg_str"].coords,
        attrs={"long_name": "Is patch irrigated?"},
    )

    # How much area is irrigated?
    this_ds["IRRIGATED_AREA_CFT"] = this_ds["IRRIGATED"] * this_ds["AREA_CFT"]
    this_ds["IRRIGATED_AREA_CFT"].attrs = {
        "long name": "CFT area (irrigated types only)",
        "units": "m^2",
    }
    this_ds["IRRIGATED_AREA_GRID"] = (
        this_ds["IRRIGATED_AREA_CFT"]
        .groupby(this_ds["patches1d_gi"])
        .sum()
        .rename({"patches1d_gi": "gridcell"})
    )
    this_ds["IRRIGATED_AREA_GRID"].attrs = {
        "long name": "Irrigated area in gridcell",
        "units": "m^2",
    }

    return this_ds


def check_v0_le_v1(this_ds, var_list, msg_txt=" ", both_nan_ok=False, throw_error=False):
    """
    Make sure that, e.g., GDDACCUM_PERHARV is always <= HUI_PERHARV
    """
    var0 = var_list[0]
    var1 = var_list[1]
    gdd_lt_hui = this_ds[var0] <= this_ds[var1]
    if both_nan_ok:
        gdd_lt_hui = gdd_lt_hui | (np.isnan(this_ds[var0]) & np.isnan(this_ds[var1]))
    if np.all(gdd_lt_hui):
        any_bad = False
        print(f"✅{msg_txt}{var0} always <= {var1}")
    else:
        any_bad = True
        msg = f"❌{msg_txt}{var0} *not* always <= {var1}"
        gdd_lt_hui_vals = gdd_lt_hui.values
        patch_index = np.where(~gdd_lt_hui_vals)[0][0]
        msg = (
            msg
            + f"\ne.g., patch {patch_index}: {this_ds.patches1d_itype_veg_str.values[patch_index]},"
            + f" lon {this_ds.patches1d_lon.values[patch_index]} lat "
            + f"{this_ds.patches1d_lat.values[patch_index]}:"
        )
        msg = msg + f"\n{this_ds[var0].values[patch_index,:]}"
        msg = msg + f"\n{this_ds[var1].values[patch_index,:]}"
        if throw_error:
            print(msg)
        else:
            raise RuntimeError(msg)
    return any_bad


def get_gs_len_da(this_da):
    """
    Get growing season lengths from a DataArray of hdate-sdate
    """
    tmp = this_da.values
    with np.errstate(invalid="ignore"):
        tmp_lt_0 = tmp < 0
    tmp[tmp_lt_0] = 365 + tmp[tmp_lt_0]
    this_da.values = tmp
    this_da.attrs["units"] = "days"
    return this_da


def import_max_gs_length(paramfile_dir, my_clm_ver, my_clm_subver):
    """
    Import maximum growing season length
    """
    # Get parameter file
    pattern = os.path.join(paramfile_dir, f"*{my_clm_ver}_params.{my_clm_subver}.nc")
    paramfile = glob.glob(pattern)
    if len(paramfile) != 1:
        raise RuntimeError(f"Expected to find 1 match of {pattern}; found {len(paramfile)}")
    paramfile_ds = xr.open_dataset(paramfile[0])

    # Import max growing season length (stored in netCDF as nanoseconds!)
    paramfile_mxmats = paramfile_ds["mxmat"].values / np.timedelta64(1, "D")

    # Import PFT name list
    paramfile_pftnames = [
        x.decode("UTF-8").replace(" ", "") for x in paramfile_ds["pftname"].values
    ]

    # Build dict
    mxmat_dict = {}
    for i, pftname in enumerate(paramfile_pftnames):
        mxmat = paramfile_mxmats[i]
        if not np.isnan(mxmat):
            mxmat_dict[pftname] = int(mxmat)
        else:
            mxmat_dict[pftname] = np.inf

    return mxmat_dict


def unexpected_negative_rx_gdd(data_array):
    """
    Return True if there's a negative value not matching the designated missing value
    """
    return np.any((data_array.values < 0) & (data_array.values != MISSING_RX_GDD_VAL))


def import_rx_dates(var_prefix, date_infile, dates_ds, set_neg1_to_nan=True):
    """
    Import prescribed sowing/harvest dates

    E.g. import_rx_dates("sdate", sdates_rx_file, dates_ds0_orig)
    """
    # Get run info:
    # Max number of growing seasons per year
    if "mxsowings" in dates_ds:
        mxsowings = dates_ds.dims["mxsowings"]
    else:
        mxsowings = 1

    # Which vegetation types were simulated?
    itype_veg_to_import = np.unique(dates_ds.patches1d_itype_veg)

    date_varlist = []
    for i in itype_veg_to_import:
        for j in np.arange(mxsowings):
            this_var = f"{var_prefix}{j+1}_{i}"
            date_varlist = date_varlist + [this_var]

    this_ds = import_ds(date_infile, my_vars=date_varlist)

    did_warn = False
    for var in this_ds:
        v_new = var.replace(var_prefix, "gs")
        this_ds = this_ds.rename({var: v_new})

        # Set GDD values matching MISSING_RX_GDD_VAL to NaN. Only warn the first time.
        if (
            set_neg1_to_nan
            and var_prefix == "gdd"
            and v_new != var
            and np.any(this_ds[v_new].values < 0)
        ):
            if unexpected_negative_rx_gdd(this_ds[v_new]):
                raise RuntimeError(f"Unexpected negative value in {var}")
            if not did_warn:
                print(f"Setting {MISSING_RX_GDD_VAL} rx GDD values to NaN")
                did_warn = True
            this_ds[v_new] = this_ds[v_new].where(this_ds[v_new] != MISSING_RX_GDD_VAL)

    return this_ds


def check_no_negative(this_ds_in, varlist_no_negative, which_file, verbose):
    """
    In import_output(), check that there are no unexpected negative values.
    """
    tiny_neg_ok = 1e-12
    this_ds = this_ds_in.copy()
    for var in this_ds:
        if not any(x in var for x in varlist_no_negative):
            continue
        the_min = np.nanmin(this_ds[var].values)
        if the_min < 0:
            if np.abs(the_min) <= tiny_neg_ok:
                if verbose:
                    print(
                        f"Tiny negative value(s) in {var} (abs <= {tiny_neg_ok}) being set to 0"
                        + f" ({which_file})"
                    )
            else:
                print(
                    f"WARNING: Unexpected negative value(s) in {var}; minimum {the_min} "
                    + f"({which_file})"
                )
            values = this_ds[var].copy().values
            with np.errstate(invalid="ignore"):
                do_setto_0 = (values < 0) & (values >= -tiny_neg_ok)
            values[np.where(do_setto_0)] = 0
            this_ds[var] = xr.DataArray(
                values,
                coords=this_ds[var].coords,
                dims=this_ds[var].dims,
                attrs=this_ds[var].attrs,
            )

        elif verbose:
            print(f"No negative value(s) in {var}; min {the_min} ({which_file})")
    return this_ds


def check_no_zeros(this_ds, varlist_no_zero, which_file, verbose):
    """
    In import_output(), check that there are no unexpected zeros.
    """
    for var in this_ds:
        if not any(x in var for x in varlist_no_zero):
            continue
        if np.any(this_ds[var].values == 0):
            print(f"WARNING: Unexpected zero(s) in {var} ({which_file})")
        elif verbose:
            print(f"No zero value(s) in {var} ({which_file})")


def import_output(
    filename,
    my_vars,
    *,
    year_1=None,
    year_n=None,
    my_vegtypes=utils.define_mgdcrop_list_withgrasses(),
    sdates_rx_ds=None,
    gdds_rx_ds=None,
    verbose=False,
    throw_errors=True,
):
    """
    Import CLM output
    """
    any_bad = False

    # Import
    this_ds = import_ds(filename, my_vars=my_vars, my_vegtypes=my_vegtypes)

    # Trim to years of interest (do not include extra year needed for finishing last growing season)
    if year_1 and year_n:
        this_ds = check_and_trim_years(year_1, year_n, this_ds)
    else:  # Assume including all growing seasons except last complete one are "of interest"
        year_1 = this_ds.time.values[0].year
        year_n = this_ds.time.values[-1].year - 2
        this_ds = check_and_trim_years(year_1, year_n, this_ds)

    # What vegetation types are included?
    vegtype_list = [
        x for x in this_ds.vegtype_str.values if x in this_ds.patches1d_itype_veg_str.values
    ]

    # Check for consistency among sowing/harvest date/year info
    date_vars = ["SDATES_PERHARV", "SYEARS_PERHARV", "HDATES", "HYEARS"]
    all_nan = np.full(this_ds[date_vars[0]].shape, True)
    all_nonpos = np.full(this_ds[date_vars[0]].shape, True)
    all_pos = np.full(this_ds[date_vars[0]].shape, True)
    for var in date_vars:
        all_nan = all_nan & np.isnan(this_ds[var].values)
        with np.errstate(invalid="ignore"):
            all_nonpos = all_nonpos & (this_ds[var].values <= 0)
            all_pos = all_pos & (this_ds[var].values > 0)
    if np.any(np.bitwise_not(all_nan | all_nonpos | all_pos)):
        raise RuntimeError("Inconsistent missing/present values on mxharvests axis")

    # When doing transient runs, it's somehow possible for crops in newly-active patches to be
    # *already alive*. They even have a sowing date (idop)! This will of course not show up in
    # SDATES, but it does show up in SDATES_PERHARV.
    # I could put the SDATES_PERHARV dates into where they "should" be, but instead I'm just going
    # to invalidate those "seasons."
    this_ds = handle_zombie_crops(this_ds)

    # Check for no zero values where there shouldn't be
    varlist_no_zero = ["DATE", "YEAR"]
    check_no_zeros(this_ds, varlist_no_zero, "original file", verbose)

    # Convert time*mxharvests axes to growingseason axis
    this_ds_gs = convert_axis_time2gs(this_ds, verbose=verbose, incl_orig=False)

    # These are needed for calculating yield later
    this_ds_gs["GRAINC_TO_FOOD_PERHARV"] = this_ds["GRAINC_TO_FOOD_PERHARV"]
    this_ds_gs["GDDHARV_PERHARV"] = this_ds["GDDHARV_PERHARV"]

    # Get growing season length
    this_ds["GSLEN_PERHARV"] = get_gs_len_da(this_ds["HDATES"] - this_ds["SDATES_PERHARV"])
    this_ds_gs["GSLEN"] = get_gs_len_da(this_ds_gs["HDATES"] - this_ds_gs["SDATES"])
    this_ds_gs["GSLEN_PERHARV"] = this_ds["GSLEN_PERHARV"]

    # Get HUI accumulation as fraction of required
    this_ds_gs["HUIFRAC"] = this_ds_gs["HUI"] / this_ds_gs["GDDHARV"]
    this_ds_gs["HUIFRAC_PERHARV"] = this_ds["HUI_PERHARV"] / this_ds["GDDHARV_PERHARV"]
    for var in ["HUIFRAC", "HUIFRAC_PERHARV"]:
        this_ds_gs[var].attrs["units"] = "Fraction of required"

    # Avoid tiny negative values
    varlist_no_negative = ["GRAIN", "REASON", "GDD", "HUI", "YEAR", "DATE", "GSLEN"]
    this_ds_gs = check_no_negative(this_ds_gs, varlist_no_negative, "new file", verbose)

    # Check for no zero values where there shouldn't be
    varlist_no_zero = ["REASON", "DATE"]
    check_no_zeros(this_ds_gs, varlist_no_zero, "new file", verbose)

    # Check that e.g., GDDACCUM <= HUI
    for var_list in [["GDDACCUM", "HUI"], ["SYEARS", "HYEARS"]]:
        if all(v in this_ds_gs for v in var_list):
            any_bad = any_bad or check_v0_le_v1(
                this_ds_gs, var_list, both_nan_ok=True, throw_error=throw_errors
            )

    # Check that prescribed calendars were obeyed
    if sdates_rx_ds:
        check_rx_obeyed(vegtype_list, sdates_rx_ds, this_ds, "this_ds", "SDATES")
    if gdds_rx_ds:
        check_rx_obeyed(
            vegtype_list,
            gdds_rx_ds,
            this_ds,
            "this_ds",
            "GDDHARV",
            gdd_min=DEFAULT_GDD_MIN,
        )

    # Convert time axis to integer year, saving original as 'cftime'
    this_ds_gs = convert_time_to_int_year(filename, this_ds, this_ds_gs)

    # Get number of harvests
    this_ds_gs["NHARVESTS"] = (this_ds_gs["GDDHARV_PERHARV"] > 0).sum(dim="mxharvests")
    # Get number of harvests that would be missed if only seeing max 1 per calendar year
    if np.any(this_ds_gs["NHARVESTS"] > 2):
        raise RuntimeError("How to get NHARVEST_DISCREP for NHARVESTS > 2?")
    this_ds_gs["NHARVEST_DISCREP"] = (this_ds_gs["NHARVESTS"] == 2).astype(int)

    return this_ds_gs, any_bad


def convert_time_to_int_year(filename, this_ds, this_ds_gs):
    """
    Convert time axis to integer year, saving original as 'cftime'
    """
    if "time_bounds" in this_ds:
        # Always true before PR #2838, when even files with all instantaneous variables got
        # time_bounds saved. After that PR (and before the segregation of instantaneous and other
        # variables onto separate files), files with an instantaneous variable first in their list
        # do not get time_bounds saved.
        this_ds_gs = this_ds_gs.assign_coords({"cftime": this_ds["time_bounds"].isel({"nbnd": 0})})
        this_ds_gs = this_ds_gs.assign_coords(
            {"time": [t.year for t in this_ds_gs["cftime"].values]}
        )
    elif this_ds["time"].attrs["long_name"] == "time at end of time step":
        # This is an "instantaneous file."
        this_ds_gs = this_ds_gs.assign_coords({"cftime": this_ds["time"]})
        this_ds_gs = this_ds_gs.assign_coords(
            {"time": [t.year - 1 for t in this_ds_gs["cftime"].values]}
        )
    else:
        raise RuntimeError(
            f"{filename} is neither an instantaneous nor a combined/non-instantaneous file."
        )

    return this_ds_gs


def handle_zombie_crops(this_ds):
    """
    When doing transient runs, it's somehow possible for crops in newly-active patches to be
    *already alive*. They even have a sowing date (idop)! This will of course not show up in
    SDATES, but it does show up in SDATES_PERHARV.
    I could put the SDATES_PERHARV dates into where they "should" be, but instead I'm just going
    to invalidate those "seasons."
    """
    # In all but the last calendar year, which patches had no sowing?
    no_sowing_yp = np.all(np.isnan(this_ds.SDATES.values[:-1, :, :]), axis=1)
    # In all but the first calendar year, which harvests' jdays are < their sowings' jdays?
    # (Indicates sowing the previous calendar year.)
    with np.errstate(invalid="ignore"):
        hsdate1_gt_hdate1_yp = (
            this_ds.SDATES_PERHARV.values[1:, 0, :] > this_ds.HDATES.values[1:, 0, :]
        )
    # Where both, we have the problem.
    falsely_alive_yp = no_sowing_yp & hsdate1_gt_hdate1_yp
    if np.any(falsely_alive_yp):
        print(
            f"Warning: {np.sum(falsely_alive_yp)} patch-seasons being ignored: Seemingly sown the "
            + "year before harvest, but no sowings occurred that year."
        )
        falsely_alive_yp = np.concatenate(
            (np.full((1, this_ds.dims["patch"]), False), falsely_alive_yp), axis=0
        )
        falsely_alive_y1p = np.expand_dims(falsely_alive_yp, axis=1)
        dummy_false_y1p = np.expand_dims(np.full_like(falsely_alive_yp, False), axis=1)
        falsely_alive_yhp = np.concatenate((falsely_alive_y1p, dummy_false_y1p), axis=1)
        for var in this_ds.data_vars:
            if this_ds[var].dims != ("time", "mxharvests", "patch"):
                continue
            this_ds[var] = this_ds[var].where(~falsely_alive_yhp)
    return this_ds
