"""
Helper functions for various crop calendar stuff
"""
# pylint: disable=too-many-lines

import warnings
import sys
import os
import glob
import numpy as np
import xarray as xr

# Import the CTSM Python utilities.
# sys.path.insert() is necessary for RXCROPMATURITY to work. The fact that it's calling this script
# in the RUN phase seems to require the python/ directory to be manually added to path.
_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)
import ctsm.crop_calendars.cropcal_utils as utils  # pylint: disable=wrong-import-position

try:
    import pandas as pd
except ModuleNotFoundError:
    pass


# Define conversion multipliers, {from: {to1, to2, ...}, ...}
multiplier_dict = {
    # Mass
    "g": {
        "Mt": 1e-12,
    },
    "t": {
        "Mt": 1e-6,
    },
    # Volume
    "m3": {
        "km3": 1e-9,
    },
    # Yield
    "g/m2": {
        "t/ha": 1e-6 * 1e4,
    },
}

# Minimum harvest threshold allowed in PlantCrop()
# Was 50 before cropcal runs 2023-01-28
DEFAULT_GDD_MIN = 1.0


def check_and_trim_years(year_1, year_n, ds_in):
    """
    After importing a file, restrict it to years of interest.
    """
    ### In annual outputs, file with name Y is actually results from year Y-1.
    ### Note that time values refer to when it was SAVED. So 1981-01-01 is for year 1980.

    def get_year_from_cftime(cftime_date):
        # Subtract 1 because the date for annual files is when it was SAVED
        return cftime_date.year - 1

    # Check that all desired years are included
    if get_year_from_cftime(ds_in.time.values[0]) > year_1:
        raise RuntimeError(
            f"Requested year_1 is {year_1} but first year in outputs is "
            + f"{get_year_from_cftime(ds_in.time.values[0])}"
        )
    if get_year_from_cftime(ds_in.time.values[-1]) < year_1:
        raise RuntimeError(
            f"Requested year_n is {year_n} but last year in outputs is "
            + f"{get_year_from_cftime(ds_in.time.values[-1])}"
        )

    # Remove years outside range of interest
    ### Include an extra year at the end to finish out final seasons.
    ds_in = utils.safer_timeslice(ds_in, slice(f"{year_1+1}-01-01", f"{year_n+2}-01-01"))

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


def check_constant_vars(
    this_ds, case, ignore_nan, const_growing_seasons=None, verbose=True, throw_error=True
):
    """
    For variables that should stay constant, make sure they are
    """
    if isinstance(case, str):
        const_vars = [case]
    elif isinstance(case, list):
        const_vars = case
    elif isinstance(case, dict):
        const_vars = case["const_vars"]
    else:
        raise TypeError(f"case must be str or dict, not {type(case)}")

    if not const_vars:
        return None

    if const_growing_seasons:
        gs_0 = this_ds.gs.values[0]
        gs_n = this_ds.gs.values[-1]
        if const_growing_seasons.start > gs_0 or const_growing_seasons.stop < gs_n:
            print(
                f"❗ Only checking const_vars over {const_growing_seasons.start}-"
                + f"{const_growing_seasons.stop} (run includes {gs_0}-{gs_n})"
            )
        this_ds = this_ds.sel(gs=const_growing_seasons)

    any_bad = False
    any_bad_before_checking_rx = False
    if throw_error:
        emojus = "❌"
    else:
        emojus = "❗"
    if not isinstance(const_vars, list):
        const_vars = [const_vars]

    for var in const_vars:
        everything_ok = True

        if "gs" in this_ds[var].dims:
            time_coord = "gs"
        elif "time" in this_ds[var].dims:
            time_coord = "time"
        else:
            raise RuntimeError(f"Which of these is the time coordinate? {this_ds[var].dims}")
        i_time_coord = this_ds[var].dims.index(time_coord)

        this_da = this_ds[var]
        ra_sp = np.moveaxis(this_da.copy().values, i_time_coord, 0)
        incl_patches = []
        bad_patches = np.array([])
        str_list = []

        # Read prescription file, if needed
        rx_ds = None
        if isinstance(case, dict):
            if var == "GDDHARV" and "rx_gdds_file" in case:
                rx_ds = import_rx_dates(
                    "gdd", case["rx_gdds_file"], this_ds, set_neg1_to_nan=False
                ).squeeze()

        for time_1 in np.arange(this_ds.dims[time_coord] - 1):
            condn = ~np.isnan(ra_sp[time_1, ...])
            if time_1 > 0:
                condn = np.bitwise_and(condn, np.all(np.isnan(ra_sp[:time_1, ...]), axis=0))
            these_patches = np.where(condn)[0]
            if these_patches.size == 0:
                continue
            these_patches = list(np.where(condn)[0])
            incl_patches += these_patches
            # print(f't1 {t1}: {thesePatches}')

            t1_yr = this_ds[time_coord].values[time_1]
            t1_vals = np.squeeze(this_da.isel({time_coord: time_1, "patch": these_patches}).values)

            for timestep in np.arange(time_1 + 1, this_ds.dims[time_coord]):
                t_yr = this_ds[time_coord].values[timestep]
                t_vals = np.squeeze(
                    this_da.isel({time_coord: timestep, "patch": these_patches}).values
                )
                ok_p = t1_vals == t_vals

                # If allowed, ignore where either t or t1 is NaN. Should only be used for runs where
                # land use varies over time.
                if ignore_nan:
                    ok_p = np.squeeze(np.bitwise_or(ok_p, np.isnan(t1_vals + t_vals)))

                if not np.all(ok_p):
                    any_bad_before_checking_rx = True
                    bad_patches_this_time = list(np.where(np.bitwise_not(ok_p))[0])
                    bad_patches = np.concatenate(
                        (bad_patches, np.array(these_patches)[bad_patches_this_time])
                    )
                    if rx_ds:
                        found_in_rx = np.array([False for x in bad_patches])
                    vary_patches = list(np.array(these_patches)[bad_patches_this_time])
                    vary_lons = this_ds.patches1d_lon.values[bad_patches_this_time]
                    vary_lats = this_ds.patches1d_lat.values[bad_patches_this_time]
                    vary_crops = this_ds.patches1d_itype_veg_str.values[bad_patches_this_time]
                    vary_crops_int = this_ds.patches1d_itype_veg.values[bad_patches_this_time]

                    any_bad_any_crop = False
                    for crop_int in np.unique(vary_crops_int):
                        rx_var = f"gs1_{crop_int}"
                        vary_lons_this_crop = vary_lons[np.where(vary_crops_int == crop_int)]
                        vary_lats_this_crop = vary_lats[np.where(vary_crops_int == crop_int)]
                        these_rx_vals = np.diag(
                            rx_ds[rx_var]
                            .sel(lon=vary_lons_this_crop, lat=vary_lats_this_crop)
                            .values
                        )
                        if len(these_rx_vals) != len(vary_lats_this_crop):
                            raise RuntimeError(
                                f"Expected {len(vary_lats_this_crop)} rx values; got "
                                + f"{len(these_rx_vals)}"
                            )
                        if not np.any(these_rx_vals != -1):
                            continue
                        any_bad_any_crop = True
                        break
                    if not any_bad_any_crop:
                        continue

                    # This bit is pretty inefficient, but I'm not going to optimize it until I
                    # actually need to use it.
                    for i, patch in enumerate(bad_patches_this_time):
                        this_patch = vary_patches[i]
                        this_lon = vary_lons[i]
                        this_lat = vary_lats[i]
                        this_crop = vary_crops[i]
                        this_crop_int = vary_crops_int[i]

                        # If prescribed input had missing value (-1), it's fine for it to vary.
                        if rx_ds:
                            rx_var = f"gs1_{this_crop_int}"
                            if this_lon in rx_ds.lon.values and this_lat in rx_ds.lat.values:
                                rx_vals = rx_ds[rx_var].sel(lon=this_lon, lat=this_lat).values
                                n_unique = len(np.unique(rx_vals))
                                if n_unique == 1:
                                    found_in_rx[i] = True
                                    if rx_vals == -1:
                                        continue
                                elif n_unique > 1:
                                    raise RuntimeError(
                                        f"How does lon {this_lon} lat {this_lat} {this_crop} have "
                                        + f"time-varying {var}?"
                                    )
                            else:
                                raise RuntimeError(
                                    f"lon {this_lon} lat {this_lat} {this_crop} not in rx dataset?"
                                )

                        # Print info (or save to print later)
                        any_bad = True
                        if verbose:
                            this_str = (
                                f"   Patch {this_patch} (lon {this_lon} lat {this_lat}) "
                                + f"{this_crop} ({this_crop_int})"
                            )
                            if rx_ds and not found_in_rx[i]:
                                this_str = this_str.replace("(lon", "* (lon")
                            if not np.isnan(t1_vals[patch]):
                                t1_val_print = int(t1_vals[patch])
                            else:
                                t1_val_print = "NaN"
                            if not np.isnan(t_vals[patch]):
                                t_val_print = int(t_vals[patch])
                            else:
                                t_val_print = "NaN"
                            if var == "SDATES":
                                str_list.append(
                                    f"{this_str}: Sowing {t1_yr} jday {t1_val_print}, {t_yr} "
                                    + f"jday {t_val_print}"
                                )
                            else:
                                str_list.append(
                                    f"{this_str}: {t1_yr} {var} {t1_val_print}, {t_yr} {var} "
                                    + f"{t_val_print}"
                                )
                        else:
                            if everything_ok:
                                print(f"{emojus} CLM output {var} unexpectedly vary over time:")
                                everything_ok = False
                            print(f"{var} timestep {timestep} does not match timestep {time_1}")
                            break
        if verbose and any_bad:
            print(f"{emojus} CLM output {var} unexpectedly vary over time:")
            str_list.sort()
            if rx_ds and np.any(~found_in_rx):
                str_list = [
                    "*: Not found in prescribed input file (maybe minor lon/lat mismatch)"
                ] + str_list
            elif not rx_ds:
                str_list = ["(No rx file checked)"] + str_list
            print("\n".join(str_list))

        # Make sure every patch was checked once (or is all-NaN except possibly final season)
        incl_patches = np.sort(incl_patches)
        if not np.array_equal(incl_patches, np.unique(incl_patches)):
            raise RuntimeError("Patch(es) checked more than once!")
        incl_patches = list(incl_patches)
        incl_patches += list(
            np.where(
                np.all(
                    np.isnan(
                        ra_sp[
                            :-1,
                        ]
                    ),
                    axis=0,
                )
            )[0]
        )
        incl_patches = np.sort(incl_patches)
        if not np.array_equal(incl_patches, np.unique(incl_patches)):
            raise RuntimeError("Patch(es) checked but also all-NaN??")
        if not np.array_equal(incl_patches, np.arange(this_ds.dims["patch"])):
            for patch in np.arange(this_ds.dims["patch"]):
                if patch not in incl_patches:
                    break
            raise RuntimeError(
                f"Not all patches checked! E.g., {patch}: {this_da.isel(patch=patch).values}"
            )

        if not any_bad:
            if any_bad_before_checking_rx:
                print(
                    f"✅ CLM output {var} do not vary through {this_ds.dims[time_coord]} growing "
                    + "seasons of output (except for patch(es) with missing rx)."
                )
            else:
                print(
                    f"✅ CLM output {var} do not vary through {this_ds.dims[time_coord]} growing "
                    + "seasons of output."
                )

    if any_bad and throw_error:
        raise RuntimeError("Stopping due to failed check_constant_vars().")

    bad_patches = np.unique(bad_patches)
    return [int(p) for p in bad_patches]


def check_rx_obeyed(
    vegtype_list, rx_ds, dates_ds, which_ds, output_var, gdd_min=None, verbose=False
):
    """
    Check that prescribed crop calendars were obeyed
    """
    all_ok = 2
    diff_str_list = []
    gdd_tolerance = 1

    if "GDDHARV" in output_var and verbose:
        harvest_reason_da = dates_ds["HARVEST_REASON"]
        unique_harvest_reasons = np.unique(
            harvest_reason_da.values[np.where(~np.isnan(harvest_reason_da.values))]
        )
        pct_harv_at_mature = get_pct_harv_at_mature(harvest_reason_da)
        print(
            f"{which_ds} harvest reasons: {unique_harvest_reasons} ({pct_harv_at_mature}% harv at "
            + "maturity)"
        )

    for vegtype_str in vegtype_list:
        thisveg_patches = np.where(dates_ds.patches1d_itype_veg_str == vegtype_str)[0]
        if thisveg_patches.size == 0:
            continue
        ds_thisveg = dates_ds.isel(patch=thisveg_patches)
        patch_inds_lon_thisveg = ds_thisveg.patches1d_ixy.values.astype(int) - 1
        patch_inds_lat_thisveg = ds_thisveg.patches1d_jxy.values.astype(int) - 1
        patch_lons_thisveg = ds_thisveg.patches1d_lon
        patch_lats_thisveg = ds_thisveg.patches1d_lat

        vegtype_int = utils.vegtype_str2int(vegtype_str)[0]
        rx_da = rx_ds[f"gs1_{vegtype_int}"]
        rx_array = rx_da.values[patch_inds_lat_thisveg, patch_inds_lon_thisveg]
        rx_array = np.expand_dims(rx_array, axis=1)
        sim_array = ds_thisveg[output_var].values
        sim_array_dims = ds_thisveg[output_var].dims

        # Ignore patches without prescribed value
        with np.errstate(invalid="ignore"):
            rx_array[np.where(rx_array < 0)] = np.nan

        # Account for...
        if "GDDHARV" in output_var:
            # ...GDD harvest threshold minimum set in PlantCrop()
            if gdd_min is None:
                gdd_min = DEFAULT_GDD_MIN
                print(
                    f"gdd_min not provided when doing check_rx_obeyed() for {output_var}; using "
                    + f"default {gdd_min}"
                )
            with np.errstate(invalid="ignore"):
                rx_array[(rx_array >= 0) & (rx_array < gdd_min)] = gdd_min

            # ...harvest reason
            # 0: Should never happen in any simulation
            # 1: Harvesting at maturity
            # 2: Harvesting at max season length (mxmat)
            # 3: Crop was incorrectly planted in last time step of Dec. 31
            # 4: Today was supposed to be the planting day, but the previous crop still hasn't been
            #    harvested.
            # 5: Harvest the day before the next sowing date this year.
            # 6: Same as #5.
            # 7: Harvest the day before the next sowing date (today is Dec. 31 and the sowing date
            #    is Jan. 1)
            harvest_reason_da = ds_thisveg["HARVEST_REASON"]
            unique_harvest_reasons = np.unique(
                harvest_reason_da.values[np.where(~np.isnan(harvest_reason_da.values))]
            )
            pct_harv_at_mature = get_pct_harv_at_mature(harvest_reason_da)

        if np.any(sim_array != rx_array):
            diff_array = sim_array - rx_array

            # Allow negative GDDHARV values when harvest occurred because sowing was scheduled for
            # the next day
            if output_var == "GDDHARV_PERHARV":
                diff_array = np.ma.masked_array(
                    diff_array,
                    mask=(diff_array < 0) & (ds_thisveg["HARVEST_REASON_PERHARV"].values == 5),
                )
            elif output_var == "GDDHARV":
                with np.errstate(invalid="ignore"):
                    diff_lt_0 = diff_array < 0
                    harv_reason_5 = ds_thisveg["HARVEST_REASON"].values == 5
                diff_array = np.ma.masked_array(diff_array, mask=diff_lt_0 & harv_reason_5)

            with np.errstate(invalid="ignore"):
                abs_gt_0 = abs(diff_array) > 0
            if np.any(np.abs(diff_array[abs_gt_0]) > 0):
                min_diff, min_lon, min_lat, min_gs, min_rx = get_extreme_info(
                    diff_array,
                    rx_array,
                    np.nanmin,
                    sim_array_dims,
                    dates_ds.gs,
                    patch_lons_thisveg,
                    patch_lats_thisveg,
                )
                max_diff, max_lon, max_lat, max_gs, max_rx = get_extreme_info(
                    diff_array,
                    rx_array,
                    np.nanmax,
                    sim_array_dims,
                    dates_ds.gs,
                    patch_lons_thisveg,
                    patch_lats_thisveg,
                )

                diffs_eg_txt = (
                    f"{vegtype_str} ({vegtype_int}): diffs range {min_diff} (lon {min_lon}, lat "
                    + f"{min_lat}, gs {min_gs}, rx ~{min_rx}) to {max_diff} (lon {max_lon}, lat "
                    + f"{max_lat}, gs {max_gs}, rx ~{max_rx})"
                )
                if "GDDHARV" in output_var:
                    diffs_eg_txt += (
                        f"; harvest reasons: {unique_harvest_reasons} ({pct_harv_at_mature}"
                        + "% harvested at maturity)"
                    )
                if "GDDHARV" in output_var and np.nanmax(abs(diff_array)) <= gdd_tolerance:
                    if all_ok > 0:
                        all_ok = 1
                        diff_str_list.append(f"	  {diffs_eg_txt}")
                else:
                    all_ok = 0
                    if verbose:
                        print(
                            f"❌ {which_ds}: Prescribed {output_var} *not* always obeyed. E.g., "
                            + f"{diffs_eg_txt}"
                        )
                    else:
                        break

    if all_ok == 2:
        print(f"✅ {which_ds}: Prescribed {output_var} always obeyed")
    elif all_ok == 1:
        # print(f"🟨 {which_ds}: Prescribed {output_var} *not* always obeyed, but acceptable:")
        # for x in diff_str_list: print(x)
        print(
            f"🟨 {which_ds}: Prescribed {output_var} *not* always obeyed, but acceptable (diffs <= "
            + f"{gdd_tolerance})"
        )
    elif not verbose:
        print(f"❌ {which_ds}: Prescribed {output_var} *not* always obeyed. E.g., {diffs_eg_txt}")


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
        print(f"✅{msg_txt}{var0} always <= {var1}")
    else:
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


def convert_axis_time2gs(this_ds, verbose=False, my_vars=None, incl_orig=False):
    """
    Convert time*mxharvests axes to growingseason axis
    """
    # How many non-NaN patch-seasons do we expect to have once we're done organizing things?
    n_patch = this_ds.dims["patch"]
    # Because some patches will be planted in the last year but not complete, we have to ignore any
    # finalyear-planted seasons that do complete.
    n_gs = this_ds.dims["time"] - 1
    expected_valid = n_patch * n_gs

    mxharvests = this_ds.dims["mxharvests"]

    if verbose:
        print(
            f"Start: discrepancy of {np.sum(~np.isnan(this_ds.HDATES.values)) - expected_valid} "
            + "patch-seasons"
        )

    # Set all non-positive date values to NaN. These are seasons that were never harvested
    # (or never started): "non-seasons."
    if this_ds.HDATES.dims != ("time", "mxharvests", "patch"):
        raise RuntimeError(
            "This code relies on HDATES dims ('time', 'mxharvests', 'patch'), not "
            + f"{this_ds.HDATES.dims}"
        )
    hdates_ymp = this_ds.HDATES.copy().where(this_ds.HDATES > 0).values
    hdates_pym = np.transpose(hdates_ymp.copy(), (2, 0, 1))
    sdates_ymp = this_ds.SDATES_PERHARV.copy().where(this_ds.SDATES_PERHARV > 0).values
    sdates_pym = np.transpose(sdates_ymp.copy(), (2, 0, 1))
    with np.errstate(invalid="ignore"):
        hdates_pym[hdates_pym <= 0] = np.nan

    # Find years where patch was inactive
    inactive_py = np.transpose(
        np.isnan(this_ds.HDATES).all(dim="mxharvests").values
        & np.isnan(this_ds.SDATES_PERHARV).all(dim="mxharvests").values
    )
    # Find seasons that were planted while the patch was inactive
    with np.errstate(invalid="ignore"):
        sown_inactive_py = inactive_py[:, :-1] & (hdates_pym[:, 1:, 0] < sdates_pym[:, 1:, 0])
    sown_inactive_py = np.concatenate((np.full((n_patch, 1), False), sown_inactive_py), axis=1)

    # "Ignore harvests from seasons sown (a) before this output began or (b) when the crop was
    # inactive"
    with np.errstate(invalid="ignore"):
        first_season_before_first_year_p = hdates_pym[:, 0, 0] < sdates_pym[:, 0, 0]
    first_season_before_first_year_py = np.full(hdates_pym.shape[:-1], fill_value=False)
    first_season_before_first_year_py[:, 0] = first_season_before_first_year_p
    sown_prerun_or_inactive_py = first_season_before_first_year_py | sown_inactive_py
    sown_prerun_or_inactive_pym = np.concatenate(
        (
            np.expand_dims(sown_prerun_or_inactive_py, axis=2),
            np.full((n_patch, n_gs + 1, mxharvests - 1), False),
        ),
        axis=2,
    )
    where_sown_prerun_or_inactive_pym = np.where(sown_prerun_or_inactive_pym)
    hdates_pym[where_sown_prerun_or_inactive_pym] = np.nan
    sdates_pym[where_sown_prerun_or_inactive_pym] = np.nan
    if verbose:
        print(
            "After 'Ignore harvests from before this output began: discrepancy of "
            + f"{np.sum(~np.isnan(hdates_pym)) - expected_valid} patch-seasons'"
        )

    # We need to keep some non-seasons---it's possible that "the yearY growing season" never
    # happened (sowing conditions weren't met), but we still need something there so that we can
    # make an array of dimension Npatch*Ngs. We do this by changing those non-seasons from NaN to
    # -Inf before doing the filtering and reshaping, after which we'll convert them back to NaNs.

    # "In years with no sowing, pretend the first no-harvest is meaningful, unless that was
    # intentionally ignored above."
    sdates_orig_ymp = this_ds.SDATES.copy().values
    sdates_orig_pym = np.transpose(sdates_orig_ymp.copy(), (2, 0, 1))
    hdates_pym2 = hdates_pym.copy()
    sdates_pym2 = sdates_pym.copy()
    with np.errstate(invalid="ignore"):
        sdates_gt_0 = sdates_orig_pym > 0
    nosow_py = np.all(~sdates_gt_0, axis=2)
    nosow_py_1st = nosow_py & np.isnan(hdates_pym[:, :, 0])
    where_nosow_py_1st = np.where(nosow_py_1st)
    hdates_pym2[where_nosow_py_1st[0], where_nosow_py_1st[1], 0] = -np.inf
    sdates_pym2[where_nosow_py_1st[0], where_nosow_py_1st[1], 0] = -np.inf
    for harvest_index in np.arange(mxharvests - 1):
        if harvest_index == 0:
            continue
        elif harvest_index == 1:
            print("Warning: Untested with mxharvests > 2")
        where_nosow_py = np.where(
            nosow_py
            & ~np.any(np.isnan(hdates_pym[:, :, 0:harvest_index]), axis=2)
            & np.isnan(hdates_pym[:, :, harvest_index])
        )
        hdates_pym2[where_nosow_py[0], where_nosow_py[1], harvest_index + 1] = -np.inf
        sdates_pym2[where_nosow_py[0], where_nosow_py[1], harvest_index + 1] = -np.inf

    # "In years with sowing that are followed by inactive years, check whether the last sowing was
    # harvested before the patch was deactivated. If not, pretend the LAST [easier to implement!]
    # no-harvest is meaningful."
    sdates_orig_masked_pym = sdates_orig_pym.copy()
    with np.errstate(invalid="ignore"):
        sdates_le_0 = sdates_orig_masked_pym <= 0
    sdates_orig_masked_pym[np.where(sdates_le_0)] = np.nan
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")
        last_sdate_first_n_gs_py = np.nanmax(sdates_orig_masked_pym[:, :-1, :], axis=2)
        last_hdate_first_n_gs_py = np.nanmax(hdates_pym2[:, :-1, :], axis=2)
    with np.errstate(invalid="ignore"):
        hdate_lt_sdate = last_hdate_first_n_gs_py < last_sdate_first_n_gs_py
    last_sowing_not_harvested_sameyear_first_n_gs_py = hdate_lt_sdate | np.isnan(
        last_hdate_first_n_gs_py
    )
    inactive_last_n_gs_py = inactive_py[:, 1:]
    last_sowing_never_harvested_first_n_gs_py = (
        last_sowing_not_harvested_sameyear_first_n_gs_py & inactive_last_n_gs_py
    )
    last_sowing_never_harvested_py = np.concatenate(
        (last_sowing_never_harvested_first_n_gs_py, np.full((n_patch, 1), False)), axis=1
    )
    last_sowing_never_harvested_pym = np.concatenate(
        (
            np.full((n_patch, n_gs + 1, mxharvests - 1), False),
            np.expand_dims(last_sowing_never_harvested_py, axis=2),
        ),
        axis=2,
    )
    where_last_sowing_never_harvested_pym = last_sowing_never_harvested_pym
    hdates_pym3 = hdates_pym2.copy()
    sdates_pym3 = sdates_pym2.copy()
    hdates_pym3[where_last_sowing_never_harvested_pym] = -np.inf
    sdates_pym3[where_last_sowing_never_harvested_pym] = -np.inf

    # Convert to growingseason axis
    def pym_to_pg(pym_array, quiet=False):
        pg_array = np.reshape(pym_array, (pym_array.shape[0], -1))
        ok_pg = pg_array[~np.isnan(pg_array)]
        if not quiet:
            print(
                f"{ok_pg.size} included; unique N seasons = "
                + f"{np.unique(np.sum(~np.isnan(pg_array), axis=1))}"
            )
        return pg_array

    hdates_pg = pym_to_pg(hdates_pym3.copy(), quiet=~verbose)
    sdates_pg = pym_to_pg(sdates_pym3.copy(), quiet=True)
    if verbose:
        print(
            "After 'In years with no sowing, pretend the first no-harvest is meaningful: "
            + f"discrepancy of {np.sum(~np.isnan(hdates_pg)) - expected_valid} patch-seasons"
        )

    # "Ignore any harvests that were planted in the final year, because some cells will have
    # incomplete growing seasons for the final year."
    with np.errstate(invalid="ignore"):
        hdates_ge_sdates = hdates_pg[:, -mxharvests:] >= sdates_pg[:, -mxharvests:]
    lastyear_complete_season = hdates_ge_sdates | np.isinf(hdates_pg[:, -mxharvests:])

    def ignore_lastyear_complete_season(pg_array, excl, mxharvests):
        tmp_l = pg_array[:, :-mxharvests]
        tmp_r = pg_array[:, -mxharvests:]
        tmp_r[np.where(excl)] = np.nan
        pg_array = np.concatenate((tmp_l, tmp_r), axis=1)
        return pg_array

    hdates_pg2 = ignore_lastyear_complete_season(
        hdates_pg.copy(), lastyear_complete_season, mxharvests
    )
    sdates_pg2 = ignore_lastyear_complete_season(
        sdates_pg.copy(), lastyear_complete_season, mxharvests
    )
    is_valid = ~np.isnan(hdates_pg2)
    is_fake = np.isneginf(hdates_pg2)
    is_fake = np.reshape(is_fake[is_valid], (this_ds.dims["patch"], n_gs))
    discrepancy = np.sum(is_valid) - expected_valid
    unique_n_seasons = np.unique(np.sum(is_valid, axis=1))
    if verbose:
        print(
            "After 'Ignore any harvests that were planted in the final year, because other cells "
            + "will have incomplete growing seasons for the final year': discrepancy of "
            + f"{discrepancy} patch-seasons"
        )
        if "pandas" in sys.modules:
            bincount = np.bincount(np.sum(is_valid, axis=1))
            bincount = bincount[bincount > 0]
            dataframe = pd.DataFrame({"Ngs": unique_n_seasons, "Count": bincount})
            print(dataframe)
        else:
            print(f"unique N seasons = {unique_n_seasons}")
        print(" ")

    # Create Dataset with time axis as "gs" (growing season) instead of what CLM puts out
    if discrepancy == 0:
        this_ds_gs = set_up_ds_with_gs_axis(this_ds)
        for var in this_ds.data_vars:
            if this_ds[var].dims != ("time", "mxharvests", "patch") or (
                my_vars and var not in my_vars
            ):
                continue

            # Set invalid values to NaN
            da_yhp = this_ds[var].copy()
            da_yhp = da_yhp.where(~np.isneginf(da_yhp))

            # Remove the nans and reshape to patches*growingseasons
            da_pyh = da_yhp.transpose("patch", "time", "mxharvests")
            ar_pg = np.reshape(da_pyh.values, (this_ds.dims["patch"], -1))
            ar_valid_pg = np.reshape(ar_pg[is_valid], (this_ds.dims["patch"], n_gs))
            # Change -infs to nans
            ar_valid_pg[is_fake] = np.nan
            # Save as DataArray to new Dataset, stripping _PERHARV from variable name
            newname = var.replace("_PERHARV", "")
            if newname in this_ds_gs:
                raise RuntimeError(f"{newname} already in dataset!")
            da_pg = xr.DataArray(
                data=ar_valid_pg,
                coords=[this_ds_gs.coords["patch"], this_ds_gs.coords["gs"]],
                name=newname,
                attrs=da_yhp.attrs,
            )
            this_ds_gs[newname] = da_pg
            this_ds_gs[newname].attrs["units"] = this_ds[var].attrs["units"]
    else:
        # Print details about example bad patch(es)
        if min(unique_n_seasons) < n_gs:
            print(f"Too few seasons (min {min(unique_n_seasons)} < {n_gs})")
            patch_index = np.where(np.sum(~np.isnan(hdates_pg2), axis=1) == min(unique_n_seasons))[
                0
            ][0]
            print_onepatch_wrong_n_gs(
                patch_index,
                this_ds,
                sdates_ymp,
                hdates_ymp,
                sdates_pym,
                hdates_pym,
                sdates_pym2,
                hdates_pym2,
                sdates_pym3,
                hdates_pym3,
                sdates_pg,
                hdates_pg,
                sdates_pg2,
                hdates_pg2,
            )
        if max(unique_n_seasons) > n_gs:
            print(f"Too many seasons (max {max(unique_n_seasons)} > {n_gs})")
            patch_index = np.where(np.sum(~np.isnan(hdates_pg2), axis=1) == max(unique_n_seasons))[
                0
            ][0]
            print_onepatch_wrong_n_gs(
                patch_index,
                this_ds,
                sdates_ymp,
                hdates_ymp,
                sdates_pym,
                hdates_pym,
                sdates_pym2,
                hdates_pym2,
                sdates_pym3,
                hdates_pym3,
                sdates_pg,
                hdates_pg,
                sdates_pg2,
                hdates_pg2,
            )
        raise RuntimeError(
            "Can't convert time*mxharvests axes to growingseason axis: discrepancy of "
            + f"{discrepancy} patch-seasons"
        )

    # Preserve units
    for var_1 in this_ds_gs:
        var_0 = var_1
        if var_0 not in this_ds:
            var_0 += "_PERHARV"
        if var_0 not in this_ds:
            continue
        if "units" in this_ds[var_0].attrs:
            this_ds_gs[var_1].attrs["units"] = this_ds[var_0].attrs["units"]

    if incl_orig:
        return this_ds_gs, this_ds
    return this_ds_gs


def get_extreme_info(diff_array, rx_array, mxn, dims, gs_da, patches1d_lon, patches1d_lat):
    """
    Get information about extreme gridcells (for debugging)
    """
    if mxn == np.min:  # pylint: disable=comparison-with-callable
        diff_array = np.ma.masked_array(diff_array, mask=np.abs(diff_array) == 0)
    themxn = mxn(diff_array)

    # Find the first patch-gs that has the mxn value
    matching_indices = np.where(diff_array == themxn)
    first_indices = [x[0] for x in matching_indices]

    # Get the lon, lat, and growing season of that patch-gs
    patch_index = first_indices[dims.index("patch")]
    this_lon = patches1d_lon.values[patch_index]
    this_lat = patches1d_lat.values[patch_index]
    season_index = first_indices[dims.index("gs")]
    this_gs = gs_da.values[season_index]

    # Get the prescribed value for this patch-gs
    this_rx = rx_array[patch_index][0]

    return round(themxn, 3), round(this_lon, 3), round(this_lat, 3), this_gs, round(this_rx)


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


def get_pct_harv_at_mature(harvest_reason_da):
    """
    Get percentage of harvests that happened at maturity
    """
    n_harv_at_mature = len(np.where(harvest_reason_da.values == 1)[0])
    with np.errstate(invalid="ignore"):
        harv_reason_gt_0 = harvest_reason_da.values > 0
    n_harv = len(np.where(harv_reason_gt_0)[0])
    if n_harv == 0:
        return np.nan
    pct_harv_at_mature = n_harv_at_mature / n_harv * 100
    pct_harv_at_mature = np.format_float_positional(
        pct_harv_at_mature, precision=2, unique=False, fractional=False, trim="k"
    )  # Round to 2 significant digits
    return pct_harv_at_mature


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

    this_ds = utils.import_ds(date_infile, myVars=date_varlist)

    did_warn = False
    for var in this_ds:
        v_new = var.replace(var_prefix, "gs")
        this_ds = this_ds.rename({var: v_new})

        # Set -1 prescribed GDD values to NaN. Only warn the first time.
        if (
            set_neg1_to_nan
            and var_prefix == "gdd"
            and v_new != var
            and np.any(this_ds[v_new].values < 0)
        ):
            if np.any((this_ds[v_new].values < 0) & (this_ds[v_new].values != -1)):
                raise RuntimeError(f"Unexpected negative value in {var}")
            if not did_warn:
                print("Setting -1 rx GDD values to NaN")
                did_warn = True
            this_ds[v_new] = this_ds[v_new].where(this_ds[v_new] != -1)

    return this_ds


def import_output(
    filename,
    my_vars,
    year_1=None,
    year_n=None,
    my_vegtypes=utils.define_mgdcrop_list(),
    sdates_rx_ds=None,
    gdds_rx_ds=None,
    verbose=False,
):
    """
    Import CLM output
    """
    # Import
    this_ds = utils.import_ds(filename, myVars=my_vars, myVegtypes=my_vegtypes)

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
    #
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

    def check_no_negative(this_ds_in, varlist_no_negative, which_file, verbose=False):
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

    def check_no_zeros(this_ds, varlist_no_zero, which_file):
        for var in this_ds:
            if not any(x in var for x in varlist_no_zero):
                continue
            if np.any(this_ds[var].values == 0):
                print(f"WARNING: Unexpected zero(s) in {var} ({which_file})")
            elif verbose:
                print(f"No zero value(s) in {var} ({which_file})")

    # Check for no zero values where there shouldn't be
    varlist_no_zero = ["DATE", "YEAR"]
    check_no_zeros(this_ds, varlist_no_zero, "original file")

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
    this_ds_gs = check_no_negative(this_ds_gs, varlist_no_negative, "new file", verbose=verbose)

    # Check for no zero values where there shouldn't be
    varlist_no_zero = ["REASON", "DATE"]
    check_no_zeros(this_ds_gs, varlist_no_zero, "new file")

    # Check that e.g., GDDACCUM <= HUI
    for var_list in [["GDDACCUM", "HUI"], ["SYEARS", "HYEARS"]]:
        if all(v in this_ds_gs for v in var_list):
            check_v0_le_v1(this_ds_gs, var_list, both_nan_ok=True, throw_error=True)

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
    this_ds_gs = this_ds_gs.assign_coords(
        {"cftime": this_ds["time_bounds"].isel({"hist_interval": 0})}
    )
    this_ds_gs = this_ds_gs.assign_coords({"time": [t.year for t in this_ds_gs["cftime"].values]})

    # Get number of harvests
    this_ds_gs["NHARVESTS"] = (this_ds_gs["GDDHARV_PERHARV"] > 0).sum(dim="mxharvests")
    # Get number of harvests that would be missed if only seeing max 1 per calendar year
    if np.any(this_ds_gs["NHARVESTS"] > 2):
        raise RuntimeError("How to get NHARVEST_DISCREP for NHARVESTS > 2?")
    this_ds_gs["NHARVEST_DISCREP"] = (this_ds_gs["NHARVESTS"] == 2).astype(int)

    return this_ds_gs


def print_onepatch_wrong_n_gs(
    patch_index,
    this_ds_orig,
    sdates_ymp,
    hdates_ymp,
    sdates_pym,
    hdates_pym,
    sdates_pym2,
    hdates_pym2,
    sdates_pym3,
    hdates_pym3,
    sdates_pg,
    hdates_pg,
    sdates_pg2,
    hdates_pg2,
):
    """
    Print information about a patch (for debugging)
    """

    print(
        f"patch {patch_index}: {this_ds_orig.patches1d_itype_veg_str.values[patch_index]}, lon "
        f"{this_ds_orig.patches1d_lon.values[patch_index]} lat "
        f"{this_ds_orig.patches1d_lat.values[patch_index]}"
    )

    print("Original SDATES (per sowing):")
    print(this_ds_orig.SDATES.values[:, :, patch_index])

    print("Original HDATES (per harvest):")
    print(this_ds_orig.HDATES.values[:, :, patch_index])

    if "pandas" in sys.modules:

        def print_pandas_ymp(msg, cols, arrs_tuple):
            print(f"{msg} ({np.sum(~np.isnan(arrs_tuple[0]))})")
            mxharvests = arrs_tuple[0].shape[1]
            arrs_list2 = []
            cols2 = []
            for harvest_index in np.arange(mxharvests):
                for i, array in enumerate(arrs_tuple):
                    arrs_list2.append(array[:, harvest_index])
                    cols2.append(cols[i] + str(harvest_index))
            arrs_tuple2 = tuple(arrs_list2)
            dataframe = pd.DataFrame(np.stack(arrs_tuple2, axis=1))
            dataframe.columns = cols2
            print(dataframe)

        print_pandas_ymp(
            "Original",
            ["sdate", "hdate"],
            (
                this_ds_orig.SDATES_PERHARV.values[:, :, patch_index],
                this_ds_orig.HDATES.values[:, :, patch_index],
            ),
        )

        print_pandas_ymp(
            "Masked",
            ["sdate", "hdate"],
            (sdates_ymp[:, :, patch_index], hdates_ymp[:, :, patch_index]),
        )

        print_pandas_ymp(
            'After "Ignore harvests from before this output began"',
            ["sdate", "hdate"],
            (
                np.transpose(sdates_pym, (1, 2, 0))[:, :, patch_index],
                np.transpose(hdates_pym, (1, 2, 0))[:, :, patch_index],
            ),
        )

        print_pandas_ymp(
            'After "In years with no sowing, pretend the first no-harvest is meaningful"',
            ["sdate", "hdate"],
            (
                np.transpose(sdates_pym2, (1, 2, 0))[:, :, patch_index],
                np.transpose(hdates_pym2, (1, 2, 0))[:, :, patch_index],
            ),
        )

        print_pandas_ymp(
            (
                'After "In years with sowing that are followed by inactive years, check whether the'
                " last sowing was harvested before the patch was deactivated. If not, pretend the"
                ' LAST no-harvest is meaningful."'
            ),
            ["sdate", "hdate"],
            (
                np.transpose(sdates_pym3, (1, 2, 0))[:, :, patch_index],
                np.transpose(hdates_pym3, (1, 2, 0))[:, :, patch_index],
            ),
        )

        def print_pandas_pg(msg, cols, arrs_tuple):
            print(f"{msg} ({np.sum(~np.isnan(arrs_tuple[0]))})")
            arrs_list = list(arrs_tuple)
            for i, array in enumerate(arrs_tuple):
                arrs_list[i] = np.reshape(array, (-1))
            arrs_tuple2 = tuple(arrs_list)
            dataframe = pd.DataFrame(np.stack(arrs_tuple2, axis=1))
            dataframe.columns = cols
            print(dataframe)

        print_pandas_pg(
            "Same, but converted to gs axis",
            ["sdate", "hdate"],
            (sdates_pg[patch_index, :], hdates_pg[patch_index, :]),
        )

        print_pandas_pg(
            (
                'After "Ignore any harvests that were planted in the final year, because some cells'
                ' will have incomplete growing seasons for the final year"'
            ),
            ["sdate", "hdate"],
            (sdates_pg2[patch_index, :], hdates_pg2[patch_index, :]),
        )
    else:
        print("Couldn't import pandas, so not displaying example bad patch ORIGINAL.")

        def print_nopandas(array_1, array_2, msg):
            print(msg)
            if array_1.ndim == 1:
                # I don't know why these aren't side-by-side!
                print(np.stack((array_1, array_2), axis=1))
            else:
                print(np.concatenate((array_1, array_2), axis=1))

        print_nopandas(sdates_ymp[:, :, patch_index], hdates_ymp[:, :, patch_index], "Masked:")

        print_nopandas(
            np.transpose(sdates_pym, (1, 2, 0))[:, :, patch_index],
            np.transpose(hdates_pym, (1, 2, 0))[:, :, patch_index],
            'After "Ignore harvests from before this output began"',
        )

        print_nopandas(
            np.transpose(sdates_pym2, (1, 2, 0))[:, :, patch_index],
            np.transpose(hdates_pym2, (1, 2, 0))[:, :, patch_index],
            'After "In years with no sowing, pretend the first no-harvest is meaningful"',
        )

        print_nopandas(
            np.transpose(sdates_pym3, (1, 2, 0))[:, :, patch_index],
            np.transpose(hdates_pym3, (1, 2, 0))[:, :, patch_index],
            (
                'After "In years with sowing that are followed by inactive years, check whether the'
                " last sowing was harvested before the patch was deactivated. If not, pretend the"
                ' LAST [easier to implement!] no-harvest is meaningful."'
            ),
        )

        print_nopandas(
            sdates_pg[patch_index, :], hdates_pg[patch_index, :], "Same, but converted to gs axis"
        )

        print_nopandas(
            sdates_pg2[patch_index, :],
            hdates_pg2[patch_index, :],
            (
                'After "Ignore any harvests that were planted in the final year, because some cells'
                ' will have incomplete growing seasons for the final year"'
            ),
        )

    print("\n\n")


def set_up_ds_with_gs_axis(ds_in):
    """
    Set up empty Dataset with time axis as "gs" (growing season) instead of what CLM puts out.

    Includes all the same variables as the input dataset, minus any that had dimensions mxsowings or
    mxharvests.
    """
    # Get the data variables to include in the new dataset
    data_vars = {}
    for var in ds_in.data_vars:
        if not any(x in ["mxsowings", "mxharvests"] for x in ds_in[var].dims):
            data_vars[var] = ds_in[var]
    # Set up the new dataset
    gs_years = [t.year - 1 for t in ds_in.time.values[:-1]]
    coords = ds_in.coords
    coords["gs"] = gs_years
    ds_out = xr.Dataset(data_vars=data_vars, coords=coords, attrs=ds_in.attrs)
    return ds_out
