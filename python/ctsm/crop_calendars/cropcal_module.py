# Import the CTSM Python utilities
import utils

import numpy as np
import xarray as xr
from scipy import stats, signal
import warnings
import cftime
import pandas as pd
import os
import glob

try:
    import cartopy.crs as ccrs
except:
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


# After importing a file, restrict it to years of interest.
def check_and_trim_years(y1, yN, ds_in):
    ### In annual outputs, file with name Y is actually results from year Y-1.
    ### Note that time values refer to when it was SAVED. So 1981-01-01 is for year 1980.

    def get_year_from_cftime(cftime_date):
        # Subtract 1 because the date for annual files is when it was SAVED
        return cftime_date.year - 1

    # Check that all desired years are included
    if get_year_from_cftime(ds_in.time.values[0]) > y1:
        raise RuntimeError(
            f"Requested y1 is {y1} but first year in outputs is {get_year_from_cftime(ds_in.time.values[0])}"
        )
    elif get_year_from_cftime(ds_in.time.values[-1]) < y1:
        raise RuntimeError(
            f"Requested yN is {yN} but last year in outputs is {get_year_from_cftime(ds_in.time.values[-1])}"
        )

    # Remove years outside range of interest
    ### Include an extra year at the end to finish out final seasons.
    ds_in = utils.safer_timeslice(ds_in, slice(f"{y1+1}-01-01", f"{yN+2}-01-01"))

    # Make sure you have the expected number of timesteps (including extra year)
    Nyears_expected = yN - y1 + 2
    if ds_in.dims["time"] != Nyears_expected:
        raise RuntimeError(
            f"Expected {Nyears_expected} timesteps in output but got {ds_in.dims['time']}"
        )

    return ds_in


def check_constant_vars(
    this_ds, case, ignore_nan, constantGSs=None, verbose=True, throw_error=True
):
    if isinstance(case, str):
        constantVars = [case]
    elif isinstance(case, list):
        constantVars = case
    elif isinstance(case, dict):
        constantVars = case["constantVars"]
    else:
        raise TypeError(f"case must be str or dict, not {type(case)}")

    if not constantVars:
        return None

    if constantGSs:
        gs0 = this_ds.gs.values[0]
        gsN = this_ds.gs.values[-1]
        if constantGSs.start > gs0 or constantGSs.stop < gsN:
            print(
                f"❗ Only checking constantVars over {constantGSs.start}-{constantGSs.stop} (run includes {gs0}-{gsN})"
            )
        this_ds = this_ds.sel(gs=constantGSs)

    any_bad = False
    any_bad_before_checking_rx = False
    if throw_error:
        emojus = "❌"
    else:
        emojus = "❗"
    if not isinstance(constantVars, list):
        constantVars = [constantVars]

    for v in constantVars:
        ok = True

        if "gs" in this_ds[v].dims:
            time_coord = "gs"
        elif "time" in this_ds[v].dims:
            time_coord = "time"
        else:
            raise RuntimeError(f"Which of these is the time coordinate? {this_ds[v].dims}")
        i_time_coord = this_ds[v].dims.index(time_coord)

        this_da = this_ds[v]
        ra_sp = np.moveaxis(this_da.copy().values, i_time_coord, 0)
        incl_patches = []
        bad_patches = np.array([])
        strList = []

        # Read prescription file, if needed
        rx_ds = None
        if isinstance(case, dict):
            if v == "GDDHARV" and "rx_gdds_file" in case:
                rx_ds = import_rx_dates(
                    "gdd", case["rx_gdds_file"], this_ds, set_neg1_to_nan=False
                ).squeeze()

        for t1 in np.arange(this_ds.dims[time_coord] - 1):
            condn = ~np.isnan(ra_sp[t1, ...])
            if t1 > 0:
                condn = np.bitwise_and(condn, np.all(np.isnan(ra_sp[:t1, ...]), axis=0))
            thesePatches = np.where(condn)[0]
            if thesePatches.size == 0:
                continue
            thesePatches = list(np.where(condn)[0])
            incl_patches += thesePatches
            # print(f't1 {t1}: {thesePatches}')

            t1_yr = this_ds[time_coord].values[t1]
            t1_vals = np.squeeze(this_da.isel({time_coord: t1, "patch": thesePatches}).values)

            for t in np.arange(t1 + 1, this_ds.dims[time_coord]):
                t_yr = this_ds[time_coord].values[t]
                t_vals = np.squeeze(this_da.isel({time_coord: t, "patch": thesePatches}).values)
                ok_p = t1_vals == t_vals

                # If allowed, ignore where either t or t1 is NaN. Should only be used for runs where land use varies over time.
                if ignore_nan:
                    ok_p = np.squeeze(np.bitwise_or(ok_p, np.isnan(t1_vals + t_vals)))

                if not np.all(ok_p):
                    any_bad_before_checking_rx = True
                    bad_patches_thisT = list(np.where(np.bitwise_not(ok_p))[0])
                    bad_patches = np.concatenate(
                        (bad_patches, np.array(thesePatches)[bad_patches_thisT])
                    )
                    if rx_ds:
                        found_in_rx = np.array([False for x in bad_patches])
                    varyPatches = list(np.array(thesePatches)[bad_patches_thisT])
                    varyLons = this_ds.patches1d_lon.values[bad_patches_thisT]
                    varyLats = this_ds.patches1d_lat.values[bad_patches_thisT]
                    varyCrops = this_ds.patches1d_itype_veg_str.values[bad_patches_thisT]
                    varyCrops_int = this_ds.patches1d_itype_veg.values[bad_patches_thisT]

                    any_bad_anyCrop = False
                    for c in np.unique(varyCrops_int):
                        rx_var = f"gs1_{c}"
                        varyLons_thisCrop = varyLons[np.where(varyCrops_int == c)]
                        varyLats_thisCrop = varyLats[np.where(varyCrops_int == c)]
                        theseRxVals = np.diag(
                            rx_ds[rx_var].sel(lon=varyLons_thisCrop, lat=varyLats_thisCrop).values
                        )
                        if len(theseRxVals) != len(varyLats_thisCrop):
                            raise RuntimeError(
                                f"Expected {len(varyLats_thisCrop)} rx values; got {len(theseRxVals)}"
                            )
                        if not np.any(theseRxVals != -1):
                            continue
                        any_bad_anyCrop = True
                        break
                    if not any_bad_anyCrop:
                        continue

                    # This bit is pretty inefficient, but I'm not going to optimize it until I actually need to use it.
                    for i, p in enumerate(bad_patches_thisT):
                        thisPatch = varyPatches[i]
                        thisLon = varyLons[i]
                        thisLat = varyLats[i]
                        thisCrop = varyCrops[i]
                        thisCrop_int = varyCrops_int[i]

                        # If prescribed input had missing value (-1), it's fine for it to vary.
                        if rx_ds:
                            rx_var = f"gs1_{thisCrop_int}"
                            if thisLon in rx_ds.lon.values and thisLat in rx_ds.lat.values:
                                rx = rx_ds[rx_var].sel(lon=thisLon, lat=thisLat).values
                                Nunique = len(np.unique(rx))
                                if Nunique == 1:
                                    found_in_rx[i] = True
                                    if rx == -1:
                                        continue
                                elif Nunique > 1:
                                    raise RuntimeError(
                                        f"How does lon {thisLon} lat {thisLat} {thisCrop} have time-varying {v}?"
                                    )
                            else:
                                raise RuntimeError(
                                    "lon {thisLon} lat {thisLat} {thisCrop} not in rx dataset?"
                                )

                        # Print info (or save to print later)
                        any_bad = True
                        if verbose:
                            thisStr = f"   Patch {thisPatch} (lon {thisLon} lat {thisLat}) {thisCrop} ({thisCrop_int})"
                            if rx_ds and not found_in_rx[i]:
                                thisStr = thisStr.replace("(lon", "* (lon")
                            if not np.isnan(t1_vals[p]):
                                t1_val_print = int(t1_vals[p])
                            else:
                                t1_val_print = "NaN"
                            if not np.isnan(t_vals[p]):
                                t_val_print = int(t_vals[p])
                            else:
                                t_val_print = "NaN"
                            if v == "SDATES":
                                strList.append(
                                    f"{thisStr}: Sowing {t1_yr} jday {t1_val_print}, {t_yr} jday {t_val_print}"
                                )
                            else:
                                strList.append(
                                    f"{thisStr}: {t1_yr} {v} {t1_val_print}, {t_yr} {v} {t_val_print}"
                                )
                        else:
                            if ok:
                                print(f"{emojus} CLM output {v} unexpectedly vary over time:")
                                ok = False
                            print(f"{v} timestep {t} does not match timestep {t1}")
                            break
        if verbose and any_bad:
            print(f"{emojus} CLM output {v} unexpectedly vary over time:")
            strList.sort()
            if rx_ds and np.any(~found_in_rx):
                strList = [
                    "*: Not found in prescribed input file (maybe minor lon/lat mismatch)"
                ] + strList
            elif not rx_ds:
                strList = ["(No rx file checked)"] + strList
            print("\n".join(strList))

        # Make sure every patch was checked once (or is all-NaN except possibly final season)
        incl_patches = np.sort(incl_patches)
        if not np.array_equal(incl_patches, np.unique(incl_patches)):
            raise RuntimeError("Patch(es) checked more than once!")
        incl_patches = list(incl_patches)
        incl_patches += list(np.where(np.all(np.isnan(ra_sp[:-1,]), axis=0))[0])
        incl_patches = np.sort(incl_patches)
        if not np.array_equal(incl_patches, np.unique(incl_patches)):
            raise RuntimeError("Patch(es) checked but also all-NaN??")
        if not np.array_equal(incl_patches, np.arange(this_ds.dims["patch"])):
            for p in np.arange(this_ds.dims["patch"]):
                if p not in incl_patches:
                    break
            raise RuntimeError(
                f"Not all patches checked! E.g., {p}: {this_da.isel(patch=p).values}"
            )

        if not any_bad:
            if any_bad_before_checking_rx:
                print(
                    f"✅ CLM output {v} do not vary through {this_ds.dims[time_coord]} growing seasons of output (except for patch(es) with missing rx)."
                )
            else:
                print(
                    f"✅ CLM output {v} do not vary through {this_ds.dims[time_coord]} growing seasons of output."
                )

    if any_bad and throw_error:
        raise RuntimeError("Stopping due to failed check_constant_vars().")

    bad_patches = np.unique(bad_patches)
    return [int(p) for p in bad_patches]


def check_rx_obeyed(
    vegtype_list, rx_ds, dates_ds, which_ds, output_var, gdd_min=None, verbose=False
):
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
            f"{which_ds} harvest reasons: {unique_harvest_reasons} ({pct_harv_at_mature}% harv at maturity)"
        )

    for vegtype_str in vegtype_list:
        thisVeg_patches = np.where(dates_ds.patches1d_itype_veg_str == vegtype_str)[0]
        if thisVeg_patches.size == 0:
            continue
        ds_thisVeg = dates_ds.isel(patch=thisVeg_patches)
        patch_inds_lon_thisVeg = ds_thisVeg.patches1d_ixy.values.astype(int) - 1
        patch_inds_lat_thisVeg = ds_thisVeg.patches1d_jxy.values.astype(int) - 1
        patch_lons_thisVeg = ds_thisVeg.patches1d_lon
        patch_lats_thisVeg = ds_thisVeg.patches1d_lat

        vegtype_int = utils.vegtype_str2int(vegtype_str)[0]
        rx_da = rx_ds[f"gs1_{vegtype_int}"]
        rx_array = rx_da.values[patch_inds_lat_thisVeg, patch_inds_lon_thisVeg]
        rx_array = np.expand_dims(rx_array, axis=1)
        sim_array = ds_thisVeg[output_var].values
        sim_array_dims = ds_thisVeg[output_var].dims

        # Ignore patches without prescribed value
        with np.errstate(invalid="ignore"):
            rx_array[np.where(rx_array < 0)] = np.nan

        # Account for...
        if "GDDHARV" in output_var:
            # ...GDD harvest threshold minimum set in PlantCrop()
            if gdd_min == None:
                gdd_min = default_gdd_min()
                print(
                    f"gdd_min not provided when doing check_rx_obeyed() for {output_var}; using default {gdd_min}"
                )
            with np.errstate(invalid="ignore"):
                rx_array[(rx_array >= 0) & (rx_array < gdd_min)] = gdd_min

            # ...harvest reason
            # 0: Should never happen in any simulation
            # 1: Harvesting at maturity
            # 2: Harvesting at max season length (mxmat)
            # 3: Crop was incorrectly planted in last time step of Dec. 31
            # 4: Today was supposed to be the planting day, but the previous crop still hasn't been harvested.
            # 5: Harvest the day before the next sowing date this year.
            # 6: Same as #5.
            # 7: Harvest the day before the next sowing date (today is Dec. 31 and the sowing date is Jan. 1)
            harvest_reason_da = ds_thisVeg["HARVEST_REASON"]
            unique_harvest_reasons = np.unique(
                harvest_reason_da.values[np.where(~np.isnan(harvest_reason_da.values))]
            )
            pct_harv_at_mature = get_pct_harv_at_mature(harvest_reason_da)

        if np.any(sim_array != rx_array):
            diff_array = sim_array - rx_array

            # Allow negative GDDHARV values when harvest occurred because sowing was scheduled for the next day
            if output_var == "GDDHARV_PERHARV":
                diff_array = np.ma.masked_array(
                    diff_array,
                    mask=(diff_array < 0) & (ds_thisVeg["HARVEST_REASON_PERHARV"].values == 5),
                )
            elif output_var == "GDDHARV":
                with np.errstate(invalid="ignore"):
                    diff_lt_0 = diff_array < 0
                    harv_reason_5 = ds_thisVeg["HARVEST_REASON"].values == 5
                diff_array = np.ma.masked_array(diff_array, mask=diff_lt_0 & harv_reason_5)

            with np.errstate(invalid="ignore"):
                abs_gt_0 = abs(diff_array) > 0
            if np.any(np.abs(diff_array[abs_gt_0]) > 0):
                min_diff, minLon, minLat, minGS, minRx = get_extreme_info(
                    diff_array,
                    rx_array,
                    np.nanmin,
                    sim_array_dims,
                    dates_ds.gs,
                    patch_lons_thisVeg,
                    patch_lats_thisVeg,
                )
                max_diff, maxLon, maxLat, maxGS, maxRx = get_extreme_info(
                    diff_array,
                    rx_array,
                    np.nanmax,
                    sim_array_dims,
                    dates_ds.gs,
                    patch_lons_thisVeg,
                    patch_lats_thisVeg,
                )

                diffs_eg_txt = f"{vegtype_str} ({vegtype_int}): diffs range {min_diff} (lon {minLon}, lat {minLat}, gs {minGS}, rx ~{minRx}) to {max_diff} (lon {maxLon}, lat {maxLat}, gs {maxGS}, rx ~{maxRx})"
                if "GDDHARV" in output_var:
                    diffs_eg_txt += f"; harvest reasons: {unique_harvest_reasons} ({pct_harv_at_mature}% harvested at maturity)"
                if "GDDHARV" in output_var and np.nanmax(abs(diff_array)) <= gdd_tolerance:
                    if all_ok > 0:
                        all_ok = 1
                        diff_str_list.append(f"	  {diffs_eg_txt}")
                else:
                    all_ok = 0
                    if verbose:
                        print(
                            f"❌ {which_ds}: Prescribed {output_var} *not* always obeyed. E.g., {diffs_eg_txt}"
                        )
                    else:
                        break

    if all_ok == 2:
        print(f"✅ {which_ds}: Prescribed {output_var} always obeyed")
    elif all_ok == 1:
        # print(f"🟨 {which_ds}: Prescribed {output_var} *not* always obeyed, but acceptable:")
        # for x in diff_str_list: print(x)
        print(
            f"🟨 {which_ds}: Prescribed {output_var} *not* always obeyed, but acceptable (diffs <= {gdd_tolerance})"
        )
    elif not verbose:
        print(f"❌ {which_ds}: Prescribed {output_var} *not* always obeyed. E.g., {diffs_eg_txt}")


# Make sure that, e.g., GDDACCUM_PERHARV is always <= HUI_PERHARV
def check_v0_le_v1(this_ds, vars, msg_txt=" ", both_nan_ok=False, throw_error=False):
    v0 = vars[0]
    v1 = vars[1]
    gdd_lt_hui = this_ds[v0] <= this_ds[v1]
    if both_nan_ok:
        gdd_lt_hui = gdd_lt_hui | (np.isnan(this_ds[v0]) & np.isnan(this_ds[v1]))
    if np.all(gdd_lt_hui):
        print(f"✅{msg_txt}{v0} always <= {v1}")
    else:
        msg = f"❌{msg_txt}{v0} *not* always <= {v1}"
        gdd_lt_hui_vals = gdd_lt_hui.values
        p = np.where(~gdd_lt_hui_vals)[0][0]
        msg = (
            msg
            + f"\ne.g., patch {p}: {this_ds.patches1d_itype_veg_str.values[p]}, lon {this_ds.patches1d_lon.values[p]} lat {this_ds.patches1d_lat.values[p]}:"
        )
        msg = msg + f"\n{this_ds[v0].values[p,:]}"
        msg = msg + f"\n{this_ds[v1].values[p,:]}"
        if throw_error:
            print(msg)
        else:
            raise RuntimeError(msg)


# Convert time*mxharvests axes to growingseason axis
def convert_axis_time2gs(this_ds, verbose=False, myVars=None, incl_orig=False):
    # For backwards compatibility.
    if "SDATES_PERHARV" not in this_ds:
        return convert_axis_time2gs_old(this_ds, myVars=myVars)
    # Otherwise...

    # How many non-NaN patch-seasons do we expect to have once we're done organizing things?
    Npatch = this_ds.dims["patch"]
    # Because some patches will be planted in the last year but not complete, we have to ignore any finalyear-planted seasons that do complete.
    Ngs = this_ds.dims["time"] - 1
    expected_valid = Npatch * Ngs

    mxharvests = this_ds.dims["mxharvests"]

    if verbose:
        print(
            f"Start: discrepancy of {np.sum(~np.isnan(this_ds.HDATES.values)) - expected_valid} patch-seasons"
        )

    # Set all non-positive date values to NaN. These are seasons that were never harvested (or never started): "non-seasons."
    if this_ds.HDATES.dims != ("time", "mxharvests", "patch"):
        raise RuntimeError(
            f"This code relies on HDATES dims ('time', 'mxharvests', 'patch'), not {this_ds.HDATES.dims}"
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
    sown_inactive_py = np.concatenate((np.full((Npatch, 1), False), sown_inactive_py), axis=1)

    # "Ignore harvests from seasons sown (a) before this output began or (b) when the crop was inactive"
    with np.errstate(invalid="ignore"):
        first_season_before_first_year_p = hdates_pym[:, 0, 0] < sdates_pym[:, 0, 0]
    first_season_before_first_year_py = np.full(hdates_pym.shape[:-1], fill_value=False)
    first_season_before_first_year_py[:, 0] = first_season_before_first_year_p
    sown_prerun_or_inactive_py = first_season_before_first_year_py | sown_inactive_py
    sown_prerun_or_inactive_pym = np.concatenate(
        (
            np.expand_dims(sown_prerun_or_inactive_py, axis=2),
            np.full((Npatch, Ngs + 1, mxharvests - 1), False),
        ),
        axis=2,
    )
    where_sown_prerun_or_inactive_pym = np.where(sown_prerun_or_inactive_pym)
    hdates_pym[where_sown_prerun_or_inactive_pym] = np.nan
    sdates_pym[where_sown_prerun_or_inactive_pym] = np.nan
    if verbose:
        print(
            f'After "Ignore harvests from before this output began: discrepancy of {np.sum(~np.isnan(hdates_pym)) - expected_valid} patch-seasons'
        )

    # We need to keep some non-seasons---it's possible that "the yearY growing season" never happened (sowing conditions weren't met), but we still need something there so that we can make an array of dimension Npatch*Ngs. We do this by changing those non-seasons from NaN to -Inf before doing the filtering and reshaping, after which we'll convert them back to NaNs.

    # "In years with no sowing, pretend the first no-harvest is meaningful, unless that was intentionally ignored above."
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
    for h in np.arange(mxharvests - 1):
        if h == 0:
            continue
        elif h == 1:
            print("Warning: Untested with mxharvests > 2")
        where_nosow_py = np.where(
            nosow_py
            & ~np.any(np.isnan(hdates_pym[:, :, 0:h]), axis=2)
            & np.isnan(hdates_pym[:, :, h])
        )
        hdates_pym2[where_nosow_py[0], where_nosow_py[1], h + 1] = -np.inf
        sdates_pym2[where_nosow_py[0], where_nosow_py[1], h + 1] = -np.inf

    # "In years with sowing that are followed by inactive years, check whether the last sowing was harvested before the patch was deactivated. If not, pretend the LAST [easier to implement!] no-harvest is meaningful."
    sdates_orig_masked_pym = sdates_orig_pym.copy()
    with np.errstate(invalid="ignore"):
        sdates_le_0 = sdates_orig_masked_pym <= 0
    sdates_orig_masked_pym[np.where(sdates_le_0)] = np.nan
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")
        last_sdate_firstNgs_py = np.nanmax(sdates_orig_masked_pym[:, :-1, :], axis=2)
        last_hdate_firstNgs_py = np.nanmax(hdates_pym2[:, :-1, :], axis=2)
    with np.errstate(invalid="ignore"):
        hdate_lt_sdate = last_hdate_firstNgs_py < last_sdate_firstNgs_py
    last_sowing_not_harvested_sameyear_firstNgs_py = hdate_lt_sdate | np.isnan(
        last_hdate_firstNgs_py
    )
    inactive_lastNgs_py = inactive_py[:, 1:]
    last_sowing_never_harvested_firstNgs_py = (
        last_sowing_not_harvested_sameyear_firstNgs_py & inactive_lastNgs_py
    )
    last_sowing_never_harvested_py = np.concatenate(
        (last_sowing_never_harvested_firstNgs_py, np.full((Npatch, 1), False)), axis=1
    )
    last_sowing_never_harvested_pym = np.concatenate(
        (
            np.full((Npatch, Ngs + 1, mxharvests - 1), False),
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
    def pym_to_pg(pym, quiet=False):
        pg = np.reshape(pym, (pym.shape[0], -1))
        ok_pg = pg[~np.isnan(pg)]
        if not quiet:
            print(
                f"{ok_pg.size} included; unique N seasons = {np.unique(np.sum(~np.isnan(pg), axis=1))}"
            )
        return pg

    hdates_pg = pym_to_pg(hdates_pym3.copy(), quiet=~verbose)
    sdates_pg = pym_to_pg(sdates_pym3.copy(), quiet=True)
    if verbose:
        print(
            f'After "In years with no sowing, pretend the first no-harvest is meaningful: discrepancy of {np.sum(~np.isnan(hdates_pg)) - expected_valid} patch-seasons'
        )

    # "Ignore any harvests that were planted in the final year, because some cells will have incomplete growing seasons for the final year."
    with np.errstate(invalid="ignore"):
        hdates_ge_sdates = hdates_pg[:, -mxharvests:] >= sdates_pg[:, -mxharvests:]
    lastyear_complete_season = hdates_ge_sdates | np.isinf(hdates_pg[:, -mxharvests:])

    def ignore_lastyear_complete_season(pg, excl, mxharvests):
        tmp_L = pg[:, :-mxharvests]
        tmp_R = pg[:, -mxharvests:]
        tmp_R[np.where(excl)] = np.nan
        pg = np.concatenate((tmp_L, tmp_R), axis=1)
        return pg

    hdates_pg2 = ignore_lastyear_complete_season(
        hdates_pg.copy(), lastyear_complete_season, mxharvests
    )
    sdates_pg2 = ignore_lastyear_complete_season(
        sdates_pg.copy(), lastyear_complete_season, mxharvests
    )
    is_valid = ~np.isnan(hdates_pg2)
    is_fake = np.isneginf(hdates_pg2)
    is_fake = np.reshape(is_fake[is_valid], (this_ds.dims["patch"], Ngs))
    discrepancy = np.sum(is_valid) - expected_valid
    unique_Nseasons = np.unique(np.sum(is_valid, axis=1))
    if verbose:
        print(
            f'After "Ignore any harvests that were planted in the final year, because other cells will have incomplete growing seasons for the final year": discrepancy of {discrepancy} patch-seasons'
        )
        try:
            import pandas as pd

            bc = np.bincount(np.sum(is_valid, axis=1))
            bc = bc[bc > 0]
            df = pd.DataFrame({"Ngs": unique_Nseasons, "Count": bc})
            print(df)
        except:
            print(f"unique N seasons = {unique_Nseasons}")
        print(" ")

    # Create Dataset with time axis as "gs" (growing season) instead of what CLM puts out
    if discrepancy == 0:
        this_ds_gs = set_up_ds_with_gs_axis(this_ds)
        for v in this_ds.data_vars:
            if this_ds[v].dims != ("time", "mxharvests", "patch") or (myVars and v not in myVars):
                continue

            # Set invalid values to NaN
            da_yhp = this_ds[v].copy()
            da_yhp = da_yhp.where(~np.isneginf(da_yhp))

            # Remove the nans and reshape to patches*growingseasons
            da_pyh = da_yhp.transpose("patch", "time", "mxharvests")
            ar_pg = np.reshape(da_pyh.values, (this_ds.dims["patch"], -1))
            ar_valid_pg = np.reshape(ar_pg[is_valid], (this_ds.dims["patch"], Ngs))
            # Change -infs to nans
            ar_valid_pg[is_fake] = np.nan
            # Save as DataArray to new Dataset, stripping _PERHARV from variable name
            newname = v.replace("_PERHARV", "")
            if newname in this_ds_gs:
                raise RuntimeError(f"{newname} already in dataset!")
            da_pg = xr.DataArray(
                data=ar_valid_pg,
                coords=[this_ds_gs.coords["patch"], this_ds_gs.coords["gs"]],
                name=newname,
                attrs=da_yhp.attrs,
            )
            this_ds_gs[newname] = da_pg
            this_ds_gs[newname].attrs["units"] = this_ds[v].attrs["units"]
    else:
        # Print details about example bad patch(es)
        if min(unique_Nseasons) < Ngs:
            print(f"Too few seasons (min {min(unique_Nseasons)} < {Ngs})")
            p = np.where(np.sum(~np.isnan(hdates_pg2), axis=1) == min(unique_Nseasons))[0][0]
            print_onepatch_wrongNgs(
                p,
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
        if max(unique_Nseasons) > Ngs:
            print(f"Too many seasons (max {max(unique_Nseasons)} > {Ngs})")
            p = np.where(np.sum(~np.isnan(hdates_pg2), axis=1) == max(unique_Nseasons))[0][0]
            print_onepatch_wrongNgs(
                p,
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
            f"Can't convert time*mxharvests axes to growingseason axis: discrepancy of {discrepancy} patch-seasons"
        )

    # Preserve units
    for v1 in this_ds_gs:
        v0 = v1
        if v0 not in this_ds:
            v0 += "_PERHARV"
        if v0 not in this_ds:
            continue
        if "units" in this_ds[v0].attrs:
            this_ds_gs[v1].attrs["units"] = this_ds[v0].attrs["units"]

    if incl_orig:
        return this_ds_gs, this_ds
    else:
        return this_ds_gs


# Minimum harvest threshold allowed in PlantCrop()
# Was 50 before cropcal runs 2023-01-28
def default_gdd_min():
    return 1.0


# Get growing season lengths from a DataArray of hdate-sdate
def get_gs_len_da(this_da):
    tmp = this_da.values
    with np.errstate(invalid="ignore"):
        tmp_lt_0 = tmp < 0
    tmp[tmp_lt_0] = 365 + tmp[tmp_lt_0]
    this_da.values = tmp
    this_da.attrs["units"] = "days"
    return this_da


def get_pct_harv_at_mature(harvest_reason_da):
    Nharv_at_mature = len(np.where(harvest_reason_da.values == 1)[0])
    with np.errstate(invalid="ignore"):
        harv_reason_gt_0 = harvest_reason_da.values > 0
    Nharv = len(np.where(harv_reason_gt_0)[0])
    if Nharv == 0:
        return np.nan
    pct_harv_at_mature = Nharv_at_mature / Nharv * 100
    pct_harv_at_mature = np.format_float_positional(
        pct_harv_at_mature, precision=2, unique=False, fractional=False, trim="k"
    )  # Round to 2 significant digits
    return pct_harv_at_mature


def import_max_gs_length(paramfile_dir, my_clm_ver, my_clm_subver):
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


# E.g. import_rx_dates("sdate", sdates_rx_file, dates_ds0_orig)
def import_rx_dates(var_prefix, date_inFile, dates_ds, set_neg1_to_nan=True):
    # Get run info:
    # Max number of growing seasons per year
    if "mxsowings" in dates_ds:
        mxsowings = dates_ds.dims["mxsowings"]
    else:
        mxsowings = 1

    # Which vegetation types were simulated?
    itype_veg_toImport = np.unique(dates_ds.patches1d_itype_veg)

    date_varList = []
    for i in itype_veg_toImport:
        for g in np.arange(mxsowings):
            thisVar = f"{var_prefix}{g+1}_{i}"
            date_varList = date_varList + [thisVar]

    ds = utils.import_ds(date_inFile, myVars=date_varList)

    did_warn = False
    for v in ds:
        v_new = v.replace(var_prefix, "gs")
        ds = ds.rename({v: v_new})

        # Set -1 prescribed GDD values to NaN. Only warn the first time.
        if set_neg1_to_nan and var_prefix == "gdd" and v_new != v and np.any(ds[v_new].values < 0):
            if np.any((ds[v_new].values < 0) & (ds[v_new].values != -1)):
                raise RuntimeError(f"Unexpected negative value in {v}")
            if not did_warn:
                print(f"Setting -1 rx GDD values to NaN")
                did_warn = True
            ds[v_new] = ds[v_new].where(ds[v_new] != -1)

    return ds


def import_output(
    filename,
    myVars,
    y1=None,
    yN=None,
    myVegtypes=utils.define_mgdcrop_list(),
    sdates_rx_ds=None,
    gdds_rx_ds=None,
    verbose=False,
    incl_irrig=True,
):
    # Import
    this_ds = utils.import_ds(filename, myVars=myVars, myVegtypes=myVegtypes)

    # Trim to years of interest (do not include extra year needed for finishing last growing season)
    if y1 and yN:
        this_ds = check_and_trim_years(y1, yN, this_ds)
    else:  # Assume including all growing seasons except last complete one are "of interest"
        y1 = this_ds.time.values[0].year
        yN = this_ds.time.values[-1].year - 2
        this_ds = check_and_trim_years(y1, yN, this_ds)

    # What vegetation types are included?
    vegtype_list = [
        x for x in this_ds.vegtype_str.values if x in this_ds.patches1d_itype_veg_str.values
    ]

    # Check for consistency among sowing/harvest date/year info
    date_vars = ["SDATES_PERHARV", "SYEARS_PERHARV", "HDATES", "HYEARS"]
    all_nan = np.full(this_ds[date_vars[0]].shape, True)
    all_nonpos = np.full(this_ds[date_vars[0]].shape, True)
    all_pos = np.full(this_ds[date_vars[0]].shape, True)
    for v in date_vars:
        all_nan = all_nan & np.isnan(this_ds[v].values)
        with np.errstate(invalid="ignore"):
            all_nonpos = all_nonpos & (this_ds[v].values <= 0)
            all_pos = all_pos & (this_ds[v].values > 0)
    if np.any(np.bitwise_not(all_nan | all_nonpos | all_pos)):
        raise RuntimeError("Inconsistent missing/present values on mxharvests axis")

    # When doing transient runs, it's somehow possible for crops in newly-active patches to be *already alive*. They even have a sowing date (idop)! This will of course not show up in SDATES, but it does show up in SDATES_PERHARV.
    # I could put the SDATES_PERHARV dates into where they "should" be, but instead I'm just going to invalidate those "seasons."
    #
    # In all but the last calendar year, which patches had no sowing?
    no_sowing_yp = np.all(np.isnan(this_ds.SDATES.values[:-1, :, :]), axis=1)
    # In all but the first calendar year, which harvests' jdays are < their sowings' jdays? (Indicates sowing the previous calendar year.)
    with np.errstate(invalid="ignore"):
        hsdate1_gt_hdate1_yp = (
            this_ds.SDATES_PERHARV.values[1:, 0, :] > this_ds.HDATES.values[1:, 0, :]
        )
    # Where both, we have the problem.
    falsely_alive_yp = no_sowing_yp & hsdate1_gt_hdate1_yp
    if np.any(falsely_alive_yp):
        print(
            f"Warning: {np.sum(falsely_alive_yp)} patch-seasons being ignored: Seemingly sown the year before harvest, but no sowings occurred that year."
        )
        falsely_alive_yp = np.concatenate(
            (np.full((1, this_ds.dims["patch"]), False), falsely_alive_yp), axis=0
        )
        falsely_alive_y1p = np.expand_dims(falsely_alive_yp, axis=1)
        dummy_false_y1p = np.expand_dims(np.full_like(falsely_alive_yp, False), axis=1)
        falsely_alive_yhp = np.concatenate((falsely_alive_y1p, dummy_false_y1p), axis=1)
        for v in this_ds.data_vars:
            if this_ds[v].dims != ("time", "mxharvests", "patch"):
                continue
            this_ds[v] = this_ds[v].where(~falsely_alive_yhp)

    def check_no_negative(this_ds_in, varList_no_negative, which_file, verbose=False):
        tiny_negOK = 1e-12
        this_ds = this_ds_in.copy()
        for v in this_ds:
            if not any(x in v for x in varList_no_negative):
                continue
            the_min = np.nanmin(this_ds[v].values)
            if the_min < 0:
                if np.abs(the_min) <= tiny_negOK:
                    if verbose:
                        print(
                            f"Tiny negative value(s) in {v} (abs <= {tiny_negOK}) being set to 0 ({which_file})"
                        )
                else:
                    print(
                        f"WARNING: Unexpected negative value(s) in {v}; minimum {the_min} ({which_file})"
                    )
                values = this_ds[v].copy().values
                with np.errstate(invalid="ignore"):
                    do_setto_0 = (values < 0) & (values >= -tiny_negOK)
                values[np.where(do_setto_0)] = 0
                this_ds[v] = xr.DataArray(
                    values, coords=this_ds[v].coords, dims=this_ds[v].dims, attrs=this_ds[v].attrs
                )

            elif verbose:
                print(f"No negative value(s) in {v}; min {the_min} ({which_file})")
        return this_ds

    def check_no_zeros(this_ds, varList_no_zero, which_file):
        for v in this_ds:
            if not any(x in v for x in varList_no_zero):
                continue
            if np.any(this_ds[v].values == 0):
                print(f"WARNING: Unexpected zero(s) in {v} ({which_file})")
            elif verbose:
                print(f"No zero value(s) in {v} ({which_file})")

    # Check for no zero values where there shouldn't be
    varList_no_zero = ["DATE", "YEAR"]
    check_no_zeros(this_ds, varList_no_zero, "original file")

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
    for v in ["HUIFRAC", "HUIFRAC_PERHARV"]:
        this_ds_gs[v].attrs["units"] = "Fraction of required"

    # Avoid tiny negative values
    varList_no_negative = ["GRAIN", "REASON", "GDD", "HUI", "YEAR", "DATE", "GSLEN"]
    this_ds_gs = check_no_negative(this_ds_gs, varList_no_negative, "new file", verbose=verbose)

    # Check for no zero values where there shouldn't be
    varList_no_zero = ["REASON", "DATE"]
    check_no_zeros(this_ds_gs, varList_no_zero, "new file")

    # Check that e.g., GDDACCUM <= HUI
    for vars in [["GDDACCUM", "HUI"], ["SYEARS", "HYEARS"]]:
        if all(v in this_ds_gs for v in vars):
            check_v0_le_v1(this_ds_gs, vars, both_nan_ok=True, throw_error=True)

    # Check that prescribed calendars were obeyed
    if sdates_rx_ds:
        check_rx_obeyed(vegtype_list, sdates_rx_ds, this_ds, "this_ds", "SDATES")
    if gdds_rx_ds:
        check_rx_obeyed(
            vegtype_list,
            gdds_rx_ds,
            this_ds,
            "this_ds",
            "SDATES",
            "GDDHARV",
            gdd_min=default_gdd_min(),
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

    # Import irrigation data, if doing so
    if incl_irrig:
        # Monthly use and demand
        pattern = os.path.join(os.path.dirname(filename), "*.h2.*.nc")
        irrig_file_patches = glob.glob(pattern)
        if irrig_file_patches:
            if len(irrig_file_patches) > 1:
                raise RuntimeError(f"Expected at most 1 *.h2.*.nc; found {len(irrig_file_patches)}")
            irrig_file_patches = irrig_file_patches[0]
            irrig_ds_patches = utils.import_ds(
                irrig_file_patches,
                myVegtypes=myVegtypes,
                chunks={"time": 12},
                myVars=["QIRRIG_DEMAND", "QIRRIG_DRIP", "QIRRIG_SPRINKLER"],
            )
            irrig_ds_patches = time_units_and_trim_mth(irrig_ds_patches, y1, yN)

            # Combine drip + sprinkler irrigation
            irrig_ds_patches["QIRRIG_APPLIED"] = (
                irrig_ds_patches["QIRRIG_DRIP"] + irrig_ds_patches["QIRRIG_SPRINKLER"]
            )
            irrig_ds_patches["QIRRIG_APPLIED"].attrs = irrig_ds_patches["QIRRIG_DEMAND"].attrs
            irrig_ds_patches["QIRRIG_APPLIED"].attrs[
                "long_name"
            ] = "water added via drip or sprinkler irrigation"

            # Calculate unfulfilled demand
            irrig_ds_patches["QIRRIG_UNFULFILLED"] = (
                irrig_ds_patches["QIRRIG_DEMAND"] - irrig_ds_patches["QIRRIG_APPLIED"]
            )
            irrig_ds_patches["QIRRIG_UNFULFILLED"].attrs = irrig_ds_patches["QIRRIG_DEMAND"].attrs
            irrig_ds_patches["QIRRIG_UNFULFILLED"].attrs[
                "long_name"
            ] = "irrigation demand unable to be filled"

            vars_to_save = ["QIRRIG_DEMAND", "QIRRIG_APPLIED", "QIRRIG_UNFULFILLED"]

            # Append _PATCH to distinguish from eventual gridcell sums (_GRID)
            rename_dict = {}
            for v in vars_to_save:
                rename_dict[v] = v + "_PATCH"
            irrig_ds_patches = irrig_ds_patches.rename(rename_dict)
            vars_to_save = [v + "_PATCH" for v in vars_to_save]

            # Finish processing
            this_ds_gs = process_monthly_irrig(this_ds_gs, irrig_ds_patches, vars_to_save)

        # Monthly withdrawals and availability
        pattern = os.path.join(os.path.dirname(filename), "*.h3.*.nc")
        irrig_file_grid = glob.glob(pattern)
        if irrig_file_grid:
            if len(irrig_file_grid) > 1:
                raise RuntimeError(f"Expected at most 1 *.h3.*.nc; found {len(irrig_file_grid)}")
            irrig_file_grid = irrig_file_grid[0]
            irrig_ds_grid = utils.import_ds(
                irrig_file_grid,
                chunks={"time": 12},
                myVars=[
                    "area",
                    "grid1d_ixy",
                    "grid1d_jxy",
                    "grid1d_lon",
                    "grid1d_lat",
                    "QIRRIG_FROM_GW_CONFINED",
                    "QIRRIG_FROM_GW_UNCONFINED",
                    "QIRRIG_FROM_SURFACE",
                    "VOLRMCH",
                ],
            )
            irrig_ds_grid = time_units_and_trim_mth(irrig_ds_grid, y1, yN)

            # Ensure that no groundwater was used
            if np.any(
                irrig_ds_grid["QIRRIG_FROM_GW_CONFINED"]
                + irrig_ds_grid["QIRRIG_FROM_GW_UNCONFINED"]
                > 0
            ):
                raise RuntimeError("Unexpectedly found some irrigation using groundwater")

            # Rename this variable so that it has a QIRRIG prefix
            irrig_ds_grid = irrig_ds_grid.rename({"VOLRMCH": "IRRIG_SUPPLY"})
            irrig_ds_grid["IRRIG_SUPPLY"].attrs["long_name"] = (
                irrig_ds_grid["IRRIG_SUPPLY"].attrs["long_name"] + " (aka VOLRMCH)"
            )

            vars_to_save = ["QIRRIG_FROM_SURFACE", "IRRIG_SUPPLY"]

            # Append _GRID to distinguish from patch-level irrigation data
            rename_dict = {}
            for v in vars_to_save:
                rename_dict[v] = v + "_GRID"
            irrig_ds_grid = irrig_ds_grid.rename(rename_dict)
            vars_to_save = [v + "_GRID" for v in vars_to_save]

            # Finish processing
            this_ds_gs = process_monthly_irrig(this_ds_gs, irrig_ds_grid, vars_to_save)

            # Calculate irrigation as fraction of main river channel volume
            # (Do it here instead of process_monthly_irrig() because monthly doesn't add to annual.)
            for t in ["MTH", "ANN"]:
                this_ds_gs[f"QIRRIG_FROM_SURFACE_FRAC_RIVER_GRID_{t}"] = (
                    this_ds_gs[f"QIRRIG_FROM_SURFACE_GRID_{t}"]
                    / this_ds_gs[f"IRRIG_SUPPLY_GRID_{t}"]
                )
                this_ds_gs[f"QIRRIG_FROM_SURFACE_FRAC_RIVER_GRID_{t}"] = this_ds_gs[
                    f"QIRRIG_FROM_SURFACE_FRAC_RIVER_GRID_{t}"
                ].assign_coords(this_ds_gs[f"QIRRIG_FROM_SURFACE_GRID_{t}"].coords)

            # Ensure gridcells are same-ordered between patch- and gridcell-level datasets
            this_ds_gs = this_ds_gs.assign_coords({"gridcell": this_ds_gs["gridcell"] + 1})
            irrig_ds_grid = irrig_ds_grid.assign_coords({"gridcell": irrig_ds_grid["gridcell"] + 1})
            test_gi = int(this_ds_gs["gridcell"].shape[0] / 3)  # An arbitrary gridcell to test
            test_g = irrig_ds_grid["gridcell"].values[test_gi]
            test_pi = np.where(this_ds_gs["patches1d_gi"] == test_g)[0][0]
            test_p_lon = this_ds_gs["patches1d_lon"].isel(patch=test_pi)
            test_p_lat = this_ds_gs["patches1d_lat"].isel(patch=test_pi)
            test_g_lon = irrig_ds_grid["grid1d_lon"].isel(gridcell=test_gi)
            test_g_lat = irrig_ds_grid["grid1d_lat"].isel(gridcell=test_gi)
            if test_p_lon != test_g_lon or test_p_lat != test_g_lat:
                print(f"{test_p_lon} =? {test_g_lon}")
                print(f"{test_p_lat} =? {test_g_lat}")
                raise RuntimeError("Gridcell indexing mismatch")

            # Save gridcell info
            for v in irrig_ds_grid:
                if "grid" in v:
                    this_ds_gs[v] = irrig_ds_grid[v]

            # Save gridcell area
            this_ds_gs["AREA_GRID"] = ungrid(
                irrig_ds_grid["area"],
                irrig_ds_grid,
                target_dim="gridcell",
                lon="grid1d_ixy",
                lat="grid1d_jxy",
            )
            if this_ds_gs["AREA_GRID"].attrs["units"] != "m^2":
                if this_ds_gs["AREA_GRID"].attrs["units"] == "km^2":
                    to_m2 = 1e6
                else:
                    raise RuntimeError(
                        f"Unsure how to convert {this_ds_gs['AREA_GRID'].attrs['units']} to m^2"
                    )
                this_ds_gs["AREA_GRID"] *= to_m2
                this_ds_gs["AREA_GRID"].attrs["units"] = "m^2"

    return this_ds_gs


# Set up empty Dataset with time axis as "gs" (growing season) instead of what CLM puts out.
# Includes all the same variables as the input dataset, minus any that had dimensions mxsowings or mxharvests.
def set_up_ds_with_gs_axis(ds_in):
    # Get the data variables to include in the new dataset
    data_vars = dict()
    for v in ds_in.data_vars:
        if not any([x in ["mxsowings", "mxharvests"] for x in ds_in[v].dims]):
            data_vars[v] = ds_in[v]
    # Set up the new dataset
    gs_years = [t.year - 1 for t in ds_in.time.values[:-1]]
    coords = ds_in.coords
    coords["gs"] = gs_years
    ds_out = xr.Dataset(data_vars=data_vars, coords=coords, attrs=ds_in.attrs)
    return ds_out
