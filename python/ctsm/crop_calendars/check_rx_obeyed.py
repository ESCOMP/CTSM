"""
Check that prescribed crop calendars were obeyed
"""

import numpy as np

import ctsm.crop_calendars.cropcal_utils as utils
from ctsm.crop_calendars.cropcal_constants import DEFAULT_GDD_MIN

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments


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


def check_rx_obeyed_handle_gdharv(output_var, gdd_min, ds_thisveg, rx_array):
    """
    In check_rx_obeyed(), account for the GDD harvest threshold minimum set in PlantCrop()
    """
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
    return gdd_min, unique_harvest_reasons, pct_harv_at_mature


def check_rx_obeyed_setup(dates_ds, which_ds, output_var, verbose):
    """
    Various setup steps for check_rx_obeyed()
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

    return all_ok, diff_str_list, gdd_tolerance


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


def summarize_results(which_ds, output_var, verbose, all_ok, gdd_tolerance, diffs_eg_txt):
    """
    Summarize results
    """
    bad = True
    if all_ok == 2:
        bad = False
        print(f"✅ {which_ds}: Prescribed {output_var} always obeyed")
    elif all_ok == 1:
        bad = False
        print(
            f"🟨 {which_ds}: Prescribed {output_var} *not* always obeyed, but acceptable (diffs <= "
            + f"{gdd_tolerance})"
        )
    elif not verbose:
        print(f"❌ {which_ds}: Prescribed {output_var} *not* always obeyed. E.g., {diffs_eg_txt}")
    return bad


def check_rx_obeyed(
    vegtype_list, rx_ds, dates_ds, which_ds, output_var, *, gdd_min=None, verbose=False
):
    """
    Check that prescribed crop calendars were obeyed
    """
    all_ok, diff_str_list, gdd_tolerance = check_rx_obeyed_setup(
        dates_ds, which_ds, output_var, verbose
    )

    diffs_eg_txt = None
    for vegtype_str in vegtype_list:
        thisveg_patches = np.where(dates_ds.patches1d_itype_veg_str == vegtype_str)[0]
        if thisveg_patches.size == 0:
            continue
        ds_thisveg = dates_ds.isel(patch=thisveg_patches)

        vegtype_int = utils.vegtype_str2int(vegtype_str)[0]
        rx_da = rx_ds[f"gs1_{vegtype_int}"]
        rx_array = rx_da.values[
            ds_thisveg.patches1d_jxy.values.astype(int) - 1,
            ds_thisveg.patches1d_ixy.values.astype(int) - 1,
        ]
        rx_array = np.expand_dims(rx_array, axis=1)
        sim_array = ds_thisveg[output_var].values
        sim_array_dims = ds_thisveg[output_var].dims

        # Ignore patches without prescribed value
        with np.errstate(invalid="ignore"):
            rx_array[np.where(rx_array < 0)] = np.nan

        # Account for...
        if "GDDHARV" in output_var:
            # ...GDD harvest threshold minimum set in PlantCrop()
            gdd_min, unique_harvest_reasons, pct_harv_at_mature = check_rx_obeyed_handle_gdharv(
                output_var, gdd_min, ds_thisveg, rx_array
            )

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
                    ds_thisveg.patches1d_lon,
                    ds_thisveg.patches1d_lat,
                )
                max_diff, max_lon, max_lat, max_gs, max_rx = get_extreme_info(
                    diff_array,
                    rx_array,
                    np.nanmax,
                    sim_array_dims,
                    dates_ds.gs,
                    ds_thisveg.patches1d_lon,
                    ds_thisveg.patches1d_lat,
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

    bad = summarize_results(which_ds, output_var, verbose, all_ok, gdd_tolerance, diffs_eg_txt)

    return bad
