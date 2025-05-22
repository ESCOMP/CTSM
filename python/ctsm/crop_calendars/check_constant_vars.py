"""
For variables that should stay constant, make sure they are
"""

import numpy as np
from ctsm.crop_calendars.cropcal_module import import_rx_dates

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments


def check_one_constant_var_setup(this_ds, case, var):
    """
    Various setup steps for check_one_constant_var()
    """
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

    return time_coord, this_da, ra_sp, incl_patches, str_list, rx_ds, bad_patches


def loop_through_bad_patches(
    verbose,
    emojus,
    var,
    everything_ok,
    str_list,
    rx_ds,
    time_1,
    t1_yr,
    t1_vals,
    timestep,
    t_yr,
    t_vals,
    bad_patches_this_time,
    found_in_rx,
    vary_patches,
    vary_lons,
    vary_lats,
    vary_crops,
    vary_crops_int,
    any_bad,
):
    """
    Loop through and check any patches that were "bad" according to check_constant_vars().

    This is pretty inefficient, but it works.
    """
    patch = None  # In case bad_patches_this_time is empty
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
                raise RuntimeError(f"lon {this_lon} lat {this_lat} {this_crop} not in rx dataset?")

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
                    f"{this_str}: {t1_yr} {var} {t1_val_print}, {t_yr} {var} " + f"{t_val_print}"
                )
        else:
            if everything_ok:
                print(f"{emojus} CLM output {var} unexpectedly vary over time:")
                everything_ok = False
            print(f"{var} timestep {timestep} does not match timestep {time_1}")
            break
    return any_bad, patch


def ensure_all_patches_checked(this_ds, this_da, ra_sp, incl_patches):
    """
    In check_one_constant_var(), make sure every patch was checked once (or is all-NaN except
    possibly final season)
    """
    incl_patches = np.sort(incl_patches)
    if not np.array_equal(incl_patches, np.unique(incl_patches)):
        raise RuntimeError("Patch(es) checked more than once!")
    incl_patches = list(incl_patches)
    incl_patches += list(
        np.where(
            np.all(
                np.isnan(ra_sp[:-1,]),
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
                raise RuntimeError(
                    f"Not all patches checked! E.g., {patch}: {this_da.isel(patch=patch).values}"
                )


def check_one_constant_var_loop_through_timesteps(
    this_ds,
    ignore_nan,
    verbose,
    emojus,
    var,
    everything_ok,
    time_coord,
    this_da,
    str_list,
    rx_ds,
    time_1,
    these_patches,
    t1_yr,
    t1_vals,
    any_bad,
    any_bad_before_checking_rx,
    bad_patches,
):
    """
    In check_one_constant_var(), loop through timesteps
    """
    found_in_rx = None
    for timestep in np.arange(time_1 + 1, this_ds.dims[time_coord]):
        t_yr = this_ds[time_coord].values[timestep]
        t_vals = np.squeeze(this_da.isel({time_coord: timestep, "patch": these_patches}).values)
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
                    rx_ds[rx_var].sel(lon=vary_lons_this_crop, lat=vary_lats_this_crop).values
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

            # Loop through and check any patches that were "bad"
            any_bad = loop_through_bad_patches(
                verbose,
                emojus,
                var,
                everything_ok,
                str_list,
                rx_ds,
                time_1,
                t1_yr,
                t1_vals,
                timestep,
                t_yr,
                t_vals,
                bad_patches_this_time,
                found_in_rx,
                vary_patches,
                vary_lons,
                vary_lats,
                vary_crops,
                vary_crops_int,
                any_bad,
            )

    return any_bad_before_checking_rx, bad_patches, found_in_rx, any_bad


def check_one_constant_var(
    this_ds, case, ignore_nan, verbose, emojus, var, any_bad, any_bad_before_checking_rx
):
    """
    Ensure that a variable that should be constant actually is
    """
    everything_ok = True

    (
        time_coord,
        this_da,
        ra_sp,
        incl_patches,
        str_list,
        rx_ds,
        bad_patches,
    ) = check_one_constant_var_setup(this_ds, case, var)

    for time_1 in np.arange(this_ds.dims[time_coord] - 1):
        condn = ~np.isnan(ra_sp[time_1, ...])
        if time_1 > 0:
            condn = np.bitwise_and(condn, np.all(np.isnan(ra_sp[:time_1, ...]), axis=0))
        these_patches = np.where(condn)[0]
        if these_patches.size == 0:
            continue
        these_patches = list(np.where(condn)[0])
        incl_patches += these_patches

        t1_yr = this_ds[time_coord].values[time_1]
        t1_vals = np.squeeze(this_da.isel({time_coord: time_1, "patch": these_patches}).values)

        (
            any_bad_before_checking_rx,
            bad_patches,
            found_in_rx,
            any_bad,
        ) = check_one_constant_var_loop_through_timesteps(
            this_ds,
            ignore_nan,
            verbose,
            emojus,
            var,
            everything_ok,
            time_coord,
            this_da,
            str_list,
            rx_ds,
            time_1,
            these_patches,
            t1_yr,
            t1_vals,
            any_bad,
            any_bad_before_checking_rx,
            bad_patches,
        )

    if verbose and any_bad:
        print(f"{emojus} CLM output {var} unexpectedly vary over time:")
        str_list.sort()
        if found_in_rx is None:
            raise RuntimeError("Somehow any_bad True but found_in_rx None")
        if rx_ds and np.any(~found_in_rx):  # pylint: disable=invalid-unary-operand-type
            str_list = [
                "*: Not found in prescribed input file (maybe minor lon/lat mismatch)"
            ] + str_list
        elif not rx_ds:
            str_list = ["(No rx file checked)"] + str_list
        print("\n".join(str_list))

    # Make sure every patch was checked once (or is all-NaN except possibly final season)
    ensure_all_patches_checked(this_ds, this_da, ra_sp, incl_patches)

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

    return any_bad, any_bad_before_checking_rx, bad_patches


def check_constant_vars(
    this_ds, case, ignore_nan, *, const_growing_seasons=None, verbose=True, throw_error=True
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
        any_bad, any_bad_before_checking_rx, bad_patches = check_one_constant_var(
            this_ds, case, ignore_nan, verbose, emojus, var, any_bad, any_bad_before_checking_rx
        )

    if any_bad and throw_error:
        raise RuntimeError("Stopping due to failed check_constant_vars().")

    bad_patches = np.unique(bad_patches)
    return [int(p) for p in bad_patches], any_bad
