"""
Convert time*mxharvests axes to growingseason axis
"""

import warnings
import sys
import numpy as np
import xarray as xr
from ctsm.crop_calendars.cropcal_utils import get_integer_years

try:
    import pandas as pd
except ModuleNotFoundError:
    pass

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments


def pym_to_pg(pym_array, quiet=False):
    """
    In convert_axis_time2gs(), convert year x month array to growingseason axis
    """
    pg_array = np.reshape(pym_array, (pym_array.shape[0], -1))
    ok_pg = pg_array[~np.isnan(pg_array)]
    if not quiet:
        print(
            f"{ok_pg.size} included; unique N seasons = "
            + f"{np.unique(np.sum(~np.isnan(pg_array), axis=1))}"
        )
    return pg_array


def ignore_lastyear_complete_season(pg_array, excl, mxharvests):
    """
    Helper function for convert_axis_time2gs()
    """
    tmp_l = pg_array[:, :-mxharvests]
    tmp_r = pg_array[:, -mxharvests:]
    tmp_r[np.where(excl)] = np.nan
    pg_array = np.concatenate((tmp_l, tmp_r), axis=1)
    return pg_array


def convert_axis_time2gs_setup(this_ds, verbose):
    """
    Various setup steps for convert_axis_time2gs_setup()
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
    return n_patch, n_gs, expected_valid, mxharvests, hdates_ymp, hdates_pym, sdates_ymp, sdates_pym


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
    gs_years = get_integer_years(ds_in)[:-1]
    coords = ds_in.coords
    coords["gs"] = gs_years
    ds_out = xr.Dataset(data_vars=data_vars, coords=coords, attrs=ds_in.attrs)
    return ds_out


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


def handle_years_with_no_sowing(this_ds, mxharvests, hdates_pym, sdates_pym):
    """
    In years with no sowing, pretend the first no-harvest is meaningful, unless that was
    intentionally ignored earlier in convert_axis_time2gs().
    """
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
        if harvest_index == 1:
            print("Warning: Untested with mxharvests > 2")
        where_nosow_py = np.where(
            nosow_py
            & ~np.any(np.isnan(hdates_pym[:, :, 0:harvest_index]), axis=2)
            & np.isnan(hdates_pym[:, :, harvest_index])
        )
        hdates_pym2[where_nosow_py[0], where_nosow_py[1], harvest_index + 1] = -np.inf
        sdates_pym2[where_nosow_py[0], where_nosow_py[1], harvest_index + 1] = -np.inf
    return sdates_orig_pym, hdates_pym2, sdates_pym2


def handle_years_with_sowing_then_inactive(
    verbose,
    n_patch,
    n_gs,
    expected_valid,
    mxharvests,
    inactive_py,
    sdates_orig_pym,
    hdates_pym2,
    sdates_pym2,
):
    """
    In years with sowing that are followed by inactive years, check whether the last sowing was
    harvested before the patch was deactivated. If not, pretend the LAST [easier to implement!]
    no-harvest is meaningful.
    """
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

    hdates_pg = pym_to_pg(hdates_pym3.copy(), quiet=~verbose)
    sdates_pg = pym_to_pg(sdates_pym3.copy(), quiet=True)
    if verbose:
        print(
            "After 'In years with no sowing, pretend the first no-harvest is meaningful: "
            + f"discrepancy of {np.sum(~np.isnan(hdates_pg)) - expected_valid} patch-seasons"
        )

    return hdates_pym3, sdates_pym3, hdates_pg, sdates_pg


def ignore_harvests_planted_in_final_year(
    this_ds, verbose, n_gs, expected_valid, mxharvests, hdates_pg, sdates_pg
):
    """
    Ignore any harvests that were planted in the final year, because some cells will have
    incomplete growing seasons for the final year.
    """
    with np.errstate(invalid="ignore"):
        hdates_ge_sdates = hdates_pg[:, -mxharvests:] >= sdates_pg[:, -mxharvests:]
    lastyear_complete_season = hdates_ge_sdates | np.isinf(hdates_pg[:, -mxharvests:])

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
    return hdates_pg2, sdates_pg2, is_valid, is_fake, discrepancy, unique_n_seasons


def create_dataset(
    this_ds,
    my_vars,
    n_gs,
    hdates_ymp,
    hdates_pym,
    sdates_ymp,
    sdates_pym,
    hdates_pym2,
    sdates_pym2,
    hdates_pym3,
    sdates_pym3,
    hdates_pg,
    sdates_pg,
    hdates_pg2,
    sdates_pg2,
    is_valid,
    is_fake,
    discrepancy,
    unique_n_seasons,
):
    """
    Create Dataset with time axis as "gs" (growing season) instead of what CLM puts out
    """
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
    return this_ds_gs


def convert_axis_time2gs(this_ds, verbose=False, my_vars=None, incl_orig=False):
    """
    Convert time*mxharvests axes to growingseason axis
    """

    (
        n_patch,
        n_gs,
        expected_valid,
        mxharvests,
        hdates_ymp,
        hdates_pym,
        sdates_ymp,
        sdates_pym,
    ) = convert_axis_time2gs_setup(this_ds, verbose)

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
    sdates_orig_pym, hdates_pym2, sdates_pym2 = handle_years_with_no_sowing(
        this_ds, mxharvests, hdates_pym, sdates_pym
    )

    # "In years with sowing that are followed by inactive years, check whether the last sowing was
    # harvested before the patch was deactivated. If not, pretend the LAST [easier to implement!]
    # no-harvest is meaningful."
    hdates_pym3, sdates_pym3, hdates_pg, sdates_pg = handle_years_with_sowing_then_inactive(
        verbose,
        n_patch,
        n_gs,
        expected_valid,
        mxharvests,
        inactive_py,
        sdates_orig_pym,
        hdates_pym2,
        sdates_pym2,
    )

    # "Ignore any harvests that were planted in the final year, because some cells will have
    # incomplete growing seasons for the final year."
    (
        hdates_pg2,
        sdates_pg2,
        is_valid,
        is_fake,
        discrepancy,
        unique_n_seasons,
    ) = ignore_harvests_planted_in_final_year(
        this_ds, verbose, n_gs, expected_valid, mxharvests, hdates_pg, sdates_pg
    )

    # Create Dataset with time axis as "gs" (growing season) instead of what CLM puts out
    this_ds_gs = create_dataset(
        this_ds,
        my_vars,
        n_gs,
        hdates_ymp,
        hdates_pym,
        sdates_ymp,
        sdates_pym,
        hdates_pym2,
        sdates_pym2,
        hdates_pym3,
        sdates_pym3,
        hdates_pg,
        sdates_pg,
        hdates_pg2,
        sdates_pg2,
        is_valid,
        is_fake,
        discrepancy,
        unique_n_seasons,
    )

    if incl_orig:
        return this_ds_gs, this_ds
    return this_ds_gs
