"""
Functions to support generate_gdds.py
"""

# pylint: disable=too-many-lines,too-many-statements,abstract-class-instantiated
import warnings
import os
import glob
import datetime as dt
from importlib import util as importlib_util
import numpy as np
import xarray as xr

from ctsm.utils import is_instantaneous
import ctsm.crop_calendars.cropcal_utils as utils
import ctsm.crop_calendars.cropcal_module as cc
from ctsm.crop_calendars.xr_flexsel import xr_flexsel
from ctsm.crop_calendars.grid_one_variable import grid_one_variable
from ctsm.crop_calendars.import_ds import import_ds

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments

CAN_PLOT = True
try:
    # pylint: disable=wildcard-import,unused-wildcard-import
    # pylint: disable=import-error
    from ctsm.crop_calendars.cropcal_figs_module import *
    from matplotlib.transforms import Bbox
    import matplotlib.pyplot as plt

    warnings.filterwarnings(
        "ignore",
        message=(
            "__len__ for multi-part geometries is deprecated and will be removed in Shapely "
            + "2.0. Check the length of the `geoms` property instead to get the  number of "
            + "parts of a multi-part geometry."
        ),
    )
    warnings.filterwarnings(
        "ignore",
        message=(
            "Iteration over multi-part geometries is deprecated and will be removed in Shapely "
            + "2.0. Use the `geoms` property to access the constituent parts of a multi-part "
            + "geometry."
        ),
    )

    print("Will (attempt to) produce harvest requirement map figure files.")

except ModuleNotFoundError:
    print("Will NOT produce harvest requirement map figure files.")
    CAN_PLOT = False


def log(logger, string):
    """
    Simultaneously print INFO messages to console and to log file
    """
    print(string)
    logger.info(string)


def error(logger, string):
    """
    Simultaneously print ERROR messages to console and to log file
    """
    print(string)
    logger.error(string)
    raise RuntimeError(string)


def check_sdates(dates_ds, sdates_rx, outdir_figs, logger, verbose=False):
    """
    Checking that input and output sdates match
    """
    log(logger, "   Checking that input and output sdates match...")

    sdates_grid = grid_one_variable(dates_ds, "SDATES")

    all_ok = True
    any_found = False
    vegtypes_skipped = []
    vegtypes_included = []
    for i, vegtype_str in enumerate(dates_ds.vegtype_str.values):
        # Input
        vegtype_int = dates_ds.ivt.values[i]
        this_var = f"gs1_{vegtype_int}"
        if this_var not in sdates_rx:
            vegtypes_skipped = vegtypes_skipped + [vegtype_str]
            # log(logger, f"    {vt_str} ({vt}) SKIPPED...")
            continue
        vegtypes_included = vegtypes_included + [vegtype_str]
        any_found = True
        if verbose:
            log(logger, f"    {vegtype_str} ({vegtype_int})...")
        in_map = sdates_rx[this_var].squeeze(drop=True)

        # Output
        out_map = sdates_grid.sel(ivt_str=vegtype_str).squeeze(drop=True)

        # Check for differences
        diff_map = out_map - in_map
        diff_map_notnan = diff_map.values[np.invert(np.isnan(diff_map.values))]
        if np.any(diff_map_notnan):
            log(logger, f"Difference(s) found in {vegtype_str}")
            here = np.where(diff_map_notnan)
            log(logger, "in:")
            in_map_notnan = in_map.values[np.invert(np.isnan(diff_map.values))]
            log(logger, in_map_notnan[here][0:4])
            out_map_notnan = out_map.values[np.invert(np.isnan(diff_map.values))]
            log(logger, "out:")
            log(logger, out_map_notnan[here][0:4])
            log(logger, "diff:")
            log(logger, diff_map_notnan[here][0:4])
            first_diff = all_ok
            all_ok = False

            if CAN_PLOT:
                sdate_diffs_dir = os.path.join(outdir_figs, "sdate_diffs")
                if first_diff:
                    log(logger, f"Saving sdate difference figures to {sdate_diffs_dir}")
                if not os.path.exists(sdate_diffs_dir):
                    os.makedirs(sdate_diffs_dir)
                in_map.where(~np.isnan(out_map)).plot()
                plt.title(f"{vegtype_str} sdates in (masked)")
                plt.savefig(os.path.join(sdate_diffs_dir, f"{vegtype_str}.in.png"))
                plt.close()
                out_map.plot()
                plt.title(f"{vegtype_str} sdates out")
                plt.savefig(os.path.join(sdate_diffs_dir, f"{vegtype_str}.out.png"))
                plt.close()
                diff_map.plot()
                plt.title(f"{vegtype_str} sdates diff (out - in)")
                plt.savefig(os.path.join(sdate_diffs_dir, f"{vegtype_str}.diff.png"))
                plt.close()

    if not any_found:
        error(logger, "No matching variables found in sdates_rx!")

    # Sanity checks for included vegetation types
    vegtypes_skipped = np.unique([x.replace("irrigated_", "") for x in vegtypes_skipped])
    vegtypes_skipped_weird = [x for x in vegtypes_skipped if x in vegtypes_included]
    if np.array_equal(vegtypes_included, [x.replace("irrigated_", "") for x in vegtypes_included]):
        log(logger, "\nWARNING: No irrigated crops included!!!\n")
    elif vegtypes_skipped_weird:
        log(
            logger,
            "\nWarning: Some crop types had output rainfed patches but no irrigated patches: "
            + f"{vegtypes_skipped_weird}",
        )

    if all_ok:
        log(logger, "   ✅ Input and output sdates match!")
    else:
        error(logger, "   ❌ Input and output sdates differ.")


def import_rx_dates(s_or_h, date_infile, incl_patches1d_itype_veg, mxsowings, logger):
    """
    Import prescribed sowing or harvest dates
    """
    if isinstance(date_infile, xr.Dataset):
        return date_infile
    if not isinstance(date_infile, str):
        error(
            logger,
            f"Importing {s_or_h}dates_rx: Expected date_infile to be str or DataArray,"
            + f"not {type(date_infile)}",
        )

    # Which vegetation types were simulated?
    itype_veg_to_import = np.unique(incl_patches1d_itype_veg)

    date_var_list = []
    for i in itype_veg_to_import:
        for n_sowing in np.arange(mxsowings):
            this_var = f"{s_or_h}date{n_sowing+1}_{i}"
            date_var_list = date_var_list + [this_var]

    this_ds = import_ds(date_infile, my_vars=date_var_list)

    for var in this_ds:
        this_ds = this_ds.rename({var: var.replace(f"{s_or_h}date", "gs")})

    return this_ds


def this_crop_map_to_patches(lon_points, lat_points, map_ds, vegtype_int):
    """
    Given a map, get a vector of patches
    """
    # xarray pointwise indexing;
    # see https://xarray.pydata.org/en/stable/user-guide/indexing.html#more-advanced-indexing
    return (
        map_ds[f"gs1_{vegtype_int}"]
        .sel(lon=xr.DataArray(lon_points, dims="patch"), lat=xr.DataArray(lat_points, dims="patch"))
        .squeeze(drop=True)
    )


def yp_list_to_ds(yp_list, daily_ds, incl_vegtypes_str, dates_rx, longname_prefix, logger):
    """
    Get and grid mean GDDs in GGCMI growing season
    """
    # Get means
    warnings.filterwarnings(
        "ignore", message="Mean of empty slice"
    )  # Happens when you do np.nanmean() of an all-NaN array (or slice, if doing selected axis/es)
    p_list = [np.nanmean(x, axis=0) if not isinstance(x, type(None)) else x for x in yp_list]
    warnings.filterwarnings("always", message="Mean of empty slice")

    if isinstance(incl_vegtypes_str, xr.DataArray):
        incl_vegtypes_str = incl_vegtypes_str.values

    # Grid
    ds_out = xr.Dataset()
    for this_crop_int, data in enumerate(p_list):
        if isinstance(data, type(None)):
            continue
        this_crop_str = incl_vegtypes_str[this_crop_int]
        log(logger, f"   {this_crop_str}...")
        new_var = f"gdd1_{utils.ivt_str2int(this_crop_str)}"
        this_ds = daily_ds.isel(
            patch=np.where(daily_ds.patches1d_itype_veg_str.values == this_crop_str)[0]
        )
        template_da = this_ds.patches1d_itype_veg_str
        this_da = xr.DataArray(
            data=data,
            coords=template_da.coords,
            attrs={"units": "GDD", "long_name": f"{longname_prefix}{this_crop_str}"},
        )

        # Grid this crop
        this_ds["tmp"] = this_da
        da_gridded = grid_one_variable(this_ds, "tmp", vegtype=this_crop_str)
        da_gridded = da_gridded.squeeze(drop=True)

        # Add singleton time dimension and save to output Dataset
        da_gridded = da_gridded.expand_dims(time=dates_rx.time)
        ds_out[new_var] = da_gridded

    return ds_out


def import_and_process_1yr(
    year_1,
    year_n,
    year_index,
    this_year,
    sdates_rx,
    hdates_rx,
    gddaccum_yp_list,
    gddharv_yp_list,
    skip_patches_for_isel_nan_last_year,
    last_year_active_patch_indices_list,
    incorrectly_daily,
    indir,
    incl_vegtypes_str_in,
    h2_ds_file,
    mxmats,
    get_gs_len_da,
    skip_crops,
    outdir_figs,
    logger,
    h1_instantaneous,
):
    """
    Import one year of CLM output data for GDD generation
    """
    save_figs = True
    log(logger, f"netCDF year {this_year}...")
    log(logger, dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # Without dask, this can take a LONG time at resolutions finer than 2-deg
    if importlib_util.find_spec("dask"):
        chunks = {"time": 1}
    else:
        chunks = None

    # Get h1 file (list)
    h1_pattern = os.path.join(indir, "*h1.*.nc")
    h1_filelist = glob.glob(h1_pattern)
    if not h1_filelist:
        h1_pattern = os.path.join(indir, "*h1.*.nc.base")
        h1_filelist = glob.glob(h1_pattern)
        if not h1_filelist:
            error(logger, "No files found matching pattern '*h1.*.nc(.base)'")

    # Get list of crops to include
    if skip_crops is not None:
        crops_to_read = [c for c in utils.define_mgdcrop_list_withgrasses() if c not in skip_crops]
    else:
        crops_to_read = utils.define_mgdcrop_list_withgrasses()

    # Are h1 files instantaneous?
    if h1_instantaneous is None:
        h1_instantaneous = is_instantaneous(xr.open_dataset(h1_filelist[0])["time"])

    if h1_instantaneous:
        slice_year = this_year
    else:
        slice_year = this_year - 1
    dates_ds = import_ds(
        h1_filelist,
        my_vars=["SDATES", "HDATES"],
        my_vegtypes=crops_to_read,
        time_slice=slice(f"{slice_year}-01-01", f"{slice_year}-12-31"),
        chunks=chunks,
    )
    for timestep in dates_ds["time"].values:
        print(timestep)

    if dates_ds.dims["time"] > 1:
        if dates_ds.dims["time"] == 365:
            if not incorrectly_daily:
                log(
                    logger,
                    "   ℹ️ You saved SDATES and HDATES daily, but you only needed annual. Fixing.",
                )
            incorrectly_daily = True
            dates_ds = dates_ds.isel(time=-1)
    else:
        dates_ds = dates_ds.isel(time=0)

    # Make sure NaN masks match
    sdates_all_nan = (
        np.sum(~np.isnan(dates_ds.SDATES.values), axis=dates_ds.SDATES.dims.index("mxsowings")) == 0
    )
    hdates_all_nan = (
        np.sum(~np.isnan(dates_ds.HDATES.values), axis=dates_ds.HDATES.dims.index("mxharvests"))
        == 0
    )
    n_unmatched_nans = np.sum(sdates_all_nan != hdates_all_nan)
    if n_unmatched_nans > 0:
        error(logger, "Output SDATE and HDATE NaN masks do not match.")
    if np.sum(~np.isnan(dates_ds.SDATES.values)) == 0:
        error(logger, "All SDATES are NaN!")

    # Just work with non-NaN patches for now
    skip_patches_for_isel_nan = np.where(sdates_all_nan)[0]
    incl_patches_for_isel_nan = np.where(~sdates_all_nan)[0]
    different_nan_mask = year_index > 0 and not np.array_equal(
        skip_patches_for_isel_nan_last_year, skip_patches_for_isel_nan
    )
    if different_nan_mask:
        log(logger, "   Different NaN mask than last year")
        incl_thisyr_but_nan_lastyr = [
            dates_ds.patch.values[p]
            for p in incl_patches_for_isel_nan
            if p in skip_patches_for_isel_nan_last_year
        ]
    else:
        incl_thisyr_but_nan_lastyr = []
    skipping_patches_for_isel_nan = len(skip_patches_for_isel_nan) > 0
    if skipping_patches_for_isel_nan:
        log(
            logger,
            f"   Ignoring {len(skip_patches_for_isel_nan)} patches with all-NaN sowing and "
            + "harvest dates.",
        )
        dates_incl_ds = dates_ds.isel(patch=incl_patches_for_isel_nan)
    else:
        dates_incl_ds = dates_ds
    incl_patches1d_itype_veg = dates_incl_ds.patches1d_itype_veg

    if year_index == 0:
        incl_vegtypes_str = [c for c in dates_incl_ds.vegtype_str.values if c not in skip_crops]
    else:
        incl_vegtypes_str = incl_vegtypes_str_in
        if isinstance(incl_vegtypes_str, xr.DataArray):
            incl_vegtypes_str = incl_vegtypes_str.values
        if isinstance(incl_vegtypes_str, np.ndarray):
            incl_vegtypes_str = list(incl_vegtypes_str)
        if incl_vegtypes_str != list(dates_incl_ds.vegtype_str.values):
            error(
                logger,
                f"Included veg types differ. Previously {incl_vegtypes_str}, "
                + f"now {dates_incl_ds.vegtype_str.values}",
            )

    if np.sum(~np.isnan(dates_incl_ds.SDATES.values)) == 0:
        error(logger, "All SDATES are NaN after ignoring those patches!")

    # Some patches can have -1 sowing date?? Hopefully just an artifact of me incorrectly saving
    # SDATES/HDATES daily.
    mxsowings = dates_ds.dims["mxsowings"]
    mxsowings_dim = dates_ds.SDATES.dims.index("mxsowings")
    skip_patches_for_isel_sdatelt1 = np.where(dates_incl_ds.SDATES.values < 1)[1]
    skipping_patches_for_isel_sdatelt1 = len(skip_patches_for_isel_sdatelt1) > 0
    if skipping_patches_for_isel_sdatelt1:
        unique_hdates = np.unique(
            dates_incl_ds.HDATES.isel(mxharvests=0, patch=skip_patches_for_isel_sdatelt1).values
        )
        if incorrectly_daily and list(unique_hdates) == [364]:
            log(
                logger,
                f"   ❗ {len(skip_patches_for_isel_sdatelt1)} patches have SDATE < 1, but this"
                + "might have just been because of incorrectly daily outputs. Setting them to 365.",
            )
            new_sdates_ar = dates_incl_ds.SDATES.values
            if mxsowings_dim != 0:
                error(logger, "Code this up")
            new_sdates_ar[0, skip_patches_for_isel_sdatelt1] = 365
            dates_incl_ds["SDATES"] = xr.DataArray(
                data=new_sdates_ar,
                coords=dates_incl_ds["SDATES"].coords,
                attrs=dates_incl_ds["SDATES"].attrs,
            )
        else:
            error(
                logger,
                f"{len(skip_patches_for_isel_sdatelt1)} patches have SDATE < 1. "
                + f"Unique affected hdates: {unique_hdates}",
            )

    # Some patches can have -1 harvest date?? Hopefully just an artifact of me incorrectly saving
    # SDATES/HDATES daily. Can also happen if patch wasn't active last year
    mxharvests = dates_ds.dims["mxharvests"]
    mxharvests_dim = dates_ds.HDATES.dims.index("mxharvests")
    # If a patch was inactive last year but was either (a) harvested the last time it was active or
    # (b) was never active, it will have -1 as its harvest date this year. Such instances are okay.
    hdates_thisyr = dates_incl_ds.HDATES.isel(mxharvests=0)
    skip_patches_for_isel_hdatelt1 = np.where(hdates_thisyr.values < 1)[0]
    skipping_patches_for_isel_hdatelt1 = len(skip_patches_for_isel_hdatelt1) > 0
    if incl_thisyr_but_nan_lastyr and list(skip_patches_for_isel_hdatelt1):
        hdates_thisyr_where_nan_lastyr = hdates_thisyr.sel(patch=incl_thisyr_but_nan_lastyr)
        sdates_thisyr_where_nan_lastyr = dates_incl_ds.SDATES.isel(mxsowings=0).sel(
            patch=incl_thisyr_but_nan_lastyr
        )
        if np.any(hdates_thisyr_where_nan_lastyr < 1):
            new_hdates = dates_incl_ds.HDATES.values
            if mxharvests_dim != 0:
                error(logger, "Code this up")
            patch_list = list(hdates_thisyr.patch.values)
            here = [patch_list.index(x) for x in incl_thisyr_but_nan_lastyr]
            log(
                logger,
                f"   ❗ {len(here)} patches have harvest date -1 because they weren't active last"
                + "year (and were either never active or were harvested when last active). "
                + "Ignoring, but you should have done a run with patches always active if they are "
                + "ever active in the real LU timeseries.",
            )
            new_hdates[0, here] = sdates_thisyr_where_nan_lastyr.values - 1
            dates_incl_ds["HDATES"] = xr.DataArray(
                data=new_hdates,
                coords=dates_incl_ds.HDATES.coords,
                attrs=dates_incl_ds.HDATES.attrs,
            )
            # Recalculate these
            skip_patches_for_isel_hdatelt1 = np.where(
                dates_incl_ds.HDATES.isel(mxharvests=0).values < 1
            )[0]
            skipping_patches_for_isel_hdatelt1 = len(skip_patches_for_isel_hdatelt1) > 0

    # Resolve other issues
    if skipping_patches_for_isel_hdatelt1:
        unique_sdates = np.unique(
            dates_incl_ds.SDATES.isel(patch=skip_patches_for_isel_hdatelt1).values
        )
        if incorrectly_daily and list(unique_sdates) == [1]:
            log(
                logger,
                f"   ❗ {len(skip_patches_for_isel_hdatelt1)} patches have HDATE < 1??? Seems like "
                + "this might have just been because of incorrectly daily outputs; setting them to "
                + "365.",
            )
            new_hdates_ar = dates_incl_ds.HDATES.values
            if mxharvests_dim != 0:
                error(logger, "Code this up")
            new_hdates_ar[0, skip_patches_for_isel_hdatelt1] = 365
            dates_incl_ds["HDATES"] = xr.DataArray(
                data=new_hdates_ar,
                coords=dates_incl_ds["HDATES"].coords,
                attrs=dates_incl_ds["HDATES"].attrs,
            )
        else:
            error(
                logger,
                f"{len(skip_patches_for_isel_hdatelt1)} patches have HDATE < 1. Possible causes:\n"
                + "* Not using constant crop areas (e.g., flanduse_timeseries from "
                + "make_lu_for_gddgen.py)\n   * Not skipping the first 2 years of output\n"
                + f"Unique affected sdates: {unique_sdates}",
            )

    # Make sure there was only one harvest per year
    n_extra_harv = np.sum(
        np.nanmax(
            dates_incl_ds.HDATES.isel(mxharvests=slice(1, mxharvests)).values, axis=mxharvests_dim
        )
        >= 1
    )
    if n_extra_harv > 0:
        error(logger, f"{n_extra_harv} patches have >1 harvest.")

    # Make sure harvest happened the day before sowing
    sdates_clm = dates_incl_ds.SDATES.values.squeeze()
    hdates_clm = dates_incl_ds.HDATES.isel(mxharvests=0).values
    diffdates_clm = sdates_clm - hdates_clm
    diffdates_clm[(sdates_clm == 1) & (hdates_clm == 365)] = 1
    if list(np.unique(diffdates_clm)) != [1]:
        error(logger, f"Not all sdates-hdates are 1: {np.unique(diffdates_clm)}")

    # Import expected sowing dates. This will also be used as our template output file.
    imported_sdates = isinstance(sdates_rx, str)
    sdates_rx = import_rx_dates("s", sdates_rx, incl_patches1d_itype_veg, mxsowings, logger)
    check_sdates(dates_incl_ds, sdates_rx, outdir_figs, logger)

    # Import hdates, if needed
    imported_hdates = isinstance(hdates_rx, str)
    hdates_rx_orig = import_rx_dates(
        "h", hdates_rx, incl_patches1d_itype_veg, mxsowings, logger
    )  # Yes, mxsowings even when importing harvests

    # Limit growing season to CLM max growing season length, if needed
    if mxmats and (imported_sdates or imported_hdates):
        print("   Limiting growing season length...")
        hdates_rx = hdates_rx_orig.copy()
        for var in hdates_rx_orig:
            if var == "time_bounds":
                continue

            # Get max growing season length
            vegtype_int = int(
                var.split("_")[1]
            )  # netCDF variable name v should be something like gs1_17
            vegtype_str = utils.ivt_int2str(vegtype_int)
            if vegtype_str == "soybean":
                vegtype_str = "temperate_soybean"
            elif vegtype_str == "irrigated_soybean":
                vegtype_str = "irrigated_temperate_soybean"

            mxmat = mxmats[vegtype_str]
            if np.isinf(mxmat):
                print(f"      Not limiting {vegtype_str}: No mxmat value")
                continue

            # Get "prescribed" growing season length
            gs_len_rx_da = get_gs_len_da(hdates_rx_orig[var] - sdates_rx[var])
            not_ok = gs_len_rx_da.values > mxmat
            if not np.any(not_ok):
                print(f"      Not limiting {vegtype_str}: No rx season > {mxmat} days")
                continue

            hdates_limited = hdates_rx_orig[var].copy().values
            hdates_limited[np.where(not_ok)] = sdates_rx[var].values[np.where(not_ok)] + mxmat
            hdates_limited[np.where(hdates_limited > 365)] -= 365
            if np.any(hdates_limited < 1):
                raise RuntimeError("Limited hdates < 1")
            if np.any(hdates_limited > 365):
                raise RuntimeError("Limited hdates > 365")
            hdates_rx[var] = xr.DataArray(
                data=hdates_limited,
                coords=hdates_rx_orig[var].coords,
                attrs=hdates_rx_orig[var].attrs,
            )
            print(
                f"      Limited {vegtype_str} growing season length to {mxmat}. Longest was "
                + f"{int(np.max(gs_len_rx_da.values))}, now "
                + f"{int(np.max(get_gs_len_da(hdates_rx[var] - sdates_rx[var]).values))}."
            )
    else:
        hdates_rx = hdates_rx_orig

    log(logger, "   Importing accumulated GDDs...")
    clm_gdd_var = "GDDACCUM"
    my_vars = [clm_gdd_var, "GDDHARV"]
    patterns = [f"*h2.{this_year-1}-01*.nc", f"*h2.{this_year-1}-01*.nc.base"]
    for pat in patterns:
        pattern = os.path.join(indir, pat)
        h2_files = glob.glob(pattern)
        if h2_files:
            break
    if not h2_files:
        error(logger, f"No files found matching patterns: {patterns}")
    h2_ds = import_ds(
        h2_files,
        my_vars=my_vars,
        my_vegtypes=crops_to_read,
        chunks=chunks,
    )

    # Restrict to patches we're including
    if skipping_patches_for_isel_nan:
        if not np.array_equal(dates_ds.patch.values, h2_ds.patch.values):
            error(logger, "dates_ds and h2_ds don't have the same patch list!")
        h2_incl_ds = h2_ds.isel(patch=incl_patches_for_isel_nan)
    else:
        h2_incl_ds = h2_ds

    if not np.any(h2_incl_ds[clm_gdd_var].values != 0):
        error(logger, f"All {clm_gdd_var} values are zero!")

    # Get standard datetime axis for outputs
    n_years = year_n - year_1 + 1

    if len(gddaccum_yp_list) == 0:
        last_year_active_patch_indices_list = [None for vegtype_str in incl_vegtypes_str]
        gddaccum_yp_list = [None for vegtype_str in incl_vegtypes_str]
        if save_figs:
            gddharv_yp_list = [None for vegtype_str in incl_vegtypes_str]

    incl_vegtype_indices = []
    for var, vegtype_str in enumerate(incl_vegtypes_str):
        if vegtype_str in skip_crops:
            log(logger, f"      SKIPPING {vegtype_str}")
            continue

        vegtype_int = utils.vegtype_str2int(vegtype_str)[0]
        this_crop_full_patchlist = list(xr_flexsel(h2_ds, vegtype=vegtype_str).patch.values)

        # Get time series for each patch of this type
        this_crop_ds = xr_flexsel(h2_incl_ds, vegtype=vegtype_str)
        this_crop_gddaccum_da = this_crop_ds[clm_gdd_var]
        check_gddharv = False
        if save_figs:
            this_crop_gddharv_da = this_crop_ds["GDDHARV"]
            check_gddharv = True
        if not this_crop_gddaccum_da.size:
            continue
        log(logger, f"      {vegtype_str}...")
        incl_vegtype_indices = incl_vegtype_indices + [var]

        # Get prescribed harvest dates for these patches
        lon_points = this_crop_ds.patches1d_lon.values
        lat_points = this_crop_ds.patches1d_lat.values
        this_crop_hdates_rx = this_crop_map_to_patches(
            lon_points, lat_points, hdates_rx, vegtype_int
        )

        if isinstance(gddaccum_yp_list[var], type(None)):
            gddaccum_yp_list[var] = np.full((n_years + 1, len(this_crop_full_patchlist)), np.nan)
            if save_figs:
                gddharv_yp_list[var] = np.full((n_years + 1, len(this_crop_full_patchlist)), np.nan)

        # Get the accumulated GDDs at each prescribed harvest date
        gddaccum_atharv_p = np.full(this_crop_hdates_rx.shape, np.nan)
        if save_figs:
            gddharv_atharv_p = np.full(this_crop_hdates_rx.shape, np.nan)
        unique_rx_hdates = np.unique(this_crop_hdates_rx.values)
        # Build an indexing tuple
        patches = []
        i_patches = []
        i_times = []
        for hdate in unique_rx_hdates:
            here = np.where(this_crop_hdates_rx.values == hdate)[0]
            patches += list(this_crop_gddaccum_da.patch.values[here])
            i_patches += list(here)
            i_times += list(np.full((len(here),), int(hdate - 1)))
        # Sort back to correct order
        if not np.all(
            this_crop_gddaccum_da.patch.values[:-1] <= this_crop_gddaccum_da.patch.values[1:]
        ):
            error(logger, "This code depends on DataArray patch list being sorted.")
        sortorder = np.argsort(patches)
        i_patches = list(np.array(i_patches)[np.array(sortorder)])
        i_times = list(np.array(i_times)[np.array(sortorder)])
        # Select using the indexing tuple
        gddaccum_atharv_p = this_crop_gddaccum_da.values[(i_times, i_patches)]
        gddharv_atharv_p = None
        if save_figs:
            gddharv_atharv_p = this_crop_gddharv_da.values[(i_times, i_patches)]
        if np.any(np.isnan(gddaccum_atharv_p)):
            log(
                logger,
                f"         ❗ {np.sum(np.isnan(gddaccum_atharv_p))}/{len(gddaccum_atharv_p)} "
                + "NaN after extracting GDDs accumulated at harvest",
            )
        if save_figs and gddharv_atharv_p is not None and np.any(np.isnan(gddharv_atharv_p)):
            if np.all(np.isnan(gddharv_atharv_p)):
                log(logger, "         ❗ All GDDHARV are NaN; should only affect figure")
                check_gddharv = False
            else:
                log(
                    logger,
                    f"         ❗ {np.sum(np.isnan(gddharv_atharv_p))}/{len(gddharv_atharv_p)} "
                    + "NaN after extracting GDDHARV",
                )

        # Assign these to growing seasons based on whether gs crossed new year
        this_year_active_patch_indices = [
            this_crop_full_patchlist.index(x) for x in this_crop_ds.patch.values
        ]
        this_crop_sdates_rx = this_crop_map_to_patches(
            lon_points, lat_points, sdates_rx, vegtype_int
        )
        where_gs_thisyr = np.where(this_crop_sdates_rx < this_crop_hdates_rx)[0]
        tmp_gddaccum = np.full(this_crop_sdates_rx.shape, np.nan)
        tmp_gddaccum[where_gs_thisyr] = gddaccum_atharv_p[where_gs_thisyr]
        if save_figs:
            tmp_gddharv = np.full(tmp_gddaccum.shape, np.nan)
            tmp_gddharv[where_gs_thisyr] = gddharv_atharv_p[where_gs_thisyr]
        if year_index > 0:
            last_year_active_patch_indices = last_year_active_patch_indices_list[var]
            where_gs_lastyr = np.where(this_crop_sdates_rx > this_crop_hdates_rx)[0]
            active_this_year_where_gs_lastyr_indices = [
                this_year_active_patch_indices[x] for x in where_gs_lastyr
            ]
            if not np.array_equal(last_year_active_patch_indices, this_year_active_patch_indices):
                if incorrectly_daily:
                    log(
                        logger,
                        "         ❗ This year's active patch indices differ from last year's. "
                        + "Allowing because this might just be an artifact of incorrectly daily "
                        + "outputs, BUT RESULTS MUST NOT BE TRUSTED.",
                    )
                else:
                    error(logger, "This year's active patch indices differ from last year's.")
            # Make sure we're not about to overwrite any existing values.
            if np.any(
                ~np.isnan(
                    gddaccum_yp_list[var][year_index - 1, active_this_year_where_gs_lastyr_indices]
                )
            ):
                if incorrectly_daily:
                    log(
                        logger,
                        "         ❗ Unexpected non-NaN for last season's GDD accumulation. "
                        + "Allowing because this might just be an artifact of incorrectly daily "
                        + "outputs, BUT RESULTS MUST NOT BE TRUSTED.",
                    )
                else:
                    error(logger, "Unexpected non-NaN for last season's GDD accumulation")
            if save_figs and np.any(
                ~np.isnan(
                    gddharv_yp_list[var][year_index - 1, active_this_year_where_gs_lastyr_indices]
                )
            ):
                if incorrectly_daily:
                    log(
                        logger,
                        "         ❗ Unexpected non-NaN for last season's GDDHARV. Allowing "
                        + "because this might just be an artifact of incorrectly daily outputs, "
                        + "BUT RESULTS MUST NOT BE TRUSTED.",
                    )
                else:
                    error(logger, "Unexpected non-NaN for last season's GDDHARV")
            # Fill.
            gddaccum_yp_list[var][year_index - 1, active_this_year_where_gs_lastyr_indices] = (
                gddaccum_atharv_p[where_gs_lastyr]
            )
            if save_figs:
                gddharv_yp_list[var][year_index - 1, active_this_year_where_gs_lastyr_indices] = (
                    gddharv_atharv_p[where_gs_lastyr]
                )
            # Last year's season should be filled out now; make sure.
            if np.any(
                np.isnan(
                    gddaccum_yp_list[var][year_index - 1, active_this_year_where_gs_lastyr_indices]
                )
            ):
                if incorrectly_daily:
                    log(
                        logger,
                        "         ❗ Unexpected NaN for last season's GDD accumulation. Allowing "
                        + "because this might just be an artifact of incorrectly daily outputs, "
                        + "BUT RESULTS MUST NOT BE TRUSTED.",
                    )
                else:
                    error(logger, "Unexpected NaN for last season's GDD accumulation.")
            if (
                save_figs
                and check_gddharv
                and np.any(
                    np.isnan(
                        gddharv_yp_list[var][
                            year_index - 1, active_this_year_where_gs_lastyr_indices
                        ]
                    )
                )
            ):
                if incorrectly_daily:
                    log(
                        logger,
                        "         ❗ Unexpected NaN for last season's GDDHARV. Allowing because "
                        + "this might just be an artifact of incorrectly daily outputs, BUT "
                        + "RESULTS MUST NOT BE TRUSTED.",
                    )
                else:
                    error(logger, "Unexpected NaN for last season's GDDHARV.")
        gddaccum_yp_list[var][year_index, this_year_active_patch_indices] = tmp_gddaccum
        if save_figs:
            gddharv_yp_list[var][year_index, this_year_active_patch_indices] = tmp_gddharv

        # Make sure that NaN masks are the same for this year's sdates and 'filled-out' GDDs from
        # last year
        if year_index > 0:
            nanmask_output_sdates = np.isnan(
                dates_ds.SDATES.isel(
                    mxsowings=0, patch=np.where(dates_ds.patches1d_itype_veg_str == vegtype_str)[0]
                ).values
            )
            nanmask_output_gdds_lastyr = np.isnan(gddaccum_yp_list[var][year_index - 1, :])
            if not np.array_equal(nanmask_output_gdds_lastyr, nanmask_output_sdates):
                if incorrectly_daily:
                    log(
                        logger,
                        "         ❗ NaN masks differ between this year's sdates and 'filled-out' "
                        + "GDDs from last year. Allowing because this might just be an artifact of "
                        + "incorrectly daily outputs, BUT RESULTS MUST NOT BE TRUSTED.",
                    )
                else:
                    error(
                        logger,
                        "NaN masks differ between this year's sdates and 'filled-out' GDDs from "
                        + "last year",
                    )
        last_year_active_patch_indices_list[var] = this_year_active_patch_indices

    skip_patches_for_isel_nan_last_year = skip_patches_for_isel_nan

    # Could save space by only saving variables needed for gridding
    log(logger, "   Saving h2_ds...")
    h2_ds.to_netcdf(h2_ds_file)

    return (
        h2_ds,
        sdates_rx,
        hdates_rx,
        gddaccum_yp_list,
        gddharv_yp_list,
        skip_patches_for_isel_nan_last_year,
        last_year_active_patch_indices_list,
        incorrectly_daily,
        incl_vegtypes_str,
        incl_patches1d_itype_veg,
        mxsowings,
        h1_instantaneous,
    )


def get_multicrop_maps(this_ds, these_vars, crop_fracs_yx, dummy_fill, gdd_units):
    # pylint: disable=missing-function-docstring
    # Get GDDs for these crops
    da_each_cft = xr.concat((this_ds[x] for i, x in enumerate(these_vars)), dim="cft")
    if "time" in this_ds.dims:
        da_each_cft = da_each_cft.isel(time=0, drop=True)
    da_each_cft = da_each_cft.where(da_each_cft != dummy_fill)
    da_each_cft.attrs["units"] = gdd_units

    # What are the maximum differences seen between different crop types?
    if len(these_vars) > 1:
        max_diff = np.nanmax(da_each_cft.max(dim="cft") - da_each_cft.min(dim="cft"))
        if max_diff > 0:
            print(f"   Max difference among crop types: {np.round(max_diff)}")

    if crop_fracs_yx is None:
        return da_each_cft.isel(cft=0, drop=True)

    # Warn if GDD is NaN anywhere that there is area
    da_each_cft["cft"] = crop_fracs_yx["cft"]
    gdd_nan_area_pos = np.isnan(da_each_cft) & (crop_fracs_yx > 0)
    if np.any(gdd_nan_area_pos):
        total_bad_croparea = np.nansum(crop_fracs_yx.where(gdd_nan_area_pos).values)
        total_croparea = np.nansum(crop_fracs_yx.values)
        print(
            "   GDD reqt NaN but area positive "
            + f"({np.round(total_bad_croparea/total_croparea*100, 1)}% of this crop's area)"
        )

    # Get areas and weights, masking cell-crops with NaN GDDs
    crop_fracs_yx = crop_fracs_yx.where(~np.isnan(da_each_cft))
    crop_area_yx = crop_fracs_yx.sum(dim="cft")
    weights_yx = crop_fracs_yx / crop_area_yx
    weights_sum_gt0 = weights_yx.sum(dim="cft").where(weights_yx > 0)
    assert np.isclose(np.nanmin(weights_sum_gt0.values), 1.0)
    assert np.isclose(np.nanmax(weights_sum_gt0.values), 1.0)

    # Mask GDDs and weights where there is no area
    da_each_cft = da_each_cft.where(crop_fracs_yx > 0)
    if len(these_vars) == 1:
        return da_each_cft.isel(cft=0, drop=True)
    weights_yx = weights_yx.where(crop_fracs_yx > 0)
    weights_sum = weights_yx.sum(dim="cft").where(crop_area_yx > 0)
    assert np.isclose(np.nanmin(weights_sum.values), 1.0)
    assert np.isclose(np.nanmax(weights_sum.values), 1.0)

    # Ensure grid match between GDDs and weights
    if not np.array_equal(da_each_cft["lon"].values, weights_yx["lon"].values):
        raise RuntimeError("lon mismatch")
    if not np.array_equal(da_each_cft["lat"].values, weights_yx["lat"].values):
        raise RuntimeError("lat mismatch")

    # Get area-weighted mean GDD requirements for all crops
    this_da = (da_each_cft * weights_yx).sum(dim="cft")
    this_da.attrs["units"] = gdd_units
    this_da = this_da.where(crop_area_yx > 0)

    # Ensure that weighted mean is between each cell's min and max
    where_bad = (this_da < da_each_cft.min(dim="cft")) | (this_da > da_each_cft.max(dim="cft"))
    if np.any(where_bad):
        where_below_min = this_da.where(this_da < da_each_cft.min(dim="cft"))
        worst_below_min = np.min((da_each_cft.min(dim="cft") - where_below_min).values)
        where_above_max = this_da.where(this_da > da_each_cft.max(dim="cft"))
        worst_above_max = np.max((where_above_max - da_each_cft.max(dim="cft")).values)
        worst = max(worst_below_min, worst_above_max)
        tol = 1e-12
        if worst > 1e-12:
            raise RuntimeError(
                f"Some value is outside expected range by {worst} (exceeds tolerance {tol})"
            )

    return this_da


if CAN_PLOT:

    def get_bounds_ncolors(gdd_spacing, diff_map_yx):
        """
        Get information about color bar
        """
        vmax = np.floor(np.nanmax(diff_map_yx.values) / gdd_spacing) * gdd_spacing
        vmin = -vmax
        epsilon = np.nextafter(0, 1)
        bounds = list(np.arange(vmin, vmax, gdd_spacing)) + [vmax - epsilon]
        if 0 in bounds:
            bounds.remove(0)
            bounds[bounds.index(-gdd_spacing)] /= 2
            bounds[bounds.index(gdd_spacing)] /= 2
        n_colors = len(bounds) + 1
        return vmax, bounds, n_colors

    def make_gengdd_map(
        this_axis,
        this_map,
        this_title,
        vmax,
        bin_width,
        fontsize_ticklabels,
        fontsize_titles,
        bounds=None,
        extend="both",
        cmap=None,
        cbar_ticks=None,
        vmin=None,
    ):
        """
        Make maps
        """
        if bounds:
            if not cmap:
                raise RuntimeError("Calling make_map() with bounds requires cmap to be specified")
            norm = mcolors.BoundaryNorm(bounds, cmap.N, extend=extend)
            im1 = this_axis.pcolormesh(
                this_map.lon.values,
                this_map.lat.values,
                this_map,
                shading="auto",
                norm=norm,
                cmap=cmap,
            )
        else:
            if np.any(this_map.values < 0):
                gdd_spacing = 500
                vmax = np.floor(np.nanmax(this_map.values) / gdd_spacing) * gdd_spacing
                if vmin is not None:
                    raise RuntimeError("Do not specify vmin in this call of make_map()")
                vmin = -vmax
                n_colors = vmax / gdd_spacing
                if n_colors % 2 == 0:
                    n_colors += 1
                if not cmap:
                    cmap = cm.get_cmap(cropcal_colors["div_other_nonnorm"], n_colors)

                if np.any(this_map.values > vmax) and np.any(this_map.values < vmin):
                    extend = "both"
                elif np.any(this_map.values > vmax):
                    extend = "max"
                elif np.any(this_map.values < vmin):
                    extend = "min"
                else:
                    extend = "neither"

            else:
                if vmin is None:
                    vmin = 0
                else:
                    vmin = np.floor(vmin / 500) * 500
                vmax = np.floor(vmax / 500) * 500
                n_colors = int(vmax / 500)
                if not cmap:
                    cmap = cm.get_cmap(cropcal_colors["seq_other"], n_colors + 1)
                extend = "max"
                extend_color = cmap.colors[-1]
                cmap = mcolors.ListedColormap(cmap.colors[:n_colors])
                cmap.set_over(extend_color)

            im1 = this_axis.pcolormesh(
                this_map.lon.values,
                this_map.lat.values,
                this_map,
                shading="auto",
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
            )

        this_axis.set_extent([-180, 180, -63, 90], crs=ccrs.PlateCarree())
        this_axis.coastlines(linewidth=0.3)
        this_axis.set_title(this_title, fontsize=fontsize_titles, fontweight="bold", y=0.96)
        cbar = plt.colorbar(
            im1,
            orientation="horizontal",
            fraction=0.1,
            pad=0.02,
            aspect=40,
            extend=extend,
            spacing="proportional",
        )
        cbar.ax.tick_params(labelsize=fontsize_ticklabels)
        cbar.ax.set_xlabel(this_map.attrs["units"], fontsize=fontsize_ticklabels)
        cbar.ax.xaxis.set_label_coords(x=0.115, y=2.6)
        if cbar_ticks:
            cbar.ax.set_xticks(cbar_ticks)

        ticks = np.arange(-60, 91, bin_width)
        ticklabels = [str(x) for x in ticks]
        for i, tick in enumerate(ticks):
            if tick % 2:
                ticklabels[i] = ""
        plt.yticks(np.arange(-60, 91, 15), labels=ticklabels, fontsize=fontsize_ticklabels)
        plt.axis("off")

    def get_non_nans(in_da, fill_value):
        """
        Get non-NaN, non-fill values of a DataArray
        """
        in_da = in_da.where(in_da != fill_value)
        return in_da.values[~np.isnan(in_da.values)]

    def set_boxplot_props(bpl, color, linewidth):
        """
        Set boxplot properties
        """
        linewidth = 1.5
        plt.setp(bpl["boxes"], color=color, linewidth=linewidth)
        plt.setp(bpl["whiskers"], color=color, linewidth=linewidth)
        plt.setp(bpl["caps"], color=color, linewidth=linewidth)
        plt.setp(bpl["medians"], color=color, linewidth=linewidth)
        plt.setp(
            bpl["fliers"],
            markeredgecolor=color,
            markersize=6,
            linewidth=linewidth,
            markeredgewidth=linewidth / 2,
        )

    def make_plot(data, offset, linewidth):
        """
        Make boxplot
        """
        offset = 0.4 * offset
        bpl = plt.boxplot(
            data,
            positions=np.array(range(len(data))) * 2.0 + offset,
            widths=0.6,
            boxprops={"linewidth": linewidth},
            whiskerprops={"linewidth": linewidth},
            capprops={"linewidth": linewidth},
            medianprops={"linewidth": linewidth},
            flierprops={"markeredgewidth": 0.5},
        )
        return bpl

    def make_figures(
        first_land_use_year,
        last_land_use_year,
        land_use_file,
        run1_name,
        run2_name,
        logger,
        this_dir=None,
        gdd_maps_ds=None,
        gddharv_maps_ds=None,
        outdir_figs=None,
        linewidth=1.5,
    ):
        """
        Make map-and-boxplot figures
        """
        if not gdd_maps_ds:
            if not this_dir:
                error(
                    logger,
                    "If not providing gdd_maps_ds, you must provide thisDir (location of "
                    + "gdd_maps.nc)",
                )
            gdd_maps_ds = xr.open_dataset(this_dir + "gdd_maps.nc")
        if not gddharv_maps_ds:
            if not this_dir:
                error(
                    logger,
                    "If not providing gddharv_maps_ds, you must provide thisDir (location of "
                    + "gddharv_maps.nc)",
                )
            gddharv_maps_ds = xr.open_dataset(this_dir + "gdd_maps.nc")

        # Get info
        incl_vegtypes_str = gdd_maps_ds.attrs["incl_vegtypes_str"]
        if incl_vegtypes_str is None:
            incl_vegtypes_str = []
        elif isinstance(incl_vegtypes_str, np.ndarray):
            incl_vegtypes_str = list(incl_vegtypes_str)
        dummy_fill = gdd_maps_ds.attrs["dummy_fill"]
        if not outdir_figs:
            outdir_figs = gdd_maps_ds.attrs["outdir_figs"]
        try:
            year_1 = gdd_maps_ds.attrs["y1"]
            year_n = gdd_maps_ds.attrs["yN"]
        # Backwards compatibility with a bug (fixed 2023-01-03)
        except KeyError:
            year_1 = gdd_maps_ds.attrs["first_season"]
            year_n = gdd_maps_ds.attrs["last_season"]
        # Import LU data, if doing so
        if land_use_file:
            year_1_lu = year_1 if first_land_use_year is None else first_land_use_year
            year_n_lu = year_n if last_land_use_year is None else last_land_use_year
            lu_ds = cc.open_lu_ds(land_use_file, year_1_lu, year_n_lu, gdd_maps_ds, ungrid=False)
            lu_years_text = f" (masked by {year_1_lu}-{year_n_lu} area)"
            lu_years_file = f"_mask{year_1_lu}-{year_n_lu}"
        else:
            lu_ds = None
            lu_years_text = ""
            lu_years_file = ""

        # layout = "3x1"
        # layout = "2x2"
        layout = "3x2"
        bin_width = 15
        lat_bin_edges = np.arange(0, 91, bin_width)

        fontsize_titles = 12
        fontsize_axislabels = 12
        fontsize_ticklabels = 12

        n_bins = len(lat_bin_edges) - 1
        bin_names = ["All"]
        for this_bin in np.arange(n_bins):
            lower = lat_bin_edges[this_bin]
            upper = lat_bin_edges[this_bin + 1]
            bin_names.append(f"{lower}–{upper}")

        color_old = cropcal_colors_cases(run1_name)
        if color_old is None:
            color_old = "#beaed4"
        color_new = cropcal_colors_cases(run2_name)
        if color_new is None:
            color_new = "#7fc97f"
        gdd_units = "GDD (°C • day)"

        # Maps
        nplot_y = 3
        nplot_x = 1
        log(logger, "Making before/after maps...")
        vegtype_list = incl_vegtypes_str
        if land_use_file:
            vegtype_list += ["Corn", "Cotton", "Rice", "Soybean", "Sugarcane", "Wheat"]
        for vegtype_str in vegtype_list:
            print(f"{vegtype_str}...")

            # Get component types
            if vegtype_str in incl_vegtypes_str:
                vegtypes_str = [vegtype_str]
            elif not lu_ds:
                raise RuntimeError(f"If mapping {vegtype_str}, you must provide land use dataset")
            else:
                vegtypes_str = [x for x in incl_vegtypes_str if vegtype_str.lower() in x]
            vegtypes_int = [utils.vegtype_str2int(x)[0] for x in vegtypes_str]

            # Crop fraction map (for masking and weighting)
            if lu_ds:
                crop_fracs_yx = (
                    lu_ds.LANDFRAC_PFT * lu_ds.PCT_CROP * lu_ds.PCT_CFT.sel(cft=vegtypes_int)
                ).sum(dim="time")
                if np.sum(crop_fracs_yx) == 0:
                    print(f"Skipping {vegtype_str} (no area)")
                    continue
            else:
                crop_fracs_yx = None

            these_vars = [f"gdd1_{x}" for x in vegtypes_int]
            gddharv_map_yx = get_multicrop_maps(
                gddharv_maps_ds, these_vars, crop_fracs_yx, dummy_fill, gdd_units
            )
            gdd_map_yx = get_multicrop_maps(
                gdd_maps_ds, these_vars, crop_fracs_yx, dummy_fill, gdd_units
            )

            # Get figure title
            if len(vegtypes_str) > 1:
                vegtype_str_title = vegtype_str
            else:
                vegtype_str_title = vegtype_str.replace("_", " ")
                if "irrigated" not in vegtype_str:
                    vegtype_str_title = "rainfed " + vegtype_str_title
                vegtype_str_title = vegtype_str_title.capitalize()

            vmin = min(np.min(gdd_map_yx), np.min(gddharv_map_yx)).values
            vmax = max(np.max(gdd_map_yx), np.max(gddharv_map_yx)).values

            # Set up figure and first subplot
            this_axis = None
            if layout == "3x1":
                fig = plt.figure(figsize=(7.5, 14))
                this_axis = fig.add_subplot(nplot_y, nplot_x, 1, projection=ccrs.PlateCarree())
            elif layout == "2x2":
                fig = plt.figure(figsize=(12, 6))
                spec = fig.add_gridspec(nrows=2, ncols=2, width_ratios=[0.4, 0.6])
                this_axis = fig.add_subplot(spec[0, 0], projection=ccrs.PlateCarree())
            elif layout == "3x2":
                fig = plt.figure(figsize=(14, 9))
                spec = fig.add_gridspec(nrows=3, ncols=2, width_ratios=[0.5, 0.5], wspace=0.2)
                this_axis = fig.add_subplot(spec[0, 0], projection=ccrs.PlateCarree())
            else:
                error(logger, f"layout {layout} not recognized")

            gddharv_all_nan = np.all(np.isnan(gddharv_map_yx.values))
            if gddharv_all_nan:
                this_title = f"{run1_name} (GDDHARV all NaN?)"
            else:
                this_min = int(np.round(np.nanmin(gddharv_map_yx)))
                this_max = int(np.round(np.nanmax(gddharv_map_yx)))
                this_title = f"{run1_name} (range {this_min}–{this_max})"
            make_gengdd_map(
                this_axis,
                gddharv_map_yx,
                this_title,
                vmax,
                bin_width,
                fontsize_ticklabels,
                fontsize_titles,
                vmin=vmin,
            )

            if layout == "3x1":
                this_axis = fig.add_subplot(nplot_y, nplot_x, 2, projection=ccrs.PlateCarree())
            elif layout in ["2x2", "3x2"]:
                this_axis = fig.add_subplot(spec[1, 0], projection=ccrs.PlateCarree())
            else:
                error(logger, f"layout {layout} not recognized")
            this_min = int(np.round(np.nanmin(gdd_map_yx)))
            this_max = int(np.round(np.nanmax(gdd_map_yx)))
            this_title = f"{run2_name} (range {this_min}–{this_max})"
            make_gengdd_map(
                this_axis,
                gdd_map_yx,
                this_title,
                vmax,
                bin_width,
                fontsize_ticklabels,
                fontsize_titles,
                vmin=vmin,
            )

            # Difference
            if not gddharv_all_nan and layout == "3x2":
                this_axis = fig.add_subplot(spec[2, 0], projection=ccrs.PlateCarree())
                this_min = int(np.round(np.nanmin(gdd_map_yx)))
                this_max = int(np.round(np.nanmax(gdd_map_yx)))
                this_title = f"{run2_name} minus {run1_name}"
                diff_map_yx = gdd_map_yx - gddharv_map_yx
                diff_map_yx.attrs["units"] = gdd_units

                gdd_spacing = 500
                vmax, bounds, n_colors = get_bounds_ncolors(gdd_spacing, diff_map_yx)
                if n_colors < 9:
                    gdd_spacing = 250
                    vmax, bounds, n_colors = get_bounds_ncolors(gdd_spacing, diff_map_yx)

                cmap = cm.get_cmap(cropcal_colors["div_other_nonnorm"], n_colors)
                cbar_ticks = []
                include_0bin_ticks = n_colors <= 13
                if vmax <= 3000:
                    tick_spacing = gdd_spacing * 2
                elif vmax <= 5000:
                    tick_spacing = 1500
                else:
                    tick_spacing = 2000
                previous = -np.inf
                for bound in bounds:
                    if (not include_0bin_ticks) and (previous < 0 < bound):
                        cbar_ticks.append(0)
                    if bound % tick_spacing == 0 or (
                        include_0bin_ticks and abs(bound) == gdd_spacing / 2
                    ):
                        cbar_ticks.append(bound)
                    previous = bound

                make_gengdd_map(
                    this_axis,
                    diff_map_yx,
                    this_title,
                    vmax,
                    bin_width,
                    fontsize_ticklabels,
                    fontsize_titles,
                    bounds=bounds,
                    extend="both",
                    cmap=cmap,
                    cbar_ticks=cbar_ticks,
                )

            # Boxplots #####################

            gdd_vector = get_non_nans(gdd_map_yx, dummy_fill)
            gddharv_vector = get_non_nans(gddharv_map_yx, dummy_fill)

            lat_abs = np.abs(gdd_map_yx.lat.values)
            gdd_bybin_old = [gddharv_vector]
            gdd_bybin_new = [gdd_vector]
            for this_bin in np.arange(n_bins):
                lower = lat_bin_edges[this_bin]
                upper = lat_bin_edges[this_bin + 1]
                lat_inds = np.where((lat_abs >= lower) & (lat_abs < upper))[0]
                this_bin_gdd_vector = get_non_nans(gdd_map_yx[lat_inds, :], dummy_fill)
                this_bin_gddharv_vector = get_non_nans(gddharv_map_yx[lat_inds, :], dummy_fill)
                gdd_bybin_old.append(this_bin_gddharv_vector)
                gdd_bybin_new.append(this_bin_gdd_vector)

            if layout == "3x1":
                this_axis = fig.add_subplot(nplot_y, nplot_x, 3)
            elif layout in ["2x2", "3x2"]:
                this_axis = fig.add_subplot(spec[:, 1])
            else:
                error(logger, f"layout {layout} not recognized")

            # Shift bottom of plot up to make room for legend
            ax_pos = this_axis.get_position()
            this_axis.set_position(Bbox.from_extents(ax_pos.x0, 0.19, ax_pos.x1, ax_pos.y1))
            # Define legend position
            legend_bbox_to_anchor = (0, -0.15, 1, 0.2)

            bpl = make_plot(gdd_bybin_old, -1, linewidth)
            bpr = make_plot(gdd_bybin_new, 1, linewidth)
            set_boxplot_props(bpl, color_old, linewidth)
            set_boxplot_props(bpr, color_new, linewidth)

            # draw temporary lines to create a legend
            plt.plot([], c=color_old, label=run1_name, linewidth=linewidth)
            plt.plot([], c=color_new, label=run2_name, linewidth=linewidth)
            plt.legend(
                fontsize=fontsize_titles,
                bbox_to_anchor=legend_bbox_to_anchor,
                ncol=2,
                loc="lower left",
                mode="expand",
            )

            plt.xticks(range(0, len(bin_names) * 2, 2), bin_names, fontsize=fontsize_ticklabels)
            plt.yticks(fontsize=fontsize_ticklabels)
            this_axis.spines["right"].set_visible(False)
            this_axis.spines["top"].set_visible(False)

            plt.xlabel("Latitude zone (absolute value)", fontsize=fontsize_axislabels)
            plt.ylabel(gdd_units, fontsize=fontsize_axislabels)
            this_axis.yaxis.set_label_coords(-0.11, 0.5)
            plt.title("Zonal changes", fontsize=fontsize_titles, fontweight="bold")

            plt.suptitle(
                f"Maturity requirements: {vegtype_str_title}" + lu_years_text,
                fontsize=fontsize_titles * 1.2,
                fontweight="bold",
                y=0.95,
            )

            if vegtype_str in incl_vegtypes_str:
                outfile = os.path.join(
                    outdir_figs,
                    f"{these_vars[0]}_{vegtype_str}_gs{year_1}-{year_n}{lu_years_file}.png",
                )
            else:
                outfile = os.path.join(
                    outdir_figs, f"{vegtype_str}_gs{year_1}-{year_n}{lu_years_file}.png"
                )
            plt.savefig(outfile, dpi=300, transparent=False, facecolor="white", bbox_inches="tight")
            plt.close()

        log(logger, "Done.")
