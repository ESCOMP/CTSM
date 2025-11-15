"""
Generate maturity requirements (GDD) from outputs of a GDD-generating run
"""

import os
import sys
import pickle
import datetime as dt
import argparse
import logging
import numpy as np
import xarray as xr

# Import the CTSM Python utilities.
# sys.path.insert() is necessary for RXCROPMATURITY to work. The fact that it's calling this script
# in the RUN phase seems to require the python/ directory to be manually added to path.
_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)
from ctsm.ctsm_logging import log, error  # pylint: disable=wrong-import-position
import ctsm.crop_calendars.cropcal_module as cc  # pylint: disable=wrong-import-position
import ctsm.crop_calendars.generate_gdds_functions as gddfn  # pylint: disable=wrong-import-position
from ctsm.crop_calendars.import_ds import (  # pylint: disable=wrong-import-position
    get_files_in_time_slice,  # pylint: disable=wrong-import-position
)  # pylint: disable=wrong-import-position

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments


def _get_max_growing_season_lengths(max_season_length_from_hdates_file, paramfile, cushion):
    """
    Import maximum growing season lengths from paramfile, if doing so.
    """
    if max_season_length_from_hdates_file:
        return None

    mxmats = cc.import_max_gs_length(paramfile)

    if cushion:
        mxmats = cc.cushion_gs_length(mxmats, cushion)

    return mxmats


def _get_history_yr_range(first_season, last_season):
    """
    Get a range object that can be used for looping over all years we need to process timestamps
    from.
    """
    # Saving at the end of a year receive the timestamp of the END of the year's final timestep,
    # which means it will actually be 00:00 of Jan. 1 of the next year.
    first_history_yr = first_season + 1

    # Same deal for the last history timestep, but we have to read an extra year in that case,
    # because in some places the last growing season won't complete until the year after it was
    # planted.
    last_history_yr = last_season + 2

    # last_history_yr + 1 because range() will iterate up to but not including the second value.
    history_yr_range = range(first_history_yr, last_history_yr + 1)

    return history_yr_range


def _get_time_slice_lists(first_season, last_season):
    """
    Given the requested first and last seasons, get the list of time slices that the script should
    look for. The assumption here, as in _get_file_lists() and as instructed in the docs, is
    that the user is saving instantaneous files.
    """

    # Input checks
    if not all(isinstance(i, int) for i in [first_season, last_season]):
        raise TypeError("_get_time_slice_list() arguments must be integers")
    if first_season > last_season:
        raise ValueError(f"first_season ({first_season}) > last_season ({last_season})")

    slice_lists_list = [None, None]
    for i, h in enumerate([1, 2]):
        slice_list = []
        for history_yr in _get_history_yr_range(first_season, last_season):
            if h == 1:
                # Annual timesteps
                slice_start = f"{history_yr}-01-01"
                slice_stop = f"{history_yr}-01-01"
            elif h == 2:
                # Daily timesteps of instantaneous variables will go from Jan. 2 through Jan. 1
                # because they will get the time at the end of each timestep.
                slice_start = f"{history_yr}-01-02"
                slice_stop = f"{history_yr + 1}-01-01"
            else:
                raise NotImplementedError(f"What frequency are h{h}i files saved at?")
            slice_list.append(slice(slice_start, slice_stop))

        # We should be reading one more than the total number of years in
        # [first_season, last_season].
        assert len(slice_list) == last_season - first_season + 2

        # Save
        slice_lists_list[i] = slice_list

    return tuple(slice_lists_list)


def _get_file_lists(input_dir, time_slice_lists_list, logger):
    """
    For each time slice in a list, find the file(s) that need to be read to get all history
    timesteps in the slice. Returns both h1i and h2i file lists.
    """
    output_file_lists_list = [None, None]
    for i, h in enumerate([1, 2]):
        time_slice_list = time_slice_lists_list[i]
        all_h_files = gddfn.find_inst_hist_files(input_dir, h=h, logger=logger)
        h_file_lists = []
        for time_slice in time_slice_list:
            try:
                h_file_lists.append(
                    get_files_in_time_slice(all_h_files, time_slice, logger=logger, quiet=True)
                )
            except FileNotFoundError as e:
                raise FileNotFoundError(f"No h{h} timesteps found in {time_slice}") from e
        output_file_lists_list[i] = h_file_lists
    h1_file_lists, h2_file_lists = tuple(output_file_lists_list)
    return h1_file_lists, h2_file_lists


def main(
    *,
    input_dir=None,
    first_season=None,
    last_season=None,
    sdates_file=None,
    hdates_file=None,
    output_dir=None,
    save_figs=True,
    only_make_figs=False,
    run1_name=None,
    run2_name=None,
    land_use_file=None,
    first_land_use_year=None,
    last_land_use_year=None,
    max_season_length_from_hdates_file=False,
    skip_crops=None,
    logger=None,
    no_pickle=None,
    paramfile=None,
    max_season_length_cushion=None,
):  # pylint: disable=missing-function-docstring,too-many-statements
    # Directories to save output files and figures
    if not output_dir:
        if only_make_figs:
            output_dir = input_dir
        else:
            output_dir = os.path.join(input_dir, "generate_gdds")
            if not max_season_length_from_hdates_file:
                output_dir += ".mxmat"
            output_dir += "." + dt.datetime.now().strftime("%Y-%m-%d-%H%M%S")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    outdir_figs = os.path.join(output_dir, "figs")

    # Set up log file and function, if needed
    if logger is None:
        logging.basicConfig(
            level=logging.INFO,
            format="",
            filename=os.path.join(output_dir, "generate_gdds.log"),
            filemode="a",
        )
        logger = logging.getLogger("")

    # Disable plotting if any plotting module is unavailable
    if save_figs:
        try:
            # pylint: disable=import-outside-toplevel,unused-import,import-error
            import cartopy
            import matplotlib
        except ModuleNotFoundError as exc:
            if only_make_figs:
                raise RuntimeError(
                    "only_make_figs True but not all plotting modules are available"
                ) from exc
            log(logger, "Not all plotting modules are available; disabling save_figs")
            save_figs = False

    # Print some info
    log(logger, f"Saving to {output_dir}")

    # Parse list of crops to skip
    if "," in skip_crops:
        skip_crops = skip_crops.split(",")
    else:
        skip_crops = skip_crops.split(" ")

    ##########################
    ### Import and process ###
    ##########################

    if not only_make_figs:
        # This script uses pickle to save work in progress. In case of interruption, when the script
        # is resumed, it will look for a pickle file. It will resume from the year after
        # pickle_year, which is the last processed year in the pickle file.
        pickle_file = os.path.join(output_dir, f"{first_season}-{last_season}.pickle")
        h2_ds_file = os.path.join(output_dir, f"{first_season}-{last_season}.h2_ds.nc")
        if os.path.exists(pickle_file) and not no_pickle:
            with open(pickle_file, "rb") as file:
                (
                    first_season,
                    last_season,
                    pickle_season,
                    gddaccum_yp_list,
                    gddharv_yp_list,
                    skip_patches_for_isel_nan_lastyear,
                    lastyear_active_patch_indices_list,
                    save_figs,
                    incl_vegtypes_str,
                    incl_patches1d_itype_veg,
                    mxsowings,
                    skip_crops,
                ) = pickle.load(file)
            log(logger, f"Will resume import at season {pickle_season+1}")
            h2_ds = None
        else:
            skip_patches_for_isel_nan_lastyear = np.ndarray([])
            pickle_season = -np.inf
            gddaccum_yp_list = []
            gddharv_yp_list = []
            incl_vegtypes_str = None
            lastyear_active_patch_indices_list = None
        sdates_rx = sdates_file
        hdates_rx = hdates_file

        mxmats = _get_max_growing_season_lengths(
            max_season_length_from_hdates_file, paramfile, max_season_length_cushion
        )

        # Get lists of history timesteps and files to read
        h1_time_slices, h2_time_slices = _get_time_slice_lists(first_season, last_season)
        h1_file_lists, h2_file_lists = _get_file_lists(
            input_dir, (h1_time_slices, h2_time_slices), logger
        )

        for y, history_yr in enumerate(_get_history_yr_range(first_season, last_season)):
            # If resuming from a pickled file, we continue until we reach a year that hasn't yet
            # been processed.
            if history_yr <= pickle_season:
                continue
            log(logger, f"History year {history_yr}...")

            # Get time slice and files to read for this year
            h1_time_slice = h1_time_slices[y]  # pylint: disable=unsubscriptable-object
            h2_time_slice = h2_time_slices[y]  # pylint: disable=unsubscriptable-object
            h1_file_list = h1_file_lists[y]  # pylint: disable=unsubscriptable-object
            h2_file_list = h2_file_lists[y]  # pylint: disable=unsubscriptable-object

            (
                h2_ds,
                sdates_rx,
                hdates_rx,
                gddaccum_yp_list,
                gddharv_yp_list,
                skip_patches_for_isel_nan_lastyear,
                lastyear_active_patch_indices_list,
                incl_vegtypes_str,
                incl_patches1d_itype_veg,
                mxsowings,
            ) = gddfn.import_and_process_1yr(
                first_season,
                last_season,
                y,
                sdates_rx,
                hdates_rx,
                gddaccum_yp_list,
                gddharv_yp_list,
                skip_patches_for_isel_nan_lastyear,
                lastyear_active_patch_indices_list,
                incl_vegtypes_str,
                h2_ds_file,
                mxmats,
                cc.get_gs_len_da,
                skip_crops,
                outdir_figs,
                logger,
                h1_file_list,
                h2_file_list,
                h1_time_slice,
                h2_time_slice,
            )

            log(logger, f"   Saving pickle file ({pickle_file})...")
            with open(pickle_file, "wb") as file:
                pickle.dump(
                    [
                        first_season,
                        last_season,
                        history_yr,
                        gddaccum_yp_list,
                        gddharv_yp_list,
                        skip_patches_for_isel_nan_lastyear,
                        lastyear_active_patch_indices_list,
                        save_figs,
                        incl_vegtypes_str,
                        incl_patches1d_itype_veg,
                        mxsowings,
                        skip_crops,
                    ],
                    file,
                    protocol=-1,
                )

        if isinstance(incl_vegtypes_str, list):
            incl_vegtypes_str = np.array(incl_vegtypes_str)
        plot_vegtypes_str = incl_vegtypes_str[
            [i for i, c in enumerate(gddaccum_yp_list) if not isinstance(c, type(None))]
        ]

        log(logger, "Done")

        if not h2_ds:
            log(logger, f"Opening h2_ds: {h2_ds_file}")
            h2_ds = xr.open_dataset(h2_ds_file)

    ######################################################
    ### Get and grid mean GDDs in GGCMI growing season ###
    ######################################################

    if not only_make_figs:
        longname_prefix = "GDD harvest target for "

        # Could skip this by saving sdates_rx['time_bounds']
        sdates_rx = gddfn.import_rx_dates(
            "s", sdates_rx, incl_patches1d_itype_veg, mxsowings, logger
        )

        log(logger, "Getting and gridding mean GDDs...")
        gdd_maps_ds = gddfn.yp_list_to_ds(
            gddaccum_yp_list, h2_ds, incl_vegtypes_str, sdates_rx, longname_prefix, logger
        )
        gddharv_maps_ds = gddfn.yp_list_to_ds(
            gddharv_yp_list, h2_ds, incl_vegtypes_str, sdates_rx, longname_prefix, logger
        )

        # Fill NAs with dummy values
        dummy_fill = -1
        gdd_maps_ds = gdd_maps_ds.fillna(dummy_fill)
        log(logger, "Done getting and gridding means.")

        # Add dummy variables for crops not actually simulated
        log(logger, "Adding dummy variables...")
        # Unnecessary?
        template_ds = xr.open_dataset(sdates_file, decode_times=True)
        all_vars = [v.replace("sdate", "gdd") for v in template_ds if "sdate" in v]
        all_longnames = [
            template_ds[v].attrs["long_name"].replace("Planting day ", longname_prefix) + " (dummy)"
            for v in template_ds
            if "sdate" in v
        ]
        dummy_vars = []
        dummy_longnames = []
        for var_index, this_var in enumerate(all_vars):
            if this_var not in gdd_maps_ds:
                dummy_vars.append(this_var)
                dummy_longnames.append(all_longnames[var_index])

        def make_dummy(this_crop_gridded, addend):
            dummy_gridded = this_crop_gridded
            dummy_gridded.values = dummy_gridded.values * 0 + addend
            return dummy_gridded

        for var in gdd_maps_ds:
            this_crop_gridded = gdd_maps_ds[var].copy()
            break
        dummy_gridded = make_dummy(this_crop_gridded, -1)

        for var_index, this_var in enumerate(dummy_vars):
            if this_var in gdd_maps_ds:
                error(logger, f"{this_var} is already in gdd_maps_ds. Why overwrite it with dummy?")
            dummy_gridded.name = this_var
            dummy_gridded.attrs["long_name"] = dummy_longnames[var_index]
            gdd_maps_ds[this_var] = dummy_gridded

        # Add lon/lat attributes
        def add_lonlat_attrs(this_ds):
            this_ds.lon.attrs = {"long_name": "coordinate_longitude", "units": "degrees_east"}
            this_ds.lat.attrs = {"long_name": "coordinate_latitude", "units": "degrees_north"}
            return this_ds

        gdd_maps_ds = add_lonlat_attrs(gdd_maps_ds)
        gddharv_maps_ds = add_lonlat_attrs(gddharv_maps_ds)

        log(logger, "Done.")

    ######################
    ### Save to netCDF ###
    ######################

    if not only_make_figs:
        log(logger, "Saving...")

        # Get output file path
        datestr = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        outfile = os.path.join(output_dir, "gdds_" + datestr + ".nc")

        def save_gdds(sdates_file, hdates_file, outfile, gdd_maps_ds, sdates_rx):
            # Set up output file from template (i.e., prescribed sowing dates).
            template_ds = xr.open_dataset(sdates_file, decode_times=True)
            for var in template_ds:
                if "sdate" in var:
                    template_ds = template_ds.drop(var)
            template_ds.to_netcdf(path=outfile, format="NETCDF4_CLASSIC")
            template_ds.close()

            # Add global attributes
            comment = (
                "Derived from CLM run plus crop calendar input files "
                + f"{os.path.basename(sdates_file) and {os.path.basename(hdates_file)}}."
            )
            gdd_maps_ds.attrs = {
                "author": "Sam Rabin (sam.rabin@gmail.com)",
                "comment": comment,
                "created": dt.datetime.now().astimezone().isoformat(),
            }

            # Add time_bounds
            if "time_bounds" in sdates_rx:
                gdd_maps_ds["time_bounds"] = sdates_rx.time_bounds

            # Save cultivar GDDs
            gdd_maps_ds.to_netcdf(outfile, mode="w", format="NETCDF4_CLASSIC")

        save_gdds(sdates_file, hdates_file, outfile, gdd_maps_ds, sdates_rx)

        log(logger, "Done saving.")

    ########################################
    ### Save things needed for mapmaking ###
    ########################################

    def add_attrs_to_map_ds(
        map_ds, incl_vegtypes_str, dummy_fill, outdir_figs, first_season, last_season
    ):
        return map_ds.assign_attrs(
            {
                "incl_vegtypes_str": incl_vegtypes_str,
                "dummy_fill": dummy_fill,
                "outdir_figs": outdir_figs,
                "y1": first_season,
                "yN": last_season,
            }
        )

    if not only_make_figs:
        if not os.path.exists(outdir_figs):
            os.makedirs(outdir_figs)

        gdd_maps_ds = add_attrs_to_map_ds(
            gdd_maps_ds, plot_vegtypes_str, dummy_fill, outdir_figs, first_season, last_season
        )
        gddharv_maps_ds = add_attrs_to_map_ds(
            gddharv_maps_ds, plot_vegtypes_str, dummy_fill, outdir_figs, first_season, last_season
        )

        gdd_maps_ds.to_netcdf(os.path.join(outdir_figs, "gdd_maps.nc"))
        gddharv_maps_ds.to_netcdf(os.path.join(outdir_figs, "gddharv_maps.nc"))

    #################################################
    ### Save before/after map and boxplot figures ###
    #################################################

    if save_figs:
        if only_make_figs:
            gdd_maps_ds = xr.open_dataset(os.path.join(input_dir, "figs", "gdd_maps.nc"))
            gddharv_maps_ds = xr.open_dataset(os.path.join(input_dir, "figs", "gddharv_maps.nc"))
        gddfn.make_figures(
            first_land_use_year,
            last_land_use_year,
            land_use_file,
            run1_name,
            run2_name,
            logger,
            gdd_maps_ds=gdd_maps_ds,
            gddharv_maps_ds=gddharv_maps_ds,
            outdir_figs=outdir_figs,
        )


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        description=(
            "A script to generate maturity requirements for CLM crops in units of growing degree-"
            "days (GDDs)."
        )
    )

    # Required but mutually exclusive
    max_growing_season_length_group = parser.add_mutually_exclusive_group(required=True)
    max_growing_season_length_group.add_argument(
        "--paramfile",
        help=(
            "Path to parameter file with maximum growing season lengths (mxmat)."
            " Mutually exclusive with --max-season-length-from-hdates-file."
        ),
    )
    max_growing_season_length_group.add_argument(
        "--max-season-length-from-hdates-file",
        help=(
            "Rather than limiting growing season length based on mxmat values from a CLM parameter"
            " file, use the season lengths from --hdates-file. Not recommended unless you use the"
            "results of this script in a run with sufficiently long mxmat values!"
            " Mutually exclusive with --paramfile."
        ),
        action="store_true",
        default=False,
    )

    # Required
    parser.add_argument(
        "-i",
        "--input-dir",
        help=(
            "Directory where run outputs can be found (and where outputs will go). If "
            + "--only-make-figs, this is the directory with the preprocessed files (e.g., *.pickle "
            + "file)."
        ),
        required=True,
    )
    parser.add_argument(
        "-1",
        "--first-season",
        help="First growing season to include in calculation of mean",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-n",
        "-N",
        "--last-season",
        help="Last growing season to include in calculation of mean",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Output directory. Default is auto-generated subdir of -i/--input-dir.",
    )
    parser.add_argument(
        "-sd", "--sdates-file", help="File of prescribed sowing dates", required=True
    )
    parser.add_argument(
        "-hd", "--hdates-file", help="File of prescribed harvest dates", required=True
    )

    # Optional
    figsgroup = parser.add_mutually_exclusive_group()
    figsgroup.add_argument(
        "--dont-save-figs", help="Do not save figures", action="store_true", default=False
    )
    figsgroup.add_argument(
        "--only-make-figs",
        help="Use preprocessed files to make figures only",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--run1-name", help="Name of original values to show in figures", default="Old"
    )
    parser.add_argument("--run2-name", help="Name of new values to show in figures", default="New")
    parser.add_argument(
        "-lu",
        "--land-use-file",
        help="Path to CLM land use timeseries file, for masking figures",
        default=None,
    )
    parser.add_argument(
        "--first-land-use-year",
        help="First year in land use file to use for masking. Default --first-season.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--last-land-use-year",
        help="Last year in land use file to use for masking. Default --last-season.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--skip-crops",
        help="Skip processing of these crops. Comma- or space-separated list.",
        type=str,
        default="",
    )
    parser.add_argument(
        "--no-pickle",
        help="Don't read from existing pickle file; instead, overwrite. For troubleshooting.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--max-season-length-cushion",
        help=(
            "How much to reduce the maximum growing season length (mxmat) for each crop in the"
            " parameter file. This might be useful for helping avoid high rates of immature"
            " harvests for gridcells where the observed harvest date is longer than mxmat."
            " Incompatible with --max-season-length-from-hdates-file."
        ),
        default=0,
        type=int,
    )

    # Get arguments
    args_parsed = parser.parse_args(argv)
    for k, v in sorted(vars(args_parsed).items()):
        print(f"{k}: {v}")

    # Check arguments
    if args_parsed.max_season_length_from_hdates_file and args_parsed.max_season_length_cushion:
        raise argparse.ArgumentError(
            None,
            "--max-season-length-from-hdates-file is incompatible with --max-season-length-cushion"
            " ≠ 0.",
        )

    return args_parsed


if __name__ == "__main__":
    ###############################
    ### Process input arguments ###
    ###############################
    args = _parse_args(sys.argv[1:])

    # Call main()
    main(
        input_dir=args.input_dir,
        first_season=args.first_season,
        last_season=args.last_season,
        sdates_file=args.sdates_file,
        hdates_file=args.hdates_file,
        output_dir=args.output_dir,
        save_figs=not args.dont_save_figs,
        only_make_figs=args.only_make_figs,
        run1_name=args.run1_name,
        run2_name=args.run2_name,
        land_use_file=args.land_use_file,
        first_land_use_year=args.first_land_use_year,
        last_land_use_year=args.last_land_use_year,
        max_season_length_from_hdates_file=args.max_season_length_from_hdates_file,
        skip_crops=args.skip_crops,
        no_pickle=args.no_pickle,
        paramfile=args.paramfile,
        max_season_length_cushion=args.max_season_length_cushion,
    )
