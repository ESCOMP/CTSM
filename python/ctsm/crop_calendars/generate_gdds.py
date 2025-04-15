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
import ctsm.crop_calendars.cropcal_module as cc  # pylint: disable=wrong-import-position
import ctsm.crop_calendars.generate_gdds_functions as gddfn  # pylint: disable=wrong-import-position

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments

# Global constants
PARAMFILE_DIR = "/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/paramdata"
MY_CLM_VER = 51
MY_CLM_SUBVER = "c211112"


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
    unlimited_season_length=False,
    skip_crops=None,
    logger=None,
    no_pickle=None,
):  # pylint: disable=missing-function-docstring,too-many-statements
    # Directories to save output files and figures
    if not output_dir:
        if only_make_figs:
            output_dir = input_dir
        else:
            output_dir = os.path.join(input_dir, "generate_gdds")
            if not unlimited_season_length:
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
            gddfn.log(logger, "Not all plotting modules are available; disabling save_figs")
            save_figs = False

    # Print some info
    gddfn.log(logger, f"Saving to {output_dir}")

    # Parse list of crops to skip
    if "," in skip_crops:
        skip_crops = skip_crops.split(",")
    else:
        skip_crops = skip_crops.split(" ")

    ##########################
    ### Import and process ###
    ##########################

    if not only_make_figs:
        # Keep 1 extra year to avoid incomplete final growing season for crops
        # harvested after Dec. 31.
        yr_1_import_str = f"{first_season+1}-01-01"
        yr_n_import_str = f"{last_season+2}-01-01"

        gddfn.log(
            logger,
            f"Importing netCDF time steps {yr_1_import_str} through {yr_n_import_str} "
            + "(years are +1 because of CTSM output naming)",
        )

        pickle_file = os.path.join(output_dir, f"{first_season}-{last_season}.pickle")
        h2_ds_file = os.path.join(output_dir, f"{first_season}-{last_season}.h2_ds.nc")
        if os.path.exists(pickle_file) and not no_pickle:
            with open(pickle_file, "rb") as file:
                (
                    first_season,
                    last_season,
                    pickle_year,
                    gddaccum_yp_list,
                    gddharv_yp_list,
                    skip_patches_for_isel_nan_lastyear,
                    lastyear_active_patch_indices_list,
                    incorrectly_daily,
                    save_figs,
                    incl_vegtypes_str,
                    incl_patches1d_itype_veg,
                    mxsowings,
                    skip_crops,
                ) = pickle.load(file)
            print(f"Will resume import at {pickle_year+1}")
            h2_ds = None
        else:
            incorrectly_daily = False
            skip_patches_for_isel_nan_lastyear = np.ndarray([])
            pickle_year = -np.inf
            gddaccum_yp_list = []
            gddharv_yp_list = []
            incl_vegtypes_str = None
            lastyear_active_patch_indices_list = None
        sdates_rx = sdates_file
        hdates_rx = hdates_file

        if not unlimited_season_length:
            mxmats = cc.import_max_gs_length(PARAMFILE_DIR, MY_CLM_VER, MY_CLM_SUBVER)
        else:
            mxmats = None

        h1_instantaneous = None
        for yr_index, this_yr in enumerate(np.arange(first_season + 1, last_season + 3)):
            if this_yr <= pickle_year:
                continue

            (
                h2_ds,
                sdates_rx,
                hdates_rx,
                gddaccum_yp_list,
                gddharv_yp_list,
                skip_patches_for_isel_nan_lastyear,
                lastyear_active_patch_indices_list,
                incorrectly_daily,
                incl_vegtypes_str,
                incl_patches1d_itype_veg,
                mxsowings,
                h1_instantaneous,
            ) = gddfn.import_and_process_1yr(
                first_season,
                last_season,
                yr_index,
                this_yr,
                sdates_rx,
                hdates_rx,
                gddaccum_yp_list,
                gddharv_yp_list,
                skip_patches_for_isel_nan_lastyear,
                lastyear_active_patch_indices_list,
                incorrectly_daily,
                input_dir,
                incl_vegtypes_str,
                h2_ds_file,
                mxmats,
                cc.get_gs_len_da,
                skip_crops,
                outdir_figs,
                logger,
                h1_instantaneous,
            )

            gddfn.log(logger, f"   Saving pickle file ({pickle_file})...")
            with open(pickle_file, "wb") as file:
                pickle.dump(
                    [
                        first_season,
                        last_season,
                        this_yr,
                        gddaccum_yp_list,
                        gddharv_yp_list,
                        skip_patches_for_isel_nan_lastyear,
                        lastyear_active_patch_indices_list,
                        incorrectly_daily,
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

        gddfn.log(logger, "Done")

        if not h2_ds:
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

        gddfn.log(logger, "Getting and gridding mean GDDs...")
        gdd_maps_ds = gddfn.yp_list_to_ds(
            gddaccum_yp_list, h2_ds, incl_vegtypes_str, sdates_rx, longname_prefix, logger
        )
        gddharv_maps_ds = gddfn.yp_list_to_ds(
            gddharv_yp_list, h2_ds, incl_vegtypes_str, sdates_rx, longname_prefix, logger
        )

        # Fill NAs with dummy values
        dummy_fill = -1
        gdd_maps_ds = gdd_maps_ds.fillna(dummy_fill)
        gddfn.log(logger, "Done getting and gridding means.")

        # Add dummy variables for crops not actually simulated
        gddfn.log(logger, "Adding dummy variables...")
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
                gddfn.error(
                    logger, f"{this_var} is already in gdd_maps_ds. Why overwrite it with dummy?"
                )
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

        gddfn.log(logger, "Done.")

    ######################
    ### Save to netCDF ###
    ######################

    if not only_make_figs:
        gddfn.log(logger, "Saving...")

        # Get output file path
        datestr = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        outfile = os.path.join(output_dir, "gdds_" + datestr + ".nc")

        def save_gdds(sdates_file, hdates_file, outfile, gdd_maps_ds, sdates_rx):
            # Set up output file from template (i.e., prescribed sowing dates).
            template_ds = xr.open_dataset(sdates_file, decode_times=True)
            for var in template_ds:
                if "sdate" in var:
                    template_ds = template_ds.drop(var)
            template_ds.to_netcdf(path=outfile, format="NETCDF3_CLASSIC")
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
            gdd_maps_ds.to_netcdf(outfile, mode="w", format="NETCDF3_CLASSIC")

        save_gdds(sdates_file, hdates_file, outfile, gdd_maps_ds, sdates_rx)

        gddfn.log(logger, "Done saving.")

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


if __name__ == "__main__":
    ###############################
    ### Process input arguments ###
    ###############################
    parser = argparse.ArgumentParser(description="ADD DESCRIPTION HERE")

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
        "--unlimited-season-length",
        help="Limit mean growing season length based on CLM CFT parameter mxmat.",
        action="store_true",
        default=False,
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

    # Get arguments
    args = parser.parse_args(sys.argv[1:])
    for k, v in sorted(vars(args).items()):
        print(f"{k}: {v}")

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
        unlimited_season_length=args.unlimited_season_length,
        skip_crops=args.skip_crops,
        no_pickle=args.no_pickle,
    )
