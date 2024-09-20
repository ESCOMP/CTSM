#!/usr/bin/env python3
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This script plots an ESMF unstructured GRID (mesh file).

"""
import os
import logging
import argparse
import textwrap
from datetime import datetime
import xarray as xr
import numpy as np

from ctsm.site_and_regional.mesh_plot_type import MeshPlotType
from ctsm.utils import abort
from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)


def get_parser():
    """
    Get the parser object for mesh_plotter script.

    Returns:
        parser (ArgumentParser):
            ArgumentParser which includes all the parser information.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        "--input",
        help="Netcdf ESMF Mesh input file to plot.",
        action="store",
        dest="input",
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Prefix for the names of the plot files to create.",
        action="store",
        dest="output",
        required=False,
    )

    parser.add_argument(
        "--outdir",
        help="Output directory for plots",
        action="store",
        dest="out_dir",
        type=str,
    )

    parser.add_argument(
        "--overwrite",
        help="If plots exist, overwrite them.",
        action="store_true",
        dest="overwrite",
        required=False,
    )

    parser.add_argument(
        "--no-center-coords",
        help="Do not include red Xs at center of grid cells.",
        action="store_true",
        required=False,
    )

    default_dpi = 300
    parser.add_argument(
        "--dpi",
        help=f"Dots per square inch in output; default {default_dpi}",
        type=float,
    )

    add_logging_args(parser)
    return parser


def process_and_check_args(args):
    """Process and check the arguments"""
    if args.output and args.out_dir:
        logging.error(" Both --outdir and --output cannot be provided at the same time.")
        err_msg = textwrap.dedent(
            """
                \n ------------------------------------
                \n You have provided both --outdir and --output.
                \n Please provide only one of these options to proceed:
                \n --outdir : directory to save mesh file plots, file names assumed.
                \n --output : Absolute or relative path of the prefix for plot filenames.\n
                """
        )
        abort(err_msg)

    # -- no file name and output path:
    if not args.output and not args.out_dir:
        args.out_dir = os.path.join(os.getcwd(), "meshes")

    if not args.output:
        # -- make output path if does not exist.
        if not os.path.isdir(args.out_dir):
            os.mkdir(args.out_dir)

        today = datetime.today()
        today_string = today.strftime("%y%m%d")
        input_filename = os.path.basename(args.input)
        args.output = os.path.join(
            args.out_dir,
            os.path.splitext(input_filename)[0] + "_c" + today_string,
        )

    if not os.path.isfile(args.input):
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n Input file not found. Please make sure to provide the
                \n full path of Netcdf input mesh file
                \n ------------------------------------
                """
        )
        abort(err_msg)

    return args


def main():
    """Main function to read a mesh file and plot it"""

    setup_logging_pre_config()
    parser = get_parser()
    args = parser.parse_args()

    # --------------------------------- #
    # process logging args (i.e. debug and verbose)
    process_logging_args(args)

    args = process_and_check_args(args)

    nc_file = args.input
    mesh_out = args.output

    ds = xr.open_dataset(nc_file, mask_and_scale=False, decode_times=False).transpose()

    lon0 = np.array([120.0])
    lat0 = np.array([45.0])
    x_dim = "lon"
    y_dim = "lat"
    lons = xr.DataArray(lon0, name="lon", dims=x_dim, coords={x_dim: lon0})
    lats = xr.DataArray(lat0, name="lat", dims=y_dim, coords={y_dim: lat0})

    logging.info("Reading  mesh file from : %s", nc_file)

    this_mesh = MeshPlotType(lats, lons)
    this_mesh.read_file(ds)

    plot_regional = os.path.splitext(mesh_out)[0] + "_regional" + ".png"
    file_exists_msg = "File already exists but --overwrite not given: "
    if os.path.exists(plot_regional) and not args.overwrite:
        raise FileExistsError(file_exists_msg + plot_regional)

    plot_global = os.path.splitext(mesh_out)[0] + "_global" + ".png"
    if os.path.exists(plot_global) and not args.overwrite:
        raise FileExistsError(file_exists_msg + plot_global)

    this_mesh.make_mesh_plot(plot_regional, plot_global, args)


if __name__ == "__main__":
    main()
