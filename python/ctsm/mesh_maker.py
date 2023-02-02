#!/usr/bin/env python3
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This script creates ESMF unstructured GRID (mesh file) from a netcdf
file with valid lats and lons. Provided lats and lons can be 1D or 2D.

For example for running WRF-CTSM cases, the user can create a mesh
file for their domain :
    ./mesh_maker --input wrfinput_d01 --output my_region \
        --lat XLAT --lon XLONG --verbose

"""
import os
import sys
import logging
import argparse
import textwrap
from datetime import datetime
import xarray as xr
import numpy as np

from ctsm.site_and_regional.mesh_type import MeshType
from ctsm.utils import abort
from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)


def get_parser():
    """
    Get the parser object for mesh_maker script.

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
        help="Netcdf input file for creating ESMF mesh.",
        action="store",
        dest="input",
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Name of the ESMF mesh created.",
        action="store",
        dest="output",
        required=False,
    )

    parser.add_argument(
        "--outdir",
        help="Output directory (only if name of output mesh is not defined)",
        action="store",
        dest="out_dir",
        type=str,
    )

    parser.add_argument(
        "--lat",
        help="Name of latitude varibale on netcdf input file.",
        action="store",
        dest="lat_name",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--lon",
        help="Name of latitude varibale on netcdf input file.",
        action="store",
        dest="lon_name",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--mask",
        help="Name of mask varibale on netcdf input file."
        + " If none given, create a fake mask with values of 1.",
        action="store",
        dest="mask_name",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--no-plot",
        help="Do not do the plots of the mesh",
        action="store_true",
        dest="noplot",
        required=False,
    )
    parser.add_argument(
        "--area",
        help="Name of area variable on netcdf input file."
        + " If none given, ESMF calculates element areas automatically. ",
        action="store",
        dest="area_name",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--overwrite",
        help="If meshfile exists, overwrite the meshfile.",
        action="store_true",
        dest="overwrite",
        required=False,
    )

    add_logging_args(parser)
    return parser


def check_args(args):
    """Check the arguments"""
    if args.output and args.out_dir:
        logging.error(" Both --outdir and --output cannot be provided at the same time.")
        err_msg = textwrap.dedent(
            """
                \n ------------------------------------
                \n You have provided both --outdir and --output.
                \n Please provide only one of these options to proceed:
                \n --outdir : directory to save mesh file. mesh file name automatically created.
                \n --output : Absolute or relative path of the ESMF mesh file created.\n
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
        args.output = os.path.join(
            args.out_dir,
            os.path.splitext(args.input)[0]
            + "_ESMF_UNSTRUCTURED_MESH"
            + "_c"
            + today_string
            + ".nc",
        )

    # -- exit if mesh_out exists and --overwrite is not specified.
    if os.path.exists(args.output):
        if args.overwrite:
            os.remove(args.output)
        else:
            err_msg = (
                "output meshfile exists, please choose --overwrite to overwrite the mesh file."
            )
            abort(err_msg)


def main():
    """Main function to create a mesh file from another file"""

    setup_logging_pre_config()
    parser = get_parser()
    args = parser.parse_args()

    # --------------------------------- #
    # process logging args (i.e. debug and verbose)
    process_logging_args(args)

    check_args(args)

    nc_file = args.input
    lat_name = args.lat_name
    lon_name = args.lon_name
    mesh_out = args.output
    mask_name = args.mask_name
    area_name = args.area_name

    if os.path.isfile(nc_file):
        ds = xr.open_dataset(nc_file, mask_and_scale=False, decode_times=False).transpose()
    else:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n Input file not found. Please make sure to provide the
                \n full path of Netcdf input file for making the mesh.
                \n ------------------------------------
                """
        )
        abort(err_msg)

    if lat_name not in ds.coords and lat_name not in ds.variables:
        logging.error("Input file does not have variable named %s", lat_name)
        sys.exit()

    else:
        logging.debug(
            "- %s exist in the provided netcdf file with dimension of %s.",
            lat_name,
            len(ds[lat_name].dims).__str__(),
        )

    if lon_name not in ds.coords and lon_name not in ds.variables:
        logging.error("Input file does not have variable named %s", lon_name)
        sys.exit()
    else:
        logging.debug(
            "- %s exist in the provided netcdf file with dimension of %s.",
            lon_name,
            len(ds[lon_name].dims).__str__(),
        )

    lats = ds[lat_name].astype(np.float32)
    lons = ds[lon_name].astype(np.float32)

    if (len(lats.dims) > 2) or (len(lons.dims) > 2):
        time_dims = [dim for dim in lats.dims if "time" in dim.lower()]
        if time_dims:
            logging.debug("- time dimension found on lat %s", str(time_dims))
            logging.debug("- removing time dimensions from lats and lons. ")
            lats = lats[:, :, 0]
            lons = lons[:, :, 0]
        else:
            err_msg = (
                "latitude or longitude has more than 2 dimensions and "
                + "the third dimension cannot be detected as time."
            )
            abort(err_msg)

    logging.info("Creating mesh file from : %s", nc_file)
    logging.info("Writing mesh file to    : %s", mesh_out)

    if mask_name is not None:
        mask = ds[mask_name].values()
    else:
        mask = None

    if area_name is not None:
        area = ds[area_name].values()
    else:
        area = None

    this_mesh = MeshType(lats, lons, mask=mask)
    this_mesh.calculate_corners()
    this_mesh.create_esmf(mesh_out, area=area)


if __name__ == "__main__":
    main()
