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
        help="Name of latitude variable on netcdf input file.",
        action="store",
        dest="lat_name",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--lon",
        help="Name of latitude variable on netcdf input file.",
        action="store",
        dest="lon_name",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--mask",
        help="Name of mask variable on netcdf input file."
        + " If none given, create a fake mask with values of 1.",
        action="store",
        dest="mask_name",
        type=str,
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


def process_and_check_args(args):
    """Process and check the arguments"""
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

    if not os.path.isfile(args.input):
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n Input file not found. Please make sure to provide the
                \n full path of Netcdf input file for making the mesh.
                \n ------------------------------------
                """
        )
        abort(err_msg)

    return args


def check_input_file(args, ds):
    """Check that the input file has the variables expected"""
    if args.lat_name not in ds.coords and args.lat_name not in ds.variables:
        err_msg = "Input file does not have variable named " + args.lat_name
        abort(err_msg)

    else:
        logging.debug(
            "- %s exist in the provided netcdf file with dimension of %s.",
            args.lat_name,
            str(len(ds[args.lat_name].dims)),
        )

    if args.lon_name not in ds.coords and args.lon_name not in ds.variables:
        err_msg = "Input file does not have variable named " + args.lon_name
        abort(err_msg)
    else:
        logging.debug(
            "- %s exist in the provided netcdf file with dimension of %s.",
            args.lon_name,
            str(len(ds[args.lon_name].dims)),
        )
    if args.mask_name is not None:
        if args.mask_name not in ds.variables:
            err_msg = "Input file does not have mask variable named " + args.mask_name
            abort(err_msg)

    if args.area_name is not None:
        if args.area_name not in ds.variables:
            err_msg = "Input file does not have area variable named " + args.area_name
            abort(err_msg)
        if "units" not in ds[args.area_name].attrs:
            err_msg = "Units attribute is NOT on the area variable"
            abort(err_msg)
        areaunits = ["radians^2", "radian2", "radians2", "radian^2"]
        if not any(name in ds[args.area_name].attrs["units"] for name in areaunits):
            err_msg = "Area does NOT have the correct units of radians^2 but has " + str(
                ds[args.area_name].attrs["units"]
            )
            abort(err_msg)
    if len(ds[args.lon_name]) == 1:
        err_msg = "No need to create a mesh file for a single point grid."
        abort(err_msg)


def main():
    """Main function to create a mesh file from another file"""

    setup_logging_pre_config()
    parser = get_parser()
    args = parser.parse_args()

    # --------------------------------- #
    # process logging args (i.e. debug and verbose)
    process_logging_args(args)

    args = process_and_check_args(args)

    nc_file = args.input
    mesh_out = args.output
    mask_name = args.mask_name
    area_name = args.area_name

    ds = xr.open_dataset(nc_file, mask_and_scale=False, decode_times=False).transpose()

    check_input_file(args, ds)

    lats = ds[args.lat_name].astype(np.float32)
    lons = ds[args.lon_name].astype(np.float32)

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
        mask = ds[mask_name].astype(np.float32)
        if mask.max() > 1.0 or mask.min() < 0.0:
            abort("Mask variable is not within 0 to 1")
    else:
        mask = None

    if area_name is not None:
        area = np.array(ds[area_name].astype(np.float32))
    else:
        area = None

    this_mesh = MeshType(lats, lons, mask=mask, area=area)
    this_mesh.calculate_corners()
    this_mesh.calculate_nodes()
    this_mesh.create_esmf(mesh_out)


def read_main():
    """Main function to read a mesh file and output it again"""

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
    logging.info("Writing mesh file to    : %s", mesh_out)

    this_mesh = MeshType(lats, lons)
    this_mesh.read_file(ds)
    this_mesh.create_esmf(mesh_out)


if __name__ == "__main__":
    main()
