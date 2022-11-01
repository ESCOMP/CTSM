#!/usr/bin/env python3

from mesh_type import MeshType
import xarray as xr
import os
import argparse
import sys
import textwrap

def get_parser():
    """
    Get the parser object for mesh_maker.py script.

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
        "--lat",
        help="Name of latitude varibale on netcdf input file. If none given, looks to find variables that include 'lat'.",
        action="store",
        dest="lat_name",
        type = str,
        required=True,
    )
    parser.add_argument(
        "--lon",
        help="Name of latitude varibale on netcdf input file. If none given, looks to find variables that include 'lon'.",
        action="store",
        dest="lon_name",
        type = str,
        required=True,
    )
    parser.add_argument(
        "--mask",
        help="Name of mask varibale on netcdf input file. If none given, create a fake mask with values of 1.",
        action="store",
        dest="mask_name",
        type = str,
        required=False,
    )
    parser.add_argument(
        "--area",
        help="Name of area variable on netcdf input file. If none given, ESMF calculates element areas automatically. ",
        action="store",
        dest="area_name",
        type = str,
        required=False,
    )
    parser.add_argument(
        "--overwrite",
        help="If meshfile exists, overwrite the meshfile.",
        action="store_true",
        dest="overwrite",
        required=False,
    )
    return parser


def main ():
    oned = True
    if oned:
        ifile = "/glade/scratch/negins/this_region_4x5_new_2/surfdata_4x5_hist_78pfts_CMIP6_simyr1850_275.0-330.0_-40-15_c220705.nc"
    else:
        ifile = "/glade/scratch/negins/wrf-ctsm_production/WRF/test/em_real_sim1_noahmp/wrfinput_d01"

    if oned:
        lat_name = "lsmlat"
        lon_name = "lsmlon"
    else:
        lat_name = "XLAT"
        lon_name = "XLONG"

    parser = get_parser()
    args = parser.parse_args()

    nc_file = args.input
    lat_name = args.lat_name
    lon_name = args.lon_name
    mesh_out = args.output
    overwrite = args.overwrite
    mask_name = args.mask_name
    area_name = args.area_name

    if os.path.isfile(nc_file):
        ds = xr.open_dataset(nc_file, mask_and_scale=False, decode_times=False).transpose()
    else:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n Input file not found. Please make sure to provide the full path of Netcdf input file for making mesh.
                \n ------------------------------------
                """
        )
        raise parser.error(err_msg)

    if lat_name not in ds.coords and lat_name not in ds.variables :
        print('Input file does not have variable named {}.'.format(lat_name))
    else:
        print (lat_name, "exist in the provided netcdf file with dimension of ", len(ds[lat_name].dims))

    if lon_name not in ds.coords and lat_name not in ds.variables :
        print('Input file does not have variable named {}.'.format(lon_name))
    else:
        print (lat_name, "exist in the provided netcdf file with dimension of ", len(ds[lat_name].dims))


    lats = ds[lat_name]
    lons = ds[lon_name]

    if (len(lats.dims)>3) or (len(lons.dims)>3):
        time_dims = [dim for dim in lats.dims if 'time' in dim.lower()]
        if time_dims:
            print ('time dimension found on lat', time_dims)
            lats = lats [:,:,0]
            lats = lons [:,:,0]
        else:
            print ('latitude or longitude has more than 2 dimensions and the third dimension cannot be detected as time.')

    if mesh_out:
        if os.path.exists(mesh_out):
            if overwrite:
                os.remove(mesh_out)
            else:
                print ('output meshfile exists, please choose --overwrite to overwrite the mesh file.')

    if mask_name is not None:
        mask = ds[mask_name].values()

    if area_name is not None:
        area = ds[mask_name].values()

    this_mesh = MeshType(lats, lons, mask=None)
    this_mesh.calculate_corners()
    this_mesh.create_esmf(mesh_out, area=None)
    print ('DONE!')

if __name__ == "__main__":
    main()
