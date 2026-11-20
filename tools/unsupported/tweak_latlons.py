"""
'Tweak' the latitude and longitude coordinates to avoid ambiguous nearest neighbors
"""
import os
import sys
import contextlib
import argparse
import numpy as np
import xarray as xr
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

# -- add python/ctsm  to path (needed if we want to run this stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python")
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.mesh_maker import main as mesh_maker

COORD_LIST = ["lat", "lon"]
COORD_DATATYPE = np.float64

def get_tweak(ds_in, coord_str, init_tweak):
    """
    Get the tweak that will be applied to all datasets' lat/lon coordinates
    """
    da = ds_in[coord_str]
    coord2_orig = da.values.astype(COORD_DATATYPE)
    coord2 = coord2_orig
    tweak = init_tweak
    coord2 += tweak

    # This is necessary if precision is lower than float64
    max_tweak = 1e-2
    while np.any(coord2 == da.values):
        tweak *= 10
        if tweak > max_tweak:
            raise RuntimeError(f"Tweaking by +{max_tweak} failed to 'take'")
        coord2 = coord2_orig
        coord2 += tweak
    return tweak

def apply_tweak(ds_in, coord_str, tweak):
    # Apply tweak
    da = ds_in[coord_str]
    coord2 = da.values.astype(COORD_DATATYPE)
    coord2 += tweak
    if np.any(coord2 == da.values):
        raise RuntimeError('Tweak didn''t "take"')
    coord_tweak = np.full_like(coord2, tweak)

    # Ensure that no value is above maximum in input data. This is needed for mesh_maker to work.
    max_coord = np.max(da.values)
    where_toohigh = np.where(coord2 > max_coord)
    Ntoohigh = len(where_toohigh[0])
    if Ntoohigh != 1:
        raise RuntimeError(
            f"Expected 1 coordinate value too high; got {Ntoohigh}"
        )
    coord2[where_toohigh] = max_coord
    coord_tweak[where_toohigh] = max_coord

    # Convert to DataArray
    new_coords_dict = {coord_str: coord2}
    da2 = xr.DataArray(
        data=coord2,
        coords=new_coords_dict,
        dims=da.dims,
        attrs=da.attrs,
    )

    # Replace coordinate in dataset
    ds_in[coord_str] = da2

    # Add a variable with the amount of the tweak
    tweak_attrs = {}
    if "standard_name" in da.attrs:
        coord_name = da.attrs["standard_name"]
    elif "long_name" in da.attrs:
        coord_name = da.attrs["long_name"].replace("coordinate_", "")
    else:
        coord_name = coord_str
    tweak_attrs["standard_name"] = coord_name + "_tweak"
    tweak_attrs[
        "long_name"
    ] = f"Amount {coord_name} was shifted to avoid ambiguous nearest neighbors"
    if "units" in da.attrs:
        tweak_attrs["units"] = da.attrs["units"]
    da_tweak = xr.DataArray(
        data=coord_tweak,
        coords=new_coords_dict,
        dims=da.dims,
        attrs=tweak_attrs,
    )
    tweak_name = coord_str + "_tweak"
    ds_in[tweak_name] = da_tweak

    return ds_in

def check(ds, f0_base, ds2, f_base, var):
    if not np.array_equal(ds[var].values, ds2[var].values):
        if not np.array_equal(ds[var].shape, ds2[var].shape):
            msg = f"{var} shapes differ b/w {f0_base} ({ds[var].shape}) and {f_base} ({ds2[var].shape})"
            raise RuntimeError(msg)
        max_diff = np.max(np.abs(ds[var].values - ds2[var].values))
        msg = f"{var}s differ between {f0_base} and {f_base}; max = {max_diff}"
        type0 = type(ds[var].values[0])
        type2 = type(ds2[var].values[0])
        if type0 != type2:
            msg += f"\nTypes also differ: {type0} vs. {type2}"
        raise RuntimeError(msg)

@contextlib.contextmanager
def redirect_argv(arglist):
    """
    Preserve actual arg list while giving a new one to mesh_maker
    """
    argv_tmp = sys.argv[:]
    sys.argv = arglist
    yield
    sys.argv = argv_tmp

def main(input_files, mesh_file_in, output_files):
    """
    Apply tweak to all files
    """

    # Set up
    tweak_dict = {}
    for coord in COORD_LIST:
        tweak_dict[coord] = -np.inf
    mesh_file_out = output_files[-1]
    output_files = output_files[:-1]

    # Get tweaks
    for file_in in input_files:
        ds = xr.open_dataset(file_in)
        for coord in COORD_LIST:
            this_tweak = get_tweak(ds, coord, init_tweak=1e-6)
            if this_tweak > tweak_dict[coord]:
                tweak_dict[coord] = this_tweak
    for coord in COORD_LIST:
        print(f"Tweaking {coord} by {tweak_dict[coord]}")
    print(" ")

    # Apply tweaks
    for i, file_in in enumerate(input_files):
        ds = xr.open_dataset(file_in)

        for coord in COORD_LIST:
            ds = apply_tweak(ds, coord, tweak_dict[coord])

        # Set up for save
        file_out = output_files[i]
        with Dataset(file_in, "r") as netcdf_file:
            netcdf_format = netcdf_file.data_model

        # Make output dir, if needed
        output_dir = os.path.dirname(file_out)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Save
        print(f"Saving {file_out}")
        ds.to_netcdf(file_out, format=netcdf_format)
    print("Done")


    # Ensure all files got the same tweaks
    ds = xr.open_dataset(output_files[0])
    f0_base = os.path.basename(output_files[0])
    for file_out in output_files[1:]:
        ds2 = xr.open_dataset(file_out)
        f_base = os.path.basename(file_out)
        for coord in COORD_LIST:
            check(ds, f0_base, ds2, f_base, coord)
            check(ds, f0_base, ds2, f_base, coord + "_tweak")


    # Save new mesh file
    mesh_maker_args = [
        "mesh_maker",
        "--input",
        output_files[0],
        "--output",
        mesh_file_out,
        "--lat",
        "lat",
        "--lon",
        "lon",
        "--overwrite",
    ]
    print(f"Saving {mesh_file_out}...")
    with redirect_argv(mesh_maker_args):
        mesh_maker()

    # Change format, if needed
    with Dataset(mesh_file_in, "r") as netcdf_file:
        netcdf_format_in = netcdf_file.data_model
    with Dataset(mesh_file_out, "r") as netcdf_file:
        netcdf_format_out = netcdf_file.data_model
    if netcdf_format_in != netcdf_format_out:
        mesh_file_out_tmp = mesh_file_out + ".tmp"
        os.rename(mesh_file_out, mesh_file_out_tmp)
        ds = xr.open_dataset(mesh_file_out_tmp)
        ds.to_netcdf(mesh_file_out, format=netcdf_format_in)
        os.remove(mesh_file_out_tmp)

    print("Done")




if __name__ == "__main__":
    ###############################
    ### Process input arguments ###
    ###############################
    parser = argparse.ArgumentParser(
        description="'Tweak' the latitude and longitude coordinates to avoid ambiguous nearest neighbors",
    )

    # Required
    parser.add_argument(
        "-i",
        "--input-files",
        help="Comma-separated stream files whose coordinates need tweaking",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--mesh-file",
        help="Mesh file associated with input files",
        required=True,
    )

    # Optional
    parser.add_argument(
        "--overwrite",
        help="Overwrite any existing output files",
        action="store_true",
        default=False,
    )
    default_output_dir = os.getcwd()
    parser.add_argument(
        "-o",
        "--output-dir",
        help=f"Directory where output files should be saved. Default is current working directory: {default_output_dir}",
        default=default_output_dir,
    )

    # Get arguments
    args = parser.parse_args(sys.argv[1:])

    # Check/process input and output files
    _input_files = args.input_files.split(",")
    _output_files = []
    for file in _input_files + [args.mesh_file]:
        if not os.path.exists(file):
            raise FileNotFoundError(f"File not found: {file}")

        filename, ext = os.path.splitext(os.path.basename(file))
        output_file = os.path.join(
            args.output_dir, filename + ".tweaked_latlons" + ext
        )
        if os.path.exists(output_file) and not args.overwrite:
            raise FileExistsError(
                f"Output file exists but --overwrite not specified: {output_file}"
            )
        _output_files.append(output_file)

    main(_input_files, args.mesh_file, _output_files)
