"""
Generate stream_fldFileName_gdd20_baseline file from CTSM outputs
"""

import sys
import argparse
import os
import datetime as dt
import numpy as np
import xarray as xr
import cftime

# -- add python/ctsm  to path (needed if we want to run generate_gdd20_baseline stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm.crop_calendars.import_ds import import_ds
import ctsm.crop_calendars.cropcal_utils as utils
from ctsm.crop_calendars.grid_one_variable import grid_one_variable
from ctsm.crop_calendars.cropcal_module import MISSING_RX_GDD_VAL

GRIDDING_VAR_LIST = ["patches1d_ixy", "patches1d_jxy", "lat", "lon"]
STREAM_YEAR = 2000  # The year specified for stream_yearFirst and stream_yearLast in the call of
# shr_strdata_init_from_inline() for sdat_cropcal_gdd20_baseline
MGDCROP_LIST = utils.define_crop_list()

# Functions here were written with too many positional arguments. At some point that should be
# fixed. For now, we'll just disable the warning.
# pylint: disable=too-many-positional-arguments


def _parse_args():
    """
    Set up and parse input arguments
    """
    parser = argparse.ArgumentParser(
        description=(
            "Given a list of CTSM history files, generate stream_fldFileName_gdd20_baseline input"
            + "file from the GDD0, GDD8, and GDD10 variables."
        )
    )

    # Required
    parser.add_argument(
        "-i",
        "--input-files",
        help="Space-separated string of CTSM history files",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Path to which output file should be saved",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-a",
        "--author",
        help=(
            "String to be saved in author attribute of output files."
            + "E.g., 'Author Name (authorname@ucar.edu)'"
        ),
        type=str,
        required=True,
    )

    # Optional
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing output file, if any",
        action="store_true",
    )
    parser.add_argument(
        "-y1",
        "--first-year",
        help=("First calendar year to include"),
        type=int,
        required=False,
    )
    parser.add_argument(
        "-yN",
        "--last-year",
        help=("Last calendar year to include"),
        type=int,
        required=False,
    )
    parser.add_argument(
        "-v",
        "--variable",
        help=("Which type of variable should be processed?"),
        required=False,
        default="GDDBX",
        choices=["GDDBX", "GDDB20"],
    )

    # Get arguments
    args = parser.parse_args(sys.argv[1:])

    # Check arguments
    if os.path.exists(args.output_file) and not args.overwrite:
        raise FileExistsError("Output file exists but --overwrite is not specified")

    # Get and check input files
    args.input_files = args.input_files.split(" ")
    for filename in args.input_files:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Input file not found: {filename}")

    # Process time slice
    # Assumes CESM behavior where data for e.g. 1987 is saved as 1988-01-01.
    # It would be more robust, accounting for upcoming behavior (where timestamp for a year is the
    # middle of that year), to do slice("YEAR1-01-03", "YEARN-01-02"), but that's not compatible
    # with ctsm_pylib as of the version using python 3.7.9. See safer_timeslice() in cropcal_utils.
    if args.first_year is not None:
        date_1 = f"{args.first_year+1}-01-01"
    else:
        date_1 = "0000-01-01"
    if args.last_year is not None:
        date_n = f"{args.last_year+1}-01-01"
    else:
        date_n = "9999-12-31"
    time_slice = slice(date_1, date_n)

    return args, time_slice


def _get_cft_list(crop_list):
    """
    Given a list of strings, return a list of CFT names that contain any of those strings.
    Will include both irrigated and rainfed!

    Args:
        crop_list (list): List of crops to look for.
                          E.g.: ["corn", "cotton"]

    Returns:
        cft_str_list: List of CFTs containing any of the crop names in crop_list.
                      E.g.: ["tropical_corn", "irrigated_tropical_corn",
                             "temperate_corn", "irrigated_temperate_corn",
                             "cotton", "irrigated_cotton"]
    """

    cft_str_list = []
    for crop_str in crop_list:
        cft_str_list += [x for x in MGDCROP_LIST if crop_str in x]
    return cft_str_list


def _get_gddn_for_cft(cft_str, variable):
    """
    Given a CFT name, return the GDDN variable it uses.

    Args:
        cft_str (str): E.g., "irrigated_temperate_corn"

    Returns:
        str or None: Name of variable to use (e.g., "GDD8X"). If crop not yet handled, return None.
    """

    gddn = None
    gddn_str = None

    gdd0_list_str = ["wheat", "cotton", "rice"]
    if cft_str in _get_cft_list(gdd0_list_str):
        gddn = 0

    gdd8_list_str = ["corn", "sugarcane", "miscanthus", "switchgrass"]
    if cft_str in _get_cft_list(gdd8_list_str):
        gddn = 8

    gdd10_list_str = ["soybean"]
    if cft_str in _get_cft_list(gdd10_list_str):
        gddn = 10

    if gddn is not None:
        gddn_str = variable.replace("B", str(gddn))

    return gddn, gddn_str


def _get_output_varname(cft_str):
    cft_int = utils.vegtype_str2int(cft_str)[0]
    return f"gdd20bl_{cft_int}"


def _add_time_axis(da_in):
    """
    Adds a size-1 time axis to a DataArray. Needed because CDEPS streams code requires a time axis,
    even if the data in question is not supposed to vary over time.

    Args:
        da_in (DataArray): xarray DataArray which needs a time axis added

    Returns:
        DataArray: da_in with a new 1-step time axis
    """
    this_date = np.array(cftime.DatetimeNoLeap(STREAM_YEAR, 1, 1, 0, 0, 0, 0, has_year_zero=True))
    this_date = np.expand_dims(this_date, axis=0)
    da_time = xr.DataArray(
        data=this_date,
        dims={"time": this_date},
    )
    da_out = da_in.expand_dims(time=da_time)
    return da_out


def setup_output_dataset(input_files, author, variable, year_args, ds_in):
    """
    Set up output Dataset
    """
    data_var_dict = {}
    for gridding_var in GRIDDING_VAR_LIST:
        data_var_dict[gridding_var] = ds_in[gridding_var]
    ds_out = xr.Dataset(
        data_vars=data_var_dict,
        attrs={
            "author": author,
            "created": dt.datetime.now().astimezone().isoformat(),
            "input_year_range": f"{year_args[0]}-{year_args[1]}",
            "input_variable": variable,
        },
    )
    all_files_in_same_dir = len(np.unique([os.path.dirname(file) for file in input_files])) == 1
    if all_files_in_same_dir:
        ds_out.attrs["input_files_dir"] = os.path.dirname(input_files[0])
        ds_out.attrs["input_files"] = ", ".join([os.path.basename(file) for file in input_files])
    else:
        ds_out.attrs["input_files"] = ", ".join(input_files)
    return ds_out


def generate_gdd20_baseline(input_files, output_file, author, time_slice, variable, year_args):
    """
    Generate stream_fldFileName_gdd20_baseline file from CTSM outputs
    """

    # Define variables to process
    if variable == "GDDBX":
        suffix = "X"
    elif variable == "GDDB20":
        suffix = "20"
    else:
        raise ValueError(f"-v/--variable {variable} not recoginzed")
    var_list_in = []
    for base_temp in [0, 8, 10]:
        var_list_in.append(f"GDD{base_temp}{suffix}")

    # Get unique values and sort
    input_files = list(set(input_files))
    input_files.sort()

    # Import history files and ensure they have lat/lon dims
    ds_in = import_ds(input_files, my_vars=var_list_in + GRIDDING_VAR_LIST, time_slice=time_slice)
    if not all(x in ds_in.dims for x in ["lat", "lon"]):
        raise RuntimeError("Input files must have lat and lon dimensions")

    # If needed, find mean over time
    if "time" in ds_in.dims:
        ds_in = ds_in.mean(dim="time", skipna=True)

    # Set up a dummy DataArray to use for crops without an assigned GDDN variable
    dummy_da = xr.DataArray(
        data=np.full_like(ds_in[var_list_in[0]].values, MISSING_RX_GDD_VAL),
        dims=ds_in[var_list_in[0]].dims,
        coords=ds_in[var_list_in[0]].coords,
    )
    dummy_da = _add_time_axis(dummy_da)

    # Set up output Dataset
    ds_out = setup_output_dataset(input_files, author, variable, year_args, ds_in)

    # Process all crops
    encoding_dict = {}
    for cft_str in MGDCROP_LIST:
        cft_int = utils.vegtype_str2int(cft_str)[0]
        print(f"{cft_str} ({cft_int})")

        # Which GDDN history variable does this crop use? E.g., GDD0, GDD10
        gddn, gddn_str = _get_gddn_for_cft(cft_str, variable)

        # Fill any missing values with MISSING_RX_GDD_VAL. This will mean that gddmaturity there
        # never changes.
        if gddn_str is None:
            # Crop not handled yet? It's already filled with missing value
            this_da = dummy_da
            print("   dummy GDD20")
        else:
            this_da = ds_in[gddn_str]  # Already did ds_in.mean(dim="time") above
            this_da = _add_time_axis(this_da)
            print(f"   {gddn_str}")
            this_da = this_da.fillna(MISSING_RX_GDD_VAL)

        # Add attributes of output file
        if (gddn is None) != (gddn_str is None):
            raise RuntimeError("gddn and gddn_str must either both be None or both be not None")
        if gddn_str is None:
            long_name = "Dummy GDD20"
        else:
            long_name = f"GDD{gddn}20"
        this_da.attrs["long_name"] = long_name + f" baseline for {cft_str}"
        this_da.attrs["units"] = "°C days"

        # Copy that to ds_out
        var_out = _get_output_varname(cft_str)
        print(f"   Output variable {var_out}")
        ds_out[var_out] = this_da
        encoding_dict[var_out] = {"dtype": "float64"}

        # Grid, if needed
        if any(x not in this_da.dims for x in ["lat", "lon"]):
            ds_out[var_out] = grid_one_variable(ds_out, var_out)

    # Save
    ds_out.to_netcdf(output_file, format="NETCDF3_CLASSIC", encoding=encoding_dict)

    print("Done!")


def main():
    """
    main() function for calling generate_gdd20_baseline.py from command line.
    """
    args, time_slice = _parse_args()
    generate_gdd20_baseline(
        args.input_files,
        args.output_file,
        args.author,
        time_slice,
        args.variable,
        [args.first_year, args.last_year],
    )
