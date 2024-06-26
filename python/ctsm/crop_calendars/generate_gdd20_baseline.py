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

VAR_LIST_IN = ["GDD0X", "GDD8X", "GDD10X"]
GRIDDING_VAR_LIST = ["patches1d_ixy", "patches1d_jxy", "lat", "lon"]
MISSING_FILL = -1  # Something impossible to ensure that you can mark it as a missing value, to be
# bilinear-interpolated
STREAM_YEAR = 2000  # The year specified for stream_yearFirst and stream_yearLast in the call of
# shr_strdata_init_from_inline() for sdat_cropcal_gdd20_baseline


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
        help=(
            "First calendar year to include"
        ),
        type=int,
        required=False,
    )
    parser.add_argument(
        "-yN",
        "--last-year",
        help=(
            "Last calendar year to include"
        ),
        type=int,
        required=False,
    )

    # Get arguments
    args = parser.parse_args(sys.argv[1:])

    # Check arguments
    if os.path.exists(args.output_file) and not args.overwrite:
        raise FileExistsError("Output file exists but --overwrite is not specified")

    # Process time slice
    # Assumes CESM behavior where data for e.g. 1987 is saved as 1988-01-01
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

    mgdcrop_list = utils.define_mgdcrop_list()
    cft_str_list = []
    for crop_str in crop_list:
        cft_str_list += [x for x in mgdcrop_list if crop_str in x]
    return cft_str_list


def _get_gddn_for_cft(cft_str):
    """
    Given a CFT name, return the GDDN variable it uses.

    Args:
        cft_str (str): E.g., "irrigated_temperate_corn"

    Returns:
        str or None: Name of variable to use (e.g., "GDD8X"). If crop isn't yet handled, return None.
    """

    gddn = None

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
        gddn = f"GDD{gddn}X"

    return gddn


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


def generate_gdd20_baseline(input_files, output_file, author, time_slice):
    """
    Generate stream_fldFileName_gdd20_baseline file from CTSM outputs
    """

    # Get input file list
    input_files = input_files.split(sep=" ")
    # Get unique values and sort
    input_files = list(set(input_files))
    input_files.sort()

    # Import history files and ensure they have lat/lon dims
    ds_in = import_ds(input_files, VAR_LIST_IN + GRIDDING_VAR_LIST, time_slice=time_slice)
    if not all(x in ds_in.dims for x in ["lat", "lon"]):
        raise RuntimeError("Input files must have lat and lon dimensions")

    # If needed, find mean over time
    if "time" in ds_in.dims:
        ds_in = ds_in.mean(dim="time", skipna=True)

    # Set up a dummy DataArray to use for crops without an assigned GDDN variable
    dummy_da = xr.DataArray(
        data=MISSING_FILL * np.ones_like(ds_in[VAR_LIST_IN[0]].values),
        dims=ds_in[VAR_LIST_IN[0]].dims,
        coords=ds_in[VAR_LIST_IN[0]].coords,
    )
    dummy_da = _add_time_axis(dummy_da)

    # Process all crops
    data_var_dict = {}
    for v in GRIDDING_VAR_LIST:
        data_var_dict[v] = ds_in[v]
    ds_out = xr.Dataset(
        data_vars=data_var_dict,
        attrs={
            "author": author,
            "created": dt.datetime.now().astimezone().isoformat(),
        },
    )
    encoding_dict = {}
    for cft_str in utils.define_mgdcrop_list():
        cft_int = utils.vegtype_str2int(cft_str)[0]
        print(f"{cft_str} ({cft_int})")

        # Which GDDN history variable does this crop use? E.g., GDD0, GDD10
        gddn = _get_gddn_for_cft(cft_str)

        # Fill any missing values with MISSING_FILL. This will mean that gddmaturity in these cells
        # never changes.
        if gddn is None:
            # Crop not handled yet? Fill it entirely with missing value
            this_da = dummy_da
            long_name = "Dummy GDD20"
            print("   dummy GDD20")
        else:
            # this_da = ds_in[gddn].fillna(MISSING_FILL)
            this_da = ds_in[gddn]
            this_da = _add_time_axis(this_da)
            long_name = gddn.replace("X", "20")
            print(f"   {gddn}")

        # Add attributes
        this_da.attrs["long_name"] = long_name + f" baseline for {cft_str}"
        this_da.attrs["units"] = "°C days"
        # this_da.attrs["_FillValue"] = MISSING_FILL

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
        time_slice
    )
