"""
Interpolate a maturity requirement (GDD) file
"""

import os
import sys
import argparse
import logging
import re
import xarray as xr

# -- add python/ctsm  to path (needed if we want to run this stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import ctsm_logging  # pylint: disable=wrong-import-position
from ctsm.crop_calendars.cropcal_module import (  # pylint: disable=wrong-import-position
    unexpected_negative_rx_gdd,
)

logger = logging.getLogger(__name__)

OUTPUT_FORMAT = "NETCDF3_CLASSIC"


def _file_missing(filepath, descriptor):
    if not os.path.exists(filepath) and not os.path.exists(os.path.realpath(filepath)):
        raise FileNotFoundError(f"{descriptor} not found: {filepath}")
    return os.path.realpath(filepath)


def _setup_process_args():
    """Process input arguments

    Returns:
        argparse.ArgumentParser: Arguments/options
    """

    # set up logging allowing user control
    ctsm_logging.setup_logging_pre_config()

    parser = argparse.ArgumentParser(
        description=("Interpolate a maturity requirement (GDD) file."),
    )

    # Define arguments
    parser.add_argument(
        "-i",
        "--input-file",
        help="Maturity requirement file to interpolate",
        type=lambda x: _file_missing(x, "Input file"),
        required=True,
    )
    parser.add_argument(
        "-t",
        "--target-file",
        help="File to whose coordinates the interpolation should take place",
        type=lambda x: _file_missing(x, "Target file"),
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Where to save interpolated result",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--variable-prefix",
        help="Interpolate variables whose names start with this string",
        type=str,
        required=False,
        default="gdd1_",
    )
    parser.add_argument(
        "--overwrite",
        help="If output file exists, overwrite it.",
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "--dry-run",
        help="Check arguments but do not run.",
        action="store_true",
        required=False,
    )
    ctsm_logging.add_logging_args(parser)

    # Get arguments
    args = parser.parse_args(sys.argv[1:])
    ctsm_logging.process_logging_args(args)

    # Process arguments
    if os.path.exists(os.path.realpath(args.output_file)) and not args.overwrite:
        raise FileExistsError(f"Output file exists but --overwrite not given: {args.output_file}")
    args.output_file = os.path.realpath(args.output_file)

    # Make directory for output file, if needed
    if not args.dry_run:
        parent_dir = os.path.dirname(args.output_file)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

    return args


def interpolate_gdds(args):
    """
    Do the interpolation and save
    """

    # Open inputs
    ds_in = xr.open_dataset(args.input_file)
    ds_target = xr.open_dataset(args.target_file)

    # Interpolate
    ds_out = xr.Dataset()
    for var in ds_in:

        # Check variable
        if "lat" not in ds_in[var].dims and "lon" not in ds_in[var].dims:
            print(f"Skipping variable {var} with dimensions {ds_in[var].dims}")
            continue
        if not re.compile("^" + args.variable_prefix).match(var):
            print(f"Unexpected variable {var} on input file. Skipping.")
            continue
        if args.dry_run:
            continue

        # Interpolate
        da_out = ds_in[var].interp_like(
            ds_target,
            method="nearest",
            kwargs={"fill_value": "extrapolate"},  # Otherwise you get NaNs at edges
        )

        if unexpected_negative_rx_gdd(da_out):
            raise RuntimeError("Unexpected negative value")

        # Add to dataset
        ds_out[var] = da_out

    # Finish up
    ds_out.attrs["original"] = args.input_file
    ds_out.attrs["interpolation_target"] = args.target_file
    ds_out.attrs["interpolation_script"] = os.path.basename(__file__)
    if not args.dry_run:
        ds_out.to_netcdf(args.output_file, format=OUTPUT_FORMAT)
    else:
        print("Dry run looks good!")


def main():
    """
    Description
    -----------
    Calls function that interpolates a maturity requirement (GDD) file.
    """

    args = _setup_process_args()

    interpolate_gdds(args)


if __name__ == "__main__":
    main()
