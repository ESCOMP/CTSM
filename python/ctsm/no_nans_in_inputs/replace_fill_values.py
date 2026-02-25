#!/usr/bin/env python3
"""
Replace NaN fill values in NetCDF files based on new_fillvalues.json.

This script:
1. Reads the new_fillvalues.json file created by get_replacement_fill_values.py
2. For each file, creates an output filename with .no_nan_fill before the extension
3. Uses ncatted to modify or delete _FillValue attributes
4. Creates modified copies of the input files
"""

import argparse
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from typing import Any
import re

# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

import numpy as np  # pylint: disable=wrong-import-position
import xarray as xr  # pylint: disable=wrong-import-position

from ctsm.no_nans_in_inputs.json_io import load_progress  # pylint: disable=wrong-import-position

from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    ATTR,
    NEW_FILLVALUES_FILE,
    OPEN_DS_KWARGS,
    SEP_LENGTH,
    USER_REQ_DELETE,
    XML_FILE,
)
import ctsm.no_nans_in_inputs.namelist_utils as nlu  # pylint: disable=wrong-import-position


def get_output_filename(input_file: str) -> str:
    """
    Generate output filename by adding .no_nan_fill before the extension.

    Args:
        input_file: Path to the input file

    Returns:
        Path to the output file

    Examples:
        /path/to/file.nc -> /path/to/file.no_nan_fill.nc
        /path/to/file.tar.gz -> /path/to/file.no_nan_fill.tar.gz
    """
    # Split the path into directory, basename, and extension
    directory = os.path.dirname(input_file)
    basename = os.path.basename(input_file)

    # Find the last dot to split extension
    if "." in basename:
        name, ext = basename.rsplit(".", 1)
        output_basename = f"{name}.no_nan_fill.{ext}"
    else:
        # No extension
        output_basename = f"{basename}.no_nan_fill"

    # Reconstruct the full path
    return os.path.join(directory, output_basename)


def get_ncatted_type_code(dtype: np.dtype) -> str:
    """
    Get ncatted type code from numpy dtype.

    Args:
        dtype: numpy dtype object

    Returns:
        ncatted type code (f, d, c)

    Raises:
        ValueError: If dtype is not recognized or is an integer type
                    (NetCDF doesn't allow NaN fill values for integers)
    """
    dtype_str = str(dtype)

    # Float types
    if "float64" in dtype_str or "float_" in dtype_str:
        return "d"  # double
    if "float32" in dtype_str:
        return "f"  # float

    # Integer types - not allowed (NetCDF doesn't support NaN for integers)
    if any(x in dtype_str for x in ["int64", "int32", "int16", "int8", "int_", "byte"]):
        raise ValueError(
            f"Integer dtype detected: {dtype}. "
            "NetCDF does not allow NaN fill values for integer variables. So how'd this happen?"
        )

    # String/char
    if "str" in dtype_str or "char" in dtype_str or "U" in dtype_str or "S" in dtype_str:
        return "c"  # char

    # Unknown type - raise error
    raise ValueError(f"Unknown dtype for ncatted: {dtype}")


def build_ncatted_command(
    input_file: str, output_file: str, var_fillvalues: dict[str, Any]
) -> list[str]:
    """
    Build ncatted command to modify or delete fill values.

    Args:
        input_file: Path to input NetCDF file
        output_file: Path to output NetCDF file
        var_fillvalues: Dictionary mapping variable names to new fill values
                        (or USER_REQ_DELETE to delete the attribute)

    Returns:
        Command as list of arguments for subprocess

    Raises:
        ValueError: If input and output files are the same, or if variable not found
    """
    # Ensure input and output files are different (resolve symlinks)
    input_real = os.path.realpath(input_file)
    output_real = os.path.realpath(output_file)

    if input_real == output_real:
        raise ValueError(f"Input and output files are the same: {input_file} -> {input_real}")

    # Open the input file to get actual data types
    ds = xr.open_dataset(input_file, **OPEN_DS_KWARGS)

    cmd = ["ncatted", "-O"]  # -O flag to overwrite without prompting

    for var, fill_val in var_fillvalues.items():
        if fill_val == USER_REQ_DELETE:
            # Delete the attribute: -a attr_name,var_name,d,,
            cmd.extend(["-a", f"{ATTR},{var},d,,"])
        else:
            # Get the actual data type from the file
            if var in ds.data_vars:
                dtype = ds[var].dtype
            elif var in ds.coords:
                dtype = ds[var].dtype
            else:
                # Variable not found - raise error
                ds.close()
                raise ValueError(f"Variable '{var}' not found in {input_file}")

            # Get the appropriate type code for ncatted
            type_code = get_ncatted_type_code(dtype)

            # Modify the attribute: -a attr_name,var_name,o,type,value
            cmd.extend(["-a", f"{ATTR},{var},o,{type_code},{fill_val}"])

    # Close the dataset
    ds.close()

    # Add input and output files
    cmd.extend([input_file, output_file])

    return cmd


def process_files(
    fillvalues_file: str,
    dry_run: bool = False,
    overwrite: bool = False,
) -> int:
    """
    Process files to replace fill values.

    Args:
        fillvalues_file: Path to JSON file with new fill values
        dry_run: If True, show commands without executing (default: False)
        overwrite: If True, overwrite existing output files (default: False)

    Returns:
        Number of files successfully processed
    """
    # Load the new fill values
    print(f"Loading new fill values from {fillvalues_file}...")
    progress = load_progress(fillvalues_file)

    total_files = len(progress)
    total_vars = sum(len(vars_dict) for vars_dict in progress.values())
    print(f"Found {total_vars} variable(s) in {total_files} file(s)\n")

    # Process each file
    print("=" * SEP_LENGTH)
    print("NCATTED COMMANDS")
    print("=" * SEP_LENGTH)

    files_processed = 0

    for input_file_abs in progress:
        output_file = get_output_filename(input_file_abs)

        # Check whether we're skipping this file
        if skip_this_file(input_file_abs, output_file, overwrite):
            continue

        # Print things to do for this file
        var_fillvalues = progress[input_file_abs]["new_fill_values"]
        print(f"\nInput:  {input_file_abs}")
        print(f"Output: {output_file}")
        print(f"Variables to modify: {len(var_fillvalues)}")
        for var, fill_val in var_fillvalues.items():
            print(f"  {var}: {fill_val}")

        # Build and print the ncatted command
        cmd = build_ncatted_command(input_file_abs, output_file, var_fillvalues)
        print("\nCommand:")
        print("  " + " ".join(cmd))

        # Execute the command if not in dry-run mode
        if not dry_run:
            files_processed += execute_command(cmd)
            # Update the XML file(s) with the new output path
            for file_containing_netcdf, set_of_how_this_netcdf_appears in progress[input_file_abs][
                "found_in_files"
            ].items():
                for netcdf_path_in in set_of_how_this_netcdf_appears:
                    netcdf_path_out = get_output_filename(netcdf_path_in)
                    nlu.update_text_file_referencing_netcdf(
                        file_containing_netcdf, netcdf_path_in, netcdf_path_out
                    )
            # TODO: git commit!

    # Only print summary in dry-run mode
    if dry_run:
        print_dry_run_summary(total_files, total_vars)

    return files_processed


def skip_this_file(input_file: str, output_file: str, overwrite: bool) -> bool:
    """
    Determine whether to skip processing a file.

    Files are skipped if the output is a symlink (always) or if the output
    exists and overwrite is False.

    Args:
        input_file: Path to input file
        output_file: Path to output file
        overwrite: Whether to overwrite existing files

    Returns:
        True if file should be skipped, False otherwise
    """
    # Check if output is a symlink - never overwrite symlinks
    if os.path.islink(output_file):
        print(f"\n{'!' * SEP_LENGTH}")
        print("WARNING: Output file is a symlink - SKIPPING")
        print(f"  Input:  {input_file}")
        print(f"  Output: {output_file} -> {os.readlink(output_file)}")
        print("  Symlinks will never be overwritten for safety")
        print(f"{'!' * SEP_LENGTH}")
        return True

    # Skip if output file already exists and overwrite is not enabled
    if os.path.exists(output_file) and not overwrite:
        print(f"\nSkipping (output exists): {input_file}")
        print(f"  Output: {output_file}")
        print("  Use --overwrite to replace existing files")
        return True

    return False


def execute_command(cmd: list[str]) -> int:
    """
    Runs the ncatted command to create the output file with modified fill values.

    Args:
        cmd: ncatted command as list of arguments

    Returns:
        Number of files processed (1 on success, 0 on skip)

    Raises:
        SystemExit: If ncatted command fails or is not found
    """
    print("\nExecuting...")
    files_processed = 0
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("  ✓ Success")
        if result.stdout:
            print(f"  stdout: {result.stdout}")
        if result.stderr:
            print(f"  stderr: {result.stderr}")
        files_processed = 1

    except subprocess.CalledProcessError as e:
        print(f"  ✗ Error: ncatted failed with exit code {e.returncode}", file=sys.stderr)
        if e.stdout:
            print(f"  stdout: {e.stdout}", file=sys.stderr)
        if e.stderr:
            print(f"  stderr: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("  ✗ Error: ncatted command not found", file=sys.stderr)
        print("  Please ensure NCO (NetCDF Operators) is installed", file=sys.stderr)
        sys.exit(1)
    return files_processed


def print_dry_run_summary(total_files: int, total_vars: int) -> None:
    """
    Print summary information about files to be processed.

    Only called in dry-run mode to show what would be done.

    Args:
        total_files: Total number of files to process
        total_vars: Total number of variables to modify
    """
    print("\n" + "=" * SEP_LENGTH)
    print("\nSummary:")
    print(f"  {total_files} file(s) will be processed")
    print(f"  {total_vars} variable(s) will be modified")


def main() -> int:
    """
    Main function to replace fill values.

    Parses command-line arguments and processes files to replace NaN fill values.

    Returns:
        Exit code (0 for success)
    """

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Replace NaN fill values in NetCDF files using ncatted"
    )
    parser.add_argument(
        "--fillvalues-file",
        default=NEW_FILLVALUES_FILE,
        help=f"Path to JSON file with new fill values (default: {NEW_FILLVALUES_FILE})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without actually modifying files",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files (default: skip if output exists)",
    )
    parser.add_argument(
        "--xml-file",
        default=XML_FILE,
        help=f"Path to XML file to update with new paths (default: {XML_FILE})",
    )
    args = parser.parse_args()

    # Process the files
    process_files(
        args.fillvalues_file,
        dry_run=args.dry_run,
        overwrite=args.overwrite,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
