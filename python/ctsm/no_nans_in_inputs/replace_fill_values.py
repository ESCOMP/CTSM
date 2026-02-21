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

import xarray as xr

_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
print(_CTSM_PYTHON)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.no_nans_in_inputs.constants import ATTR, NEW_FILLVALUES_FILE, USER_REQ_DELETE


def load_new_fillvalues(fillvalues_file):
    """
    Load the new fill values from JSON file.

    Args:
        fillvalues_file: Path to the JSON file with new fill values

    Returns:
        dict: Dictionary mapping file paths to variable fill values

    Raises:
        SystemExit: If file not found or invalid JSON
    """
    if not os.path.exists(fillvalues_file):
        print(f"Error: Fill values file not found: {fillvalues_file}", file=sys.stderr)
        print("Please run get_replacement_fill_values.py first.", file=sys.stderr)
        sys.exit(1)

    try:
        with open(fillvalues_file, "r", encoding="utf-8") as f:
            return json.load(f)
    except (IOError, OSError) as e:
        print(f"Error reading fill values file: {e}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON: {e}", file=sys.stderr)
        sys.exit(1)


def get_output_filename(input_file):
    """
    Generate output filename by adding .no_nan_fill before the extension.

    Args:
        input_file: Path to the input file

    Returns:
        str: Path to the output file

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


def get_ncatted_type_code(dtype):
    """
    Get ncatted type code from numpy dtype.

    Args:
        dtype: numpy dtype object

    Returns:
        str: ncatted type code (f, d, i, s, b, c)

    Raises:
        ValueError: If dtype is not recognized
    """
    dtype_str = str(dtype)

    # Float types
    if "float64" in dtype_str or "float_" in dtype_str:
        return "d"  # double
    if "float32" in dtype_str:
        return "f"  # float

    # Integer types
    if "int64" in dtype_str or "int_" in dtype_str:
        return "i"  # int (will use 32-bit in ncatted)
    if "int32" in dtype_str:
        return "i"  # int
    if "int16" in dtype_str:
        return "s"  # short
    if "int8" in dtype_str or "byte" in dtype_str:
        return "b"  # byte

    # String/char
    if "str" in dtype_str or "char" in dtype_str or "U" in dtype_str or "S" in dtype_str:
        return "c"  # char

    # Unknown type - raise error
    raise ValueError(f"Unknown dtype for ncatted: {dtype}")


def build_ncatted_command(input_file, output_file, var_fillvalues):
    """
    Build ncatted command to modify or delete fill values.

    Args:
        input_file: Path to input NetCDF file
        output_file: Path to output NetCDF file
        var_fillvalues: Dictionary mapping variable names to new fill values
                        (or USER_REQ_DELETE to delete the attribute)

    Returns:
        list: Command as list of arguments for subprocess
    """
    # Open the input file to get actual data types
    ds = xr.open_dataset(input_file, decode_cf=False, decode_timedelta=False, decode_times=False)

    cmd = ["ncatted"]

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


def main():
    """Main function to replace fill values."""

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
    args = parser.parse_args()

    # Load the new fill values
    print(f"Loading new fill values from {args.fillvalues_file}...")
    new_fillvalues = load_new_fillvalues(args.fillvalues_file)

    total_files = len(new_fillvalues)
    total_vars = sum(len(vars_dict) for vars_dict in new_fillvalues.values())
    print(f"Found {total_vars} variable(s) in {total_files} file(s)\n")

    # Process each file
    print("=" * 80)
    print("NCATTED COMMANDS")
    print("=" * 80)

    for input_file, var_fillvalues in new_fillvalues.items():
        output_file = get_output_filename(input_file)

        print(f"\nInput:  {input_file}")
        print(f"Output: {output_file}")
        print(f"Variables to modify: {len(var_fillvalues)}")
        for var, fill_val in var_fillvalues.items():
            print(f"  {var}: {fill_val}")

        # Build and print the ncatted command
        cmd = build_ncatted_command(input_file, output_file, var_fillvalues)
        print("\nCommand:")
        print("  " + " ".join(cmd))

        # Execute the command if not in dry-run mode
        if not args.dry_run:
            print("\nExecuting...")
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                print("  ✓ Success")
                if result.stdout:
                    print(f"  stdout: {result.stdout}")
                if result.stderr:
                    print(f"  stderr: {result.stderr}")
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

    # Only print summary in dry-run mode
    if args.dry_run:
        print("\n" + "=" * 80)
        print("\nSummary:")
        print(f"  {total_files} file(s) will be processed")
        print(f"  {total_vars} variable(s) will be modified")

    return 0


if __name__ == "__main__":
    sys.exit(main())
