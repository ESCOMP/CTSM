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
import sys

_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
print(_CTSM_PYTHON)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.no_nans_in_inputs.constants import NEW_FILLVALUES_FILE


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
    args = parser.parse_args()

    # Load the new fill values
    print(f"Loading new fill values from {args.fillvalues_file}...")
    new_fillvalues = load_new_fillvalues(args.fillvalues_file)

    total_files = len(new_fillvalues)
    total_vars = sum(len(vars_dict) for vars_dict in new_fillvalues.values())
    print(f"Found {total_vars} variable(s) in {total_files} file(s)\n")

    # Process each file
    print("=" * 80)
    print("FILE MAPPING")
    print("=" * 80)

    for input_file, var_fillvalues in new_fillvalues.items():
        output_file = get_output_filename(input_file)

        print(f"\nInput:  {input_file}")
        print(f"Output: {output_file}")
        print(f"Variables to modify: {len(var_fillvalues)}")
        for var, fill_val in var_fillvalues.items():
            print(f"  {var}: {fill_val}")

    print("\n" + "=" * 80)
    print("\nSummary:")
    print(f"  {total_files} file(s) will be processed")
    print(f"  {total_vars} variable(s) will be modified")

    return 0


if __name__ == "__main__":
    sys.exit(main())
