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

import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR, NEW_FILLVALUES_FILE, USER_REQ_DELETE, XML_FILE


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
        str: ncatted type code (f, d, c)

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

    Raises:
        ValueError: If input and output files are the same, or if variable not found
    """
    # Ensure input and output files are different (resolve symlinks)
    input_real = os.path.realpath(input_file)
    output_real = os.path.realpath(output_file)

    if input_real == output_real:
        raise ValueError(f"Input and output files are the same: {input_file} -> {input_real}")

    # Open the input file to get actual data types
    ds = xr.open_dataset(input_file, decode_cf=False, decode_timedelta=False, decode_times=False)

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


def update_xml_file(xml_file, old_path, new_path):
    """
    Replace a file path in the XML file.

    Args:
        xml_file: Path to the XML file to update
        old_path: Old file path to replace (can be relative or absolute)
        new_path: New file path to use (can be relative or absolute)

    Raises:
        ValueError: If old_path not found in XML
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        replacements_made = 0

        # Iterate through all elements
        for elem in root.iter():
            if elem.text and old_path in elem.text:
                # Replace the old path with the new path
                elem.text = elem.text.replace(old_path, new_path)
                replacements_made += 1

        if replacements_made == 0:
            raise ValueError(f"Path '{old_path}' not found in {xml_file}")

        # Write the updated XML back to file
        tree.write(xml_file, encoding="utf-8", xml_declaration=True)
        print(f"  Updated {xml_file}: {replacements_made} replacement(s)")

    except ET.ParseError as e:
        raise ValueError(f"Error parsing XML file: {e}") from e
    except (IOError, OSError) as e:
        raise ValueError(f"Error updating XML file: {e}") from e


def process_files(fillvalues_file, dry_run=False, overwrite=False, xml_file=None):
    """
    Process files to replace fill values.

    Args:
        fillvalues_file: Path to JSON file with new fill values
        dry_run: If True, show commands without executing (default: False)
        overwrite: If True, overwrite existing output files (default: False)
        xml_file: Optional path to XML file to update with new paths (default: None)

    Returns:
        int: Number of files processed (0 if all skipped or dry-run)
    """
    # Load the new fill values
    print(f"Loading new fill values from {fillvalues_file}...")
    new_fillvalues = load_new_fillvalues(fillvalues_file)

    total_files = len(new_fillvalues)
    total_vars = sum(len(vars_dict) for vars_dict in new_fillvalues.values())
    print(f"Found {total_vars} variable(s) in {total_files} file(s)\n")

    # Process each file
    print("=" * 80)
    print("NCATTED COMMANDS")
    print("=" * 80)

    files_processed = 0

    for input_file, var_fillvalues in new_fillvalues.items():
        output_file = get_output_filename(input_file)

        # Check whether we're skipping this file
        if skip_this_file(input_file, output_file, overwrite):
            continue

        # Print things to do for this file
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
        if not dry_run:
            files_processed += execute_command(xml_file, input_file, output_file, cmd)

    # Only print summary in dry-run mode
    if dry_run:
        print_dry_run_summary(total_files, total_vars)

    return files_processed


def skip_this_file(input_file, output_file, overwrite):
    """Whether to skip the current file"""
    # Check if output is a symlink - never overwrite symlinks
    if os.path.islink(output_file):
        print(f"\n{'!' * 80}")
        print("WARNING: Output file is a symlink - SKIPPING")
        print(f"  Input:  {input_file}")
        print(f"  Output: {output_file} -> {os.readlink(output_file)}")
        print("  Symlinks will never be overwritten for safety")
        print(f"{'!' * 80}")
        return True

    # Skip if output file already exists and overwrite is not enabled
    if os.path.exists(output_file) and not overwrite:
        print(f"\nSkipping (output exists): {input_file}")
        print(f"  Output: {output_file}")
        print("  Use --overwrite to replace existing files")
        return True

    return False


def execute_command(xml_file, input_file, output_file, cmd):
    """Execute the command to make the new file and replace the old one in the XML"""
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

        # Update XML file if provided
        if xml_file:
            # Update the XML file with the new output path
            update_xml_file(xml_file, input_file, output_file)
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


def print_dry_run_summary(total_files, total_vars):
    """Print summary info about files to be processed (only in dry run)"""
    print("\n" + "=" * 80)
    print("\nSummary:")
    print(f"  {total_files} file(s) will be processed")
    print(f"  {total_vars} variable(s) will be modified")


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

    # Don't update XML in dry-run mode unless explicitly specified
    xml_file_to_update = args.xml_file if not args.dry_run else None

    # Process the files
    process_files(
        args.fillvalues_file,
        dry_run=args.dry_run,
        overwrite=args.overwrite,
        xml_file=xml_file_to_update,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
