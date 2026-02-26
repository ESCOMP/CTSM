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
import os
import sys

# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

from ctsm.no_nans_in_inputs.json_io import load_progress  # pylint: disable=wrong-import-position

from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    NEW_FILLVALUES_FILE,
    SEP_LENGTH,
    XML_FILE,
)
import ctsm.no_nans_in_inputs.namelist_utils as nlu  # pylint: disable=wrong-import-position
from ctsm.no_nans_in_inputs import netcdf_utils  # pylint: disable=wrong-import-position
from ctsm.no_nans_in_inputs.user_inputs import (  # pylint: disable=wrong-import-position
    confirm_continue,
)


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


def _process_one_file(
    progress: dict, input_file_abs: str, output_file: str, files_processed: list, dry_run: bool
):
    # Print things to do for this file
    var_fillvalues = progress[input_file_abs]["new_fill_values"]
    print(f"\nInput:  {input_file_abs}")
    print(f"Output: {output_file}")
    print(f"Variables to modify: {len(var_fillvalues)}")
    for var, fill_val in var_fillvalues.items():
        print(f"  {var}: {fill_val}")

    # Build and print the ncatted command
    cmd = netcdf_utils.build_ncatted_command(input_file_abs, output_file, var_fillvalues)
    print("\nCommand:")
    print("  " + " ".join(cmd))

    # Execute the command if not in dry-run mode
    if not dry_run:
        files_processed += netcdf_utils.execute_ncatted_command(cmd)
        # Update the XML file(s) with the new output path
        files_containing = []
        for file_containing_netcdf, set_of_how_this_netcdf_appears in progress[input_file_abs][
            "found_in_files"
        ].items():
            files_containing.append(file_containing_netcdf)
            for netcdf_path_in in set_of_how_this_netcdf_appears:
                netcdf_path_out = get_output_filename(netcdf_path_in)
                nlu.update_text_file_referencing_netcdf(
                    file_containing_netcdf, netcdf_path_in, netcdf_path_out
                )

        # Print message and wait for user to approve continuing
        print(
            f"Replaced in: {','.join(files_containing)}"
        )
        if not confirm_continue():
            sys.exit("Exiting.")
    return files_processed


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

        files_processed = _process_one_file(
            progress=progress,
            input_file_abs=input_file_abs,
            output_file=output_file,
            files_processed=files_processed,
            dry_run=dry_run,
        )

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
