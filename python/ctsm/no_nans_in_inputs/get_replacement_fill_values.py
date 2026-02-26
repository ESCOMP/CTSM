#!/usr/bin/env python3
"""
Find file paths from namelist_defaults_ctsm.xml that are also in inputdata_fillvalue.log.clm_bad.

This script:
1. Parses the XML file to extract all file paths
2. Converts relative paths (starting with lnd/clm2/) to absolute paths
3. Checks which of these paths appear in the bad files log
4. Prints the matching paths
"""

import argparse
import os
import sys


# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    ATTR,
    NEW_FILLVALUES_FILE,
    SEP_LENGTH,
    XML_FILE,
)
from ctsm.no_nans_in_inputs import json_io  # pylint: disable=wrong-import-position
import ctsm.no_nans_in_inputs.namelist_utils as nlu  # pylint: disable=wrong-import-position
from ctsm.no_nans_in_inputs.shared import (  # pylint: disable=wrong-import-position
    convert_to_absolute_path,
)
from ctsm.no_nans_in_inputs import user_inputs  # pylint: disable=wrong-import-position
from ctsm.no_nans_in_inputs.netcdf_utils import (  # pylint: disable=wrong-import-position
    file_has_nan_fill,
)

# File paths
DIR_TO_SEARCH_FOR_USER_NL_FILES = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir)
)


def check_write_access(file_path: str) -> bool:
    """
    Check if we have write access to create/update a file.

    Args:
        file_path: Path to the file to check

    Returns:
        True if we have write access, False otherwise
    """
    # Get the directory where the file would be created
    directory = os.path.dirname(file_path) or "."

    # Check if directory exists and is writable
    if os.path.exists(directory):
        return os.access(directory, os.W_OK)

    # If directory doesn't exist, check parent directories
    parent = os.path.dirname(directory)
    while parent and not os.path.exists(parent):
        parent = os.path.dirname(parent)

    return os.access(parent or ".", os.W_OK)


def main() -> int:
    """
    Main function to find matching file paths and collect new fill values.

    Parses command-line arguments, finds files with NaN fill values, and
    interactively collects replacement values from the user.

    Returns:
        Exit code (0 for success)
    """

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Find and fix files with NaN fill values in CTSM namelist defaults"
    )
    parser.add_argument(
        "--delete-if-none-filled",
        action="store_true",
        help="Automatically use 'delete' if variable has no filled elements (no prompt)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Print the variables that would be processed (and their defaults, if any), but don't"
            " request user input or save anything."
        ),
    )
    args = parser.parse_args()

    # Check write access to progress file before starting
    if not args.dry_run:
        print("Checking write access for progress file...")
        if not check_write_access(NEW_FILLVALUES_FILE):
            print(f"Error: No write access to create/update {NEW_FILLVALUES_FILE}", file=sys.stderr)
            dir_str = os.path.dirname(NEW_FILLVALUES_FILE) or '.'
            print(
                f"Please check permissions in directory: {dir_str}",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"✓ Write access confirmed for {NEW_FILLVALUES_FILE}\n")

    # Get list of files to search for netCDF that might have NaN fill values
    files_to_search = [XML_FILE]
    files_to_search.extend(nlu.find_user_nl_files(DIR_TO_SEARCH_FOR_USER_NL_FILES))

    netcdf_paths = set()
    files_referencing_netcdfs = []
    for file_to_search in files_to_search:
        print(f"Extracting file paths from file: {file_to_search}")
        netcdf_paths_thisfile = nlu.extract_file_paths_from_file(file_to_search)
        print(f"Found {len(netcdf_paths_thisfile)} file paths in file")
        if netcdf_paths_thisfile:
            files_referencing_netcdfs.append(file_to_search)
        netcdf_paths = netcdf_paths | netcdf_paths_thisfile

    # Load existing progress if available
    progress = json_io.init_progress()

    print("\nFinding matches...")

    files_not_found = []
    for netcdf_path in sorted(netcdf_paths):
        print(f"Finding matches for: {netcdf_path}")

        # Check that the file exists
        abs_path = convert_to_absolute_path(netcdf_path)
        if not os.path.exists(abs_path):
            # TODO: Actually handle files that weren't found, if possible.
            files_not_found.append(abs_path)
            continue
        # TODO: Check that the file is in CESM inputdata dir

        print(f"Does exist, abs path: {abs_path}")
        print("-" * SEP_LENGTH)
        print(f"In XML/user_nl file:   {netcdf_path}")
        print(f"Absolute: {abs_path}")

        # Check that the file actually has NaN _FillValue for at least one var
        any_nan_fill, vars_with_nan_fills = file_has_nan_fill(abs_path)
        if any_nan_fill:
            if abs_path not in progress:
                progress[abs_path] = json_io.create_empty_progress_dict_onefile()
            fif_dict = progress[abs_path]["found_in_files"]
            for file_to_search in files_referencing_netcdfs:
                set_of_how_this_netcdf_appears = nlu.how_netcdf_is_referenced_in_file(
                    file_to_search, netcdf_path
                )
                if set_of_how_this_netcdf_appears:
                    if file_to_search not in fif_dict:
                        fif_dict[file_to_search] = set()
                    fif_dict[file_to_search] = (
                        fif_dict[file_to_search] | set_of_how_this_netcdf_appears
                    )
            progress[abs_path]["vars_with_nan_fills"] = vars_with_nan_fills
            json_io.save_progress(progress.copy(), NEW_FILLVALUES_FILE)
        else:
            if abs_path in progress:
                raise RuntimeError(
                    f"Found no NaN fills in file but it was in progress dict: {abs_path}"
                )
            print(f"No variable in file has NaN {ATTR}; skipping")

    print("-" * SEP_LENGTH)

    # Summary
    print("\nSummary:")
    print(f"  {len(netcdf_paths)}\tTotal paths in XML and user_nl_ files")
    print(f"  {len(progress)}\tFiles with NaN {ATTR}")
    print(f"  {len(files_not_found)}\tFiles not found")
    if files_not_found:
        for f in files_not_found:
            print(f"\t* Not found: '{f}'")

    # Collect new fill values from user
    user_inputs.collect_new_fill_values(
        progress, delete_if_none_filled=args.delete_if_none_filled, dry_run=args.dry_run
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
