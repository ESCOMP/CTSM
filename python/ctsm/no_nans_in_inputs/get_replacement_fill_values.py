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
import re
import xml.etree.ElementTree as ET
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from typing import Any

import numpy as np
import xarray as xr

from ctsm.no_nans_in_inputs.constants import (
    ATTR,
    NEW_FILLVALUES_FILE,
    USER_REQ_DELETE,
    USER_REQ_QUIT,
    USER_REQ_SKIP_FILE,
    USER_REQ_SKIP_VAR,
    XML_FILE,
)

# File paths
BAD_FILES_LOG = "/glade/work/bdobbins/check_nan/inputdata_fillvalue.log"
INPUTDATA_PREFIX = "/glade/campaign/cesm/cesmdata/cseg/inputdata/"
OUR_PATH = "lnd/clm2/"  # String to be found in files we're responsible for
PROGRESS_FILE = NEW_FILLVALUES_FILE  # Alias for clarity in this script

SEP_LENGTH = 80  # Length of horizontal separators in stdout

VARSTARTS_TO_DEFAULT_NEG999 = ["fertl_", "irrig_", "crpbf_", "fharv_"]


@dataclass
class VarContext:
    """Context about the variable being processed.

    Attributes:
        var_name: Name of the variable
        target_type: Type to convert user input to (e.g., float, int)
        file_path: Optional path to the netCDF file (for ncdump on Ctrl-C)
    """

    var_name: str
    target_type: type
    file_path: str | None = None


@dataclass
class FillValueConfig:
    """Configuration for how the fill value prompt should behave.

    Attributes:
        default_value: Optional default value to use if user presses enter
        allow_delete: Whether to allow deleting the fill value attribute
        delete_if_none_filled: If True, automatically use delete when it's the default
    """

    default_value: Any = None
    allow_delete: bool = True
    delete_if_none_filled: bool = False


def extract_file_paths_from_xml(xml_file: str) -> set[str]:
    """
    Extract all file paths from the XML file.

    Args:
        xml_file: Path to the XML file

    Returns:
        Set of file paths found in the XML

    Raises:
        SystemExit: If XML parsing fails or file is not found
    """
    file_paths = set()

    try:
        tree = ET.parse(xml_file)
    except ET.ParseError as parse_error:
        print(f"Error parsing XML file: {parse_error}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"XML file not found: {xml_file}", file=sys.stderr)
        sys.exit(1)

    root = tree.getroot()

    # Iterate through all elements in the XML
    for elem in root.iter():
        # Get the text content of the element
        if elem.text and elem.text.strip():
            text = elem.text.strip()
            # Check if it looks like a file path (contains OUR_PATH)
            if OUR_PATH in text:
                # Extract just the path part (in case there's other text)
                # Split by whitespace and look for the path
                for token in text.split():
                    if OUR_PATH in token:
                        file_paths.add(token)

    return file_paths


def load_bad_files(bad_files_log: str, path_filter: str | None = None) -> set[str]:
    """
    Load the list of bad files from the log.

    Args:
        bad_files_log: Path to the bad files log
        path_filter: Optional string that must be in the file path to include it.
                     If None, all bad files are included.

    Returns:
        Set of absolute file paths from the log

    Raises:
        SystemExit: If file not found
    """
    bad_files = set()
    bad_line_contents = " : NaN_FillValue : "

    try:
        with open(bad_files_log, "r", encoding="utf-8") as f:
            for line in f:
                # Each line starts with the file path followed by " : NaN_FillValue : "
                if bad_line_contents in line:
                    file_path = line.split(bad_line_contents)[0].strip()
                    # Only add files that match our filter (if specified)
                    if path_filter is None or path_filter in file_path:
                        bad_files.add(file_path)

    except FileNotFoundError:
        print(f"Bad files log not found: {bad_files_log}", file=sys.stderr)
        sys.exit(1)

    return bad_files


def convert_to_absolute_path(relative_path: str) -> str:
    """
    Convert a relative path to an absolute path.

    Args:
        relative_path: Relative path starting with OUR_PATH, or already absolute path

    Returns:
        Absolute path
    """
    # If the path is already absolute, return it as-is
    if os.path.isabs(relative_path):
        return relative_path

    # Otherwise, convert relative path to absolute
    return os.path.join(INPUTDATA_PREFIX, relative_path)


def var_has_nan_fill(ds: xr.Dataset, var: str, attr: str = ATTR) -> bool:
    """
    Check if a variable has a NaN fill value attribute.

    Args:
        ds: xarray Dataset containing the variable
        var: Name of the variable to check
        attr: Name of the attribute to check (typically '_FillValue')

    Returns:
        bool: True if the variable has the specified attribute and its value is NaN,
              False otherwise
    """
    da = ds[var]
    if not attr in da.attrs:
        return False
    try:
        result = np.isnan(da.attrs[attr])
    except TypeError:
        return False
    return result


def var_data_has_nan(da: xr.DataArray) -> bool:
    """
    Check if a variable's data contains any NaN values.

    Args:
        da: xarray DataArray to check

    Returns:
        bool: True if the data contains any NaN values, False otherwise
    """
    try:
        return bool(da.isnull().any())
    except TypeError:
        # If isnan fails (e.g., for string data), assume no NaN
        return False


def show_ncdump_for_variable(file_path: str | None, var_name: str) -> None:
    """
    Run ncdump -h on a file and display lines matching the variable name.

    Args:
        file_path: Path to the netCDF file (None to skip)
        var_name: Name of the variable to search for in ncdump output
    """
    if not file_path:
        print("    No file path available for ncdump")
        print()
        return

    try:
        print(f"    Running: ncdump -h {file_path}")
        result = subprocess.run(
            ["ncdump", "-h", file_path], capture_output=True, text=True, check=True
        )
        # Filter lines containing the variable name
        matching_lines = [line for line in result.stdout.split("\n") if var_name in line]
        if matching_lines:
            print(f"    Lines matching '{var_name}':")
            for line in matching_lines:
                print(f"      {line}")
        else:
            print(f"    No lines found matching '{var_name}'")
    except subprocess.CalledProcessError as e:
        print(f"    Error running ncdump: {e}")
    except FileNotFoundError:
        print("    Error: ncdump command not found")

    print()  # Empty line for readability


def get_var_info(var: str, ds: xr.Dataset, abs_path: str, delete_if_none_filled: bool) -> Any:
    """
    Process a single variable to get its new fill value from the user.

    Displays variable metadata and statistics, calculates a smart default,
    and prompts the user for a new fill value.

    Args:
        var: Variable name
        da: xarray DataArray for the variable
        abs_path: Absolute path to the file (for context in defaults)
        delete_if_none_filled: If True, automatically use delete when it's the default

    Returns:
        New fill value (number or USER_REQ_DELETE string)

    Raises:
        ValueError: If user chooses to skip this variable or file
    """
    da = ds[var]

    # Get variable metadata
    long_name = da.attrs.get("long_name", "N/A")
    units = da.attrs.get("units", "N/A")
    shape = da.shape

    # Get data statistics
    nanmin = float(np.nanmin(da.values))
    nanmax = float(np.nanmax(da.values))

    # Check if data contains any NaN values
    data_has_nan = var_data_has_nan(da)

    # Calculate default fill value
    default_fill = None
    # Suggest delete if data has no NaN values
    if not data_has_nan:
        default_fill = USER_REQ_DELETE
    elif (
        nanmin >= 0
        or nanmin == -1
        or any(var.startswith(x) for x in VARSTARTS_TO_DEFAULT_NEG999)
        or ("/surfdata_map/" in abs_path and bool(re.match(r"[a-z0-9]{5}_to_[a-z0-9]{5}", var)))
    ):
        default_fill = type(nanmin)(-999)

    # Print variable summary
    print(f"\n  Variable: {var}")
    print(f"    long_name: {long_name}")
    print(f"    shape:     {shape}")
    print(f"    units:     {units}")
    print(f"    nanmin:    {nanmin}")
    print(f"    nanmax:    {nanmax}")
    if data_has_nan:
        print(f"    WARNING: Data contains NaN values - cannot delete {ATTR}")

    # Ask user for new fill value
    var_context = VarContext(var_name=var, target_type=type(nanmin), file_path=abs_path)
    config = FillValueConfig(
        default_value=default_fill,
        allow_delete=not data_has_nan,
        delete_if_none_filled=delete_if_none_filled,
    )
    new_fill_value = get_fill_value_from_user(var_context, config)

    return new_fill_value


def _handle_special_command(user_input: str, allow_delete: bool) -> Any | None:
    """
    Check if user input is a special command and handle it.

    Args:
        user_input: The user's input string (already stripped)
        allow_delete: Whether deleting the fill value is allowed

    Returns:
        USER_REQ_DELETE if delete was requested and allowed, or None if not a special command

    Raises:
        KeyboardInterrupt: If user typed 'quit'
        ValueError: If user typed 'skip' or 'skipfile'
    """
    lower_input = user_input.lower()
    if lower_input == USER_REQ_QUIT:
        raise KeyboardInterrupt("User requested quit")
    if lower_input == USER_REQ_SKIP_VAR:
        raise ValueError("SKIP_VARIABLE")
    if lower_input == USER_REQ_SKIP_FILE:
        raise ValueError("SKIP_FILE")
    if lower_input == USER_REQ_DELETE:
        if not allow_delete:
            print(f"    Error: Cannot delete {ATTR} - variable contains NaN values")
            return None
        print(f"    Will delete {ATTR} attribute")
        return USER_REQ_DELETE
    return None


def _convert_and_validate_input(user_input: str, target_type: type) -> Any | None:
    """
    Convert user input to the target type and validate it's not NaN.

    Args:
        user_input: The user's input string (already stripped)
        target_type: Type to convert the input to

    Returns:
        Converted value, or None if conversion failed (error message already printed)
    """
    try:
        converted_value = target_type(user_input)

        # Make sure it's not NaN
        try:
            converted_value_is_nan = np.isnan(converted_value)
        except TypeError:
            converted_value_is_nan = False
        if converted_value_is_nan:
            raise ValueError(f"Input '{user_input}' would produce a NaN {ATTR}")

        return converted_value
    except (ValueError, TypeError) as e:
        print(f"    Invalid input: {e}. Please enter a valid {target_type.__name__}.")
        return None


def _handle_empty_input(default_value: Any, allow_delete: bool) -> Any | None:
    """
    Handle empty user input by returning the default or printing help.

    Args:
        default_value: Default value to use, or None if no default
        allow_delete: Whether deleting the fill value is allowed (for help message)

    Returns:
        The default value if one exists, or None if no default (help message already printed)
    """
    if default_value is not None:
        print(f"    Using default: {default_value}")
        return default_value

    # Build help message based on what's allowed
    options = []
    if allow_delete:
        options.append(f"'{USER_REQ_DELETE}' to delete attribute")
    options.extend(
        [
            f"'{USER_REQ_SKIP_VAR}' to skip variable",
            f"'{USER_REQ_SKIP_FILE}' to skip file",
            f"'{USER_REQ_QUIT}' to save and exit",
        ]
    )
    print(f"    Please enter a value (or {', '.join(options)}).")
    return None


def _handle_ctrl_c(ctrl_c_count: int, user_input: str | None, var_context: VarContext) -> int:
    """
    Handle a KeyboardInterrupt (Ctrl-C) during user input.

    On the first Ctrl-C, shows ncdump output for the variable.
    On the second Ctrl-C (or if user had typed 'quit'), re-raises.

    Args:
        ctrl_c_count: Number of times Ctrl-C has been pressed (before incrementing)
        user_input: The user's input before Ctrl-C (may be None)
        var_context: Context about the variable being processed

    Returns:
        Updated ctrl_c_count

    Raises:
        KeyboardInterrupt: If this is the second Ctrl-C or user typed 'quit'
    """
    ctrl_c_count += 1

    # If this is the second Ctrl-C or the user requested quit, exit
    if ctrl_c_count >= 2 or user_input == USER_REQ_QUIT:
        if ctrl_c_count >= 2:
            print("\n    [Ctrl-C pressed again - exiting]")
        else:
            print("\n    User requested quit")
        raise KeyboardInterrupt

    # First Ctrl-C: show ncdump output for this variable
    print("\n    [Ctrl-C detected - press again to exit]")
    show_ncdump_for_variable(var_context.file_path, var_context.var_name)
    return ctrl_c_count


def get_fill_value_from_user(var_context: VarContext, config: FillValueConfig) -> Any:
    """
    Prompt user for a new fill value and convert it to the target type.

    Args:
        var_context: Context about the variable (name, type, file path)
        config: Configuration for prompt behavior (default, allow_delete, auto-delete)

    Returns:
        Converted fill value of the specified type, or USER_REQ_DELETE string

    Raises:
        KeyboardInterrupt: If user presses Ctrl-C twice or types 'quit'
        ValueError: If user types 'skip' or 'skipfile'
    """
    # If delete_if_none_filled is enabled and default is delete, use it automatically
    if config.delete_if_none_filled and config.default_value == USER_REQ_DELETE:
        print(f"    Auto-deleting {ATTR} attribute, since no elements are filled")
        return USER_REQ_DELETE

    ctrl_c_count = 0

    while True:
        user_input = None
        try:
            # Build prompt with default value if available
            if config.default_value is not None:
                prompt = f"    New fill value for '{var_context.var_name}' [default: {config.default_value}]: "
            else:
                prompt = f"    New fill value for '{var_context.var_name}': "

            user_input = input(prompt).strip()

            if user_input:
                # Check for special commands first
                special_result = _handle_special_command(user_input, config.allow_delete)
                if special_result is not None:
                    return special_result

                # Try to convert to the target type
                converted = _convert_and_validate_input(user_input, var_context.target_type)
                if converted is not None:
                    return converted
            else:
                # Empty input - use default or show help
                empty_result = _handle_empty_input(config.default_value, config.allow_delete)
                if empty_result is not None:
                    return empty_result

        except KeyboardInterrupt:
            ctrl_c_count = _handle_ctrl_c(ctrl_c_count, user_input, var_context)


def collect_new_fill_values(
    matches: list[tuple[str, str]],
    progress_file: str = PROGRESS_FILE,
    delete_if_none_filled: bool = False,
) -> dict[str, dict[str, Any]]:
    """
    Interactively collect new fill values for variables with NaN fill values.

    For each file in matches, opens the file, identifies variables with NaN fill values,
    displays their properties, and prompts the user to enter new fill values.

    Progress is automatically saved after each variable. User can type 'quit' to save and exit,
    or 'skip' to skip a variable.

    Args:
        matches: List of tuples (relative_path, absolute_path) for files to process
        progress_file: Path to save/load progress (default: PROGRESS_FILE)
        delete_if_none_filled: If True, automatically use delete when it's the default

    Returns:
        Dictionary mapping absolute file paths to dictionaries of {variable_name: new_fill_value}
    """
    print("\n" + "=" * SEP_LENGTH)
    print("COLLECTING NEW FILL VALUES")
    print("=" * SEP_LENGTH)

    # Load existing progress if available
    all_new_fill_values = load_progress(progress_file)
    if all_new_fill_values:
        print(f"\nLoaded progress from {progress_file}")
        total_vars = sum(len(vars_dict) for vars_dict in all_new_fill_values.values())
        print(f"Already processed {total_vars} variable(s) in {len(all_new_fill_values)} file(s)")
        response = input("Continue from where you left off? [Y/n]: ").strip().lower()
        if response and response not in ("y", "yes"):
            all_new_fill_values = {}
            print("Starting fresh...")

    print(
        f"\nCommands: Type a number for fill value, '{USER_REQ_DELETE}' to delete attribute, "
        f"'{USER_REQ_SKIP_VAR}' to skip variable, '{USER_REQ_SKIP_FILE}' to skip file, "
        f"'{USER_REQ_QUIT}' to save and exit"
    )

    try:
        for path_from_xml, abs_path in matches:
            print(f"\n{'=' * SEP_LENGTH}")
            print(f"Processing: {path_from_xml}")
            print(f"Full path:  {abs_path}")
            print(f"{'=' * SEP_LENGTH}")

            # Get or create dictionary for this file's fill values
            if abs_path not in all_new_fill_values:
                all_new_fill_values[abs_path] = {}
            new_fill_values = all_new_fill_values[abs_path]

            # Open the dataset
            ds = xr.open_dataset(
                abs_path, decode_cf=False, decode_timedelta=False, decode_times=False
            )

            # Get all variables (both data and coordinate variables)
            all_vars = list(ds.data_vars) + list(ds.coords)

            # Loop through all variables
            for var in all_vars:
                if not var_has_nan_fill(ds, var):
                    continue

                # Skip variables we've already processed
                if var in new_fill_values:
                    print(f"\n  Variable: {var} [already processed, skipping]")
                    continue

                new_fill_value = get_var_info(var, ds, abs_path, delete_if_none_filled)

                # Process this variable to get new fill value
                try:
                    new_fill_values[var] = new_fill_value

                    # Save progress after each variable
                    save_progress(all_new_fill_values, progress_file)
                except ValueError as e:
                    # Check if this is the skip variable signal
                    if str(e) == "SKIP_VARIABLE":
                        print(f"    Skipping variable '{var}'")
                        continue
                    # Check if this is the skip file signal
                    if str(e) == "SKIP_FILE":
                        print("    Skipping rest of file")
                        break
                    # Otherwise re-raise
                    raise

            # Close the dataset
            ds.close()

            # Print summary for this file
            print(f"\n  Collected {len(new_fill_values)} new fill value(s) for this file:")
            for var, fill_val in new_fill_values.items():
                print(f"    {var}: {fill_val}")

    except KeyboardInterrupt:
        print("Exiting.")
        sys.exit(0)

    return all_new_fill_values


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


def save_progress(all_new_fill_values: dict[str, dict[str, Any]], progress_file: str) -> None:
    """
    Save progress to a JSON file.

    Args:
        all_new_fill_values: Dictionary of collected fill values
        progress_file: Path to the progress file
    """
    try:
        with open(progress_file, "w", encoding="utf-8") as f:
            json.dump(all_new_fill_values, f, indent=2)
        print(f"  [Progress saved to {progress_file}]")
    except (IOError, OSError) as e:
        print(f"  Warning: Could not save progress: {e}", file=sys.stderr)


def load_progress(progress_file: str) -> dict[str, dict[str, Any]]:
    """
    Load progress from a JSON file if it exists.

    Args:
        progress_file: Path to the progress file

    Returns:
        Previously saved fill values, or empty dict if file doesn't exist
    """
    if not os.path.exists(progress_file):
        return {}

    try:
        with open(progress_file, "r", encoding="utf-8") as f:
            return json.load(f)
    except (IOError, OSError, json.JSONDecodeError) as e:
        print(f"Warning: Could not load progress file: {e}", file=sys.stderr)
        return {}


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
    args = parser.parse_args()

    # Check write access to progress file before starting
    print("Checking write access for progress file...")
    if not check_write_access(PROGRESS_FILE):
        print(f"Error: No write access to create/update {PROGRESS_FILE}", file=sys.stderr)
        print(
            f"Please check permissions in directory: {os.path.dirname(PROGRESS_FILE) or '.'}",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"✓ Write access confirmed for {PROGRESS_FILE}\n")

    print("Extracting file paths from XML...")
    xml_paths = extract_file_paths_from_xml(XML_FILE)
    print(f"Found {len(xml_paths)} file paths in XML")

    print("\nLoading bad files from log...")
    bad_files = load_bad_files(BAD_FILES_LOG, path_filter=OUR_PATH)
    print(f"Found {len(bad_files)} bad files in log matching '{OUR_PATH}'")

    print("\nFinding matches...")
    matches = []

    for path_from_xml in sorted(xml_paths):
        abs_path = convert_to_absolute_path(path_from_xml)
        if abs_path in bad_files:
            # Check that the file exists
            if not os.path.exists(abs_path):
                raise FileNotFoundError(abs_path)

            print("-" * SEP_LENGTH)
            print(f"In XML:   {path_from_xml}")
            print(f"Absolute: {abs_path}")

            # Check that the file actually has NaN _FillValue for at least one var
            ds = xr.open_dataset(
                abs_path, decode_cf=False, decode_timedelta=False, decode_times=False
            )
            any_nan_fill = False
            for var in ds:
                if var_has_nan_fill(ds, var):
                    any_nan_fill = True
                    break
            if not any_nan_fill:
                print(f"No variable in file has NaN {ATTR}; skipping")

            matches.append((path_from_xml, abs_path))
    print("-" * SEP_LENGTH)

    # Summary
    print("\nSummary:")
    print(f"  {len(xml_paths)}\tTotal paths in XML")
    print(f"  {len(bad_files)}\tTotal bad files matching '{OUR_PATH}'")
    print(f"  {len(matches)}\tMatching files with NaN {ATTR}")

    # Collect new fill values from user
    all_new_fill_values = collect_new_fill_values(
        matches, delete_if_none_filled=args.delete_if_none_filled
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
