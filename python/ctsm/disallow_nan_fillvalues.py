#!/usr/bin/env python3
"""
Find file paths from namelist_defaults_ctsm.xml that are also in inputdata_fillvalue.log.clm_bad.

This script:
1. Parses the XML file to extract all file paths
2. Converts relative paths (starting with lnd/clm2/) to absolute paths
3. Checks which of these paths appear in the bad files log
4. Prints the matching paths
"""

import xml.etree.ElementTree as ET
import os
import sys

import numpy as np
import xarray as xr

# File paths
XML_FILE = "bld/namelist_files/namelist_defaults_ctsm.xml"
BAD_FILES_LOG = "/glade/work/bdobbins/check_nan/inputdata_fillvalue.log"
INPUTDATA_PREFIX = "/glade/campaign/cesm/cesmdata/cseg/inputdata/"
OUR_PATH = "lnd/clm2/"  # String to be found in files we're responsible for

SEP_LENGTH = 80  # Length of horizontal separators in stdout
ATTR = "_FillValue"


def extract_file_paths_from_xml(xml_file):
    """
    Extract all file paths from the XML file.

    Args:
        xml_file: Path to the XML file

    Returns:
        set: Set of file paths found in the XML

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


def load_bad_files(bad_files_log, path_filter=None):
    """
    Load the list of bad files from the log.

    Args:
        bad_files_log: Path to the bad files log
        path_filter: Optional string that must be in the file path to include it.
                     If None, all bad files are included.

    Returns:
        set: Set of absolute file paths from the log
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


def convert_to_absolute_path(relative_path):
    """
    Convert a relative path to an absolute path.

    Args:
        relative_path: Relative path starting with OUR_PATH, or already absolute path

    Returns:
        str: Absolute path
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


def get_fill_value_from_user(var_name, target_type):
    """
    Prompt user for a new fill value and convert it to the target type.

    Args:
        var_name: Name of the variable
        target_type: Type to convert the user input to (e.g., float, int)

    Returns:
        Converted fill value of the specified type
    """
    while True:
        user_input = input(f"    New fill value for '{var_name}': ").strip()
        if user_input:
            try:
                # Convert user input to the target type
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
        else:
            print("    Please enter a value.")


def collect_new_fill_values(matches):
    """
    Interactively collect new fill values for variables with NaN fill values.

    For each file in matches, opens the file, identifies variables with NaN fill values,
    displays their properties, and prompts the user to enter new fill values.

    Args:
        matches: List of tuples (relative_path, absolute_path) for files to process

    Returns:
        dict: Dictionary mapping absolute file paths to dictionaries of
              {variable_name: new_fill_value}
    """
    print("\n" + "=" * SEP_LENGTH)
    print("COLLECTING NEW FILL VALUES")
    print("=" * SEP_LENGTH)

    # Dictionary to store new fill values for all files
    all_new_fill_values = {}

    for path_from_xml, abs_path in matches:
        print(f"\n{'=' * SEP_LENGTH}")
        print(f"Processing: {path_from_xml}")
        print(f"Full path:  {abs_path}")
        print(f"{'=' * SEP_LENGTH}")

        # Dictionary to store new fill values for this file
        new_fill_values = {}

        # Open the dataset
        ds = xr.open_dataset(abs_path, decode_cf=False, decode_timedelta=False, decode_times=False)

        # Get all variables (both data and coordinate variables)
        all_vars = list(ds.data_vars) + list(ds.coords)

        # Loop through all variables
        for var in all_vars:
            if not var_has_nan_fill(ds, var):
                continue

            da = ds[var]

            # Get variable metadata
            var_name = var
            long_name = da.attrs.get("long_name", "N/A")
            units = da.attrs.get("units", "N/A")

            # Get data statistics
            nanmin = float(np.nanmin(da.values))
            nanmax = float(np.nanmax(da.values))

            # Print variable summary
            print(f"\n  Variable: {var_name}")
            print(f"    long_name: {long_name}")
            print(f"    units:     {units}")
            print(f"    nanmin:    {nanmin}")
            print(f"    nanmax:    {nanmax}")

            # Ask user for new fill value
            new_fill_value = get_fill_value_from_user(var_name, type(nanmin))
            new_fill_values[var_name] = new_fill_value

        # Close the dataset
        ds.close()

        # Store the new fill values for this file
        if new_fill_values:
            all_new_fill_values[abs_path] = new_fill_values
            print(f"\n  Collected {len(new_fill_values)} new fill value(s) for this file:")
            for var_name, fill_val in new_fill_values.items():
                print(f"    {var_name}: {fill_val}")
        else:
            print("\n  No variables with NaN fill values found in this file.")

    return all_new_fill_values


def main():
    """Main function to find matching file paths."""

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
    all_new_fill_values = collect_new_fill_values(matches)

    return 0


if __name__ == "__main__":
    sys.exit(main())
