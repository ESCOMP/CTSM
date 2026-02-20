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


def main():
    """Main function to find matching file paths."""

    print("Extracting file paths from XML...")
    xml_paths = extract_file_paths_from_xml(XML_FILE)
    print(f"Found {len(xml_paths)} file paths in XML")

    print("\nLoading bad files from log...")
    # Only load bad files that match our path prefix to improve performance
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
            attr = "_FillValue"
            for var in ds:
                da = ds[var]
                if not attr in da.attrs:
                    continue
                if np.isnan(da.attrs[attr]):
                    any_nan_fill = True
                    break
            if not any_nan_fill:
                print(f"No variable in file has NaN {attr}; skipping")

            matches.append((path_from_xml, abs_path))
    print("-" * SEP_LENGTH)

    # Summary
    print("\nSummary:")
    print(f"  {len(xml_paths)}\tTotal paths in XML")
    print(f"  {len(bad_files)}\tTotal bad files matching '{OUR_PATH}'")
    print(f"  {len(matches)}\tMatching files with NaN {attr}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
