"""Functions for reading and manipulating namelist XML and user_nl_ files"""

from shutil import move
import glob
import os
import re
import sys
import tempfile
from typing import List
import xml.etree.ElementTree as ET

# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    INPUTDATA_PREFIX,
    ONE_OF_OUR_FILES,
    OUR_PATH,
    USERNL_NC_PATTERN,
)
from ctsm.no_nans_in_inputs.shared import (  # pylint: disable=wrong-import-position
    convert_to_absolute_path,
)


def _check_usernl_file(usernl_file) -> None:
    """Check user_nl_ file for validity"""
    # TODO: Add check of user_nl_ file validity
    # pylint: disable=unused-argument
    return


def _check_xml_file(xml_file) -> None:
    """Check XML file for validity"""
    try:
        ET.parse(xml_file)
    except IOError as e:
        print(
            f"Error: The temporary output file could not be opened: '{xml_file}'", file=sys.stderr
        )
        try:
            os.remove(xml_file)  # Because we created temp file with delete=False
        except Exception:  # pylint: disable=broad-exception-caught
            pass
        raise e
    except ET.ParseError as e:
        print("Output file is not well-formed. Contents:", file=sys.stderr)
        with open(xml_file, "r", encoding="utf8") as f:
            print(str(f.read()), file=sys.stderr)
        os.remove(xml_file)  # Because we created temp file with delete=False
        raise e
    except Exception as e:  # pylint: disable=broad-exception-caught
        os.remove(xml_file)  # Because we created temp file with delete=False
        raise e


def find_user_nl_files(dir_to_search: str) -> list[str]:
    """Find all user_nl_* files in CTSM repo"""
    pattern = os.path.join(f"{dir_to_search}", "**", "user_nl_*")
    return glob.glob(pattern, recursive=True)


def _replace_env_vars_in_netcdf_paths(netcdf_path_in: str) -> str:
    """
    Given a path to a netCDF file, replace any environment variables like $DIN_LOC_ROOT
    """
    netcdf_path_out = netcdf_path_in
    netcdf_path_out = netcdf_path_out.replace("$DIN_LOC_ROOT", INPUTDATA_PREFIX)
    netcdf_path_out = netcdf_path_out.replace("${DIN_LOC_ROOT}", INPUTDATA_PREFIX)
    return netcdf_path_out


def _extract_file_path_list_from_usernl(usernl_file: str) -> set[str]:
    """
    Extract all file paths from a user_nl file.

    Args:
        usernl_file (str): Path to the user_nl file

    Returns:
        List of file paths found in the user_nl file

    Raises:
        SystemExit: If usernl_file is not found
    """

    try:
        # Find all quoted strings in file containing OUR_PATH
        with open(usernl_file, "r", encoding="utf8") as f:
            text = f.read()
        file_paths_list = [m.group(3) for m in re.finditer(USERNL_NC_PATTERN, text, re.MULTILINE)]
    except FileNotFoundError:
        print(f"File not found: {usernl_file}", file=sys.stderr)
        sys.exit(2)
    return file_paths_list


def _extract_file_path_set_from_usernl(usernl_file: str, exact: bool = False) -> set[str]:
    """
    Extract all unique file paths from a user_nl file.

    Args:
        usernl_file (str): Path to the user_nl file
        exact (bool): Whether returned file paths should be exactly as they were in the usernl file.
                      Default False.

    Returns:
        Set of file paths found in the user_nl file

    Raises:
        SystemExit: If usernl_file is not found
    """
    file_paths = set()

    try:
        file_paths_list = _extract_file_path_list_from_usernl(usernl_file)

        # Add those strings to our set of found paths, replacing env vars if needed
        for f in file_paths_list:
            if not exact:
                f = _replace_env_vars_in_netcdf_paths(f)
            file_paths.add(f)
    except FileNotFoundError:
        print(f"File not found: {usernl_file}", file=sys.stderr)
        sys.exit(3)
    return file_paths


def _extract_file_paths_from_xml(xml_file: str) -> set[str]:
    """
    Extract all file paths from an XML file.

    Args:
        xml_file: Path to the XML file

    Returns:
        Set of file paths found in the XML

    Raises:
        SystemExit: If parsing fails or file is not found
    """
    file_paths = set()

    try:
        tree = ET.parse(xml_file)
    except ET.ParseError as parse_error:
        print(f"Error parsing XML file: {parse_error}", file=sys.stderr)
        sys.exit(4)
    except FileNotFoundError:
        print(f"XML file not found: {xml_file}", file=sys.stderr)
        sys.exit(5)

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


def extract_file_paths_from_file(file_to_search: str, exact: bool = False) -> set[str]:
    """
    Extract all file paths from a file.

    Args:
        file_to_search (str): Path to the file to search
        exact (bool): Whether returned file paths should be exactly as they were in the usernl file.
                      Default False.

    Returns:
        Set of file paths found in the file

    Raises:
        NotImplementedError: If no function exists to process this file
    """

    basename = os.path.basename(file_to_search)
    _, ext = os.path.splitext(basename)
    if ext == ".xml":
        # Doesn't need "exact" arg because no replacement happens
        file_paths = _extract_file_paths_from_xml(file_to_search)
    elif basename.startswith("user_nl"):
        file_paths = _extract_file_path_set_from_usernl(file_to_search, exact)
    else:
        raise NotImplementedError(f"Not sure how to get file paths from file: '{file_to_search}'")
    return file_paths


def how_netcdf_is_referenced_in_file(file_to_search: str, netcdf_path: str) -> List[str]:
    """
    Get list of ways a given netCDF file is referenced in a given text file

    Args:
        file_to_search (str): Path to text file we're searching
        netcdf_path (str): Path (relative or absolute) of netCDF file we're looking for

    Returns:
        Set[str]: Unique ways that netcdf_path is referenced in file_to_search
    """

    # Convert netcdf_path to absolute
    netcdf_path = convert_to_absolute_path(netcdf_path)

    # Check whether that's referenced in file_to_search
    netcdf_files_in_file = extract_file_paths_from_file(file_to_search, exact=True)
    set_of_how_this_netcdf_appears = set()
    for netcdf_file_in_file in netcdf_files_in_file:
        netcdf_file_in_file_abs = convert_to_absolute_path(
            _replace_env_vars_in_netcdf_paths(netcdf_file_in_file)
        )
        if netcdf_path == netcdf_file_in_file_abs:
            set_of_how_this_netcdf_appears = set_of_how_this_netcdf_appears | {netcdf_file_in_file}
    return set_of_how_this_netcdf_appears


def _update_xml_file(xml_file: str, old_path: str, new_path: str) -> None:
    """
    Replace a file path in an XML file.

    Replaces all occurrences of old_path with new_path throughout the file.

    Args:
        xml_file: Path to the XML file to update
        old_path: Old file path to replace (can be relative or absolute)
        new_path: New file path to use (can be relative or absolute)

    Raises:
        ValueError: If old_path not found in XML, or if XML parsing/writing fails
    """
    with open(xml_file, "r", encoding="utf-8") as f:
        xml_file_contents = f.read()

    # Regex explanation:
    #
    # <!--.*?-->
    #     Matches an entire XML comment block (non-greedy).
    #     Because we use re.DOTALL, '.' also matches newlines.
    #
    # |
    #
    # (>\s*)
    #     Group 1: matches a '>' followed by any whitespace.
    #     This ensures the path is inside a tag's text content.
    #
    # lnd/clm2/paramdata/test_params\.nc
    #     The exact filename we want to replace.
    #     The dot is escaped because '.' is special in regex.
    #
    # (\s*<)
    #     Group 2: matches optional whitespace followed by '<'.
    #     This ensures we only replace when the filename is
    #     between '>' and '<'.
    #
    # The alternation ensures that comment blocks are matched
    # as a whole and passed through unchanged.
    old_path = re.escape(old_path)
    pattern = re.compile(
        # rf"<!--.*?-->|(>\s*){old_path}(\s*<)",
        rf"<!--.*?-->|(>\s*){old_path}(\s*<)",
        re.DOTALL,
    )

    # Function to replace any quoted instances of old_path with new_path
    def replacer(match: re.Match):
        """
        If group(1) is None, the match was the comment branch.
        In that case, return the original text unchanged.

        Otherwise, reconstruct the match with NEW_PATH while
        preserving the original surrounding whitespace.
        """
        if match.group(1) is None:
            return match.group(0)  # leave comment untouched

        return f"{match.group(1)}{new_path}{match.group(2)}"

    # Perform the substitution
    new_text, n_repl = pattern.subn(replacer, xml_file_contents)
    if n_repl == 0:
        raise ValueError(f"No matches for '{old_path=}' found in '{xml_file}'")

    # Write the updated XML back to file, first stopping in a temp file for checks
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf8", delete=False) as tmp_path:
        tmp_path.write(new_text)
        xml_file_tmp = tmp_path.name
    # This will error if it's not well-formed, deleting the temporary file and causing this
    # function to stop.
    _check_xml_file(xml_file_tmp)
    move(xml_file_tmp, xml_file)
    print(f"  Updated {xml_file}")


def _update_usernl_file(usernl_file: str, old_path: str, new_path: str) -> None:
    """
    Replace a file path in a user_nl file.

    Replaces all occurrences of old_path with new_path throughout the file.

    Args:
        usernl_file: Path to the user_nl file to update
        old_path: Old file path to replace (can be relative or absolute)
        new_path: New file path to use (can be relative or absolute)
    """
    with open(usernl_file, "r", encoding="utf8") as f:
        # Get the existing file contents
        file_contents = f.read()

    # Function to replace any quoted instances of old_path with new_path
    def replacer(match):
        return (
            f"{match.group(1)}"  # Everything before opening apostrophe/quote (varname = )
            f"{match.group(2)}"  # Opening apostrophe/quote
            f"{new_path}"  # Our new path
            f"{match.group(4)}"  # Closing apostrophe/quote
            f"{match.group(5)}"  # Everything after closing apostrophe/quote (e.g., comments)
        )

    # Get the pattern we're going to replace. USERNL_NC_PATTERN is built for finding ANY of our
    # files; here, we're replacing just one specific path.
    pattern = USERNL_NC_PATTERN.replace(ONE_OF_OUR_FILES, re.escape(old_path))

    # Replace matching paths with our new one
    file_contents = re.sub(pattern, replacer, file_contents, flags=re.MULTILINE)

    # Write the updated user_nl_ back to file, first stopping in a temp file for checks
    with tempfile.NamedTemporaryFile(mode="wb", suffix=".xml", delete=False) as tmp_path:
        usernl_file_tmp = tmp_path.name
        with open(usernl_file_tmp, "w", encoding="utf8") as f:
            f.write(file_contents)
    # This will error if it's not well-formed, deleting the temporary file and causing this
    # function to stop.
    _check_usernl_file(usernl_file_tmp)
    move(usernl_file_tmp, usernl_file)


def update_text_file_referencing_netcdf(text_file: str, old_netcdf: str, new_netcdf: str) -> None:
    """
    Replace path to a netCDF file in a text file

    Args:
        text_file: Path to the text file to update
        old_path: Old file path to replace (can be relative or absolute)
        new_path: New file path to use (can be relative or absolute)
    """
    basename = os.path.basename(text_file)
    _, ext = os.path.splitext(basename)
    if ext == ".xml":
        _update_xml_file(text_file, old_netcdf, new_netcdf)
    elif basename.startswith("user_nl"):
        _update_usernl_file(text_file, old_netcdf, new_netcdf)
    else:
        raise NotImplementedError(f"Not sure how to replace file paths in file: '{text_file}'")
