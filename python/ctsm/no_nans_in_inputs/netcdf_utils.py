"""Functions for working with netCDF files"""

import os
import sys
import re
import subprocess
from typing import List, Tuple

import numpy as np
import xarray as xr

# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    ATTR,
    USER_REQ_DELETE,
)

from ctsm.no_nans_in_inputs.shared import (  # pylint: disable=wrong-import-position
    FillValueConfig,
    VarContext,
)
from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    VARSTARTS_TO_DEFAULT_NEG999,
)


def file_has_nan_fill(abs_path: str) -> Tuple[bool, List[str]]:
    """
    Check if a netCDF file has any variable with NaN fill value attribute.

    Args:
        abs_path: Absolute path to file

    Returns:
        bool: True if the file has any variable with NaN fill value attribute, False otherwise
        List[str]: Variables with NaN fill value attributes
    """
    vars_with_nan_fills = get_vars_with_nan_fills(abs_path)
    return bool(vars_with_nan_fills), vars_with_nan_fills


def get_var_info(
    var: str, ds: xr.Dataset, abs_path: str, delete_if_none_filled: bool, dry_run: bool
) -> Tuple[VarContext, FillValueConfig]:
    """
    Process a single variable to get information to be used as settings.

    Displays variable metadata and statistics and calculates a smart default.

    Args:
        var: Variable name
        da: xarray DataArray for the variable
        abs_path: Absolute path to the file (for context in defaults)
        delete_if_none_filled: If True, automatically use delete when it's the default
        dry_run: If true, just print vars to process (and defaults, if any).

    Returns:
        VarContext: Information about the variable
        FillValueConfig: Information about the fill value
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

    # Save and return info
    var_context = VarContext(
        var_name=var, target_type=type(nanmin), file_path=abs_path, dry_run=dry_run
    )
    config = FillValueConfig(
        default_value=default_fill,
        allow_delete=not data_has_nan,
        delete_if_none_filled=delete_if_none_filled,
    )

    return var_context, config


def get_vars_with_nan_fills(abs_path: str) -> List[str]:
    """
    Given a file, get variables with NaN fill value attribute (if any).

    Args:
        abs_path: Absolute path to file

    Returns:
        bool: List of variables with NaN fill value attribute
    """
    ncdump_results = subprocess.check_output(["ncdump", "-h", abs_path], text=True)

    # Regex breakdown:
    # ^\s* : Start of line and any leading whitespace
    # (\S+)       : Capture one or more non-whitespace characters (the variable name)
    # :{ATTR}} : The attribute where fill value is stored
    # \s*=\s* : The equals sign with flexible surrounding whitespace
    # NaNf?\s*;    : The NaN/NaNf value and the closing semicolon
    regex_pattern = rf"^\s*(\S+):{ATTR}\s*=\s*NaNf?\s*;"

    # Use re.MULTILINE to treat each line in the string as a new start
    vars_with_nan_fills = re.findall(regex_pattern, ncdump_results, re.MULTILINE)
    vars_with_nan_fills.sort()
    return vars_with_nan_fills


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
