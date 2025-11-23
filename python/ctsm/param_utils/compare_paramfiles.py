"""
Tool for comparing two CTSM paramfiles.
"""

import os
import argparse
import warnings
import numpy as np
import xarray as xr

from ctsm.param_utils.paramfile_shared import open_paramfile, get_pft_names
from ctsm.param_utils.paramfile_shared import are_paramfile_dataarrays_identical
from ctsm.netcdf_utils import get_netcdf_format
from ctsm.args_utils import comma_separated_list

INDENT = "   "


def check_arguments(args) -> bool:
    """
    Validate command-line arguments for compare_paramfiles.

    Checks that both input files exist and are not the same file.
    Prints a message if both input paths point to the same file.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    bool
        True if comparison should proceed, False otherwise.

    Raises
    ------
    FileNotFoundError
        If either file0 or file1 does not exist.
    """
    if not os.path.exists(args.file0):
        raise FileNotFoundError(args.file0)
    if not os.path.exists(args.file1):
        raise FileNotFoundError(args.file1)

    if os.path.realpath(args.file0) == os.path.realpath(args.file1):
        print("These are the same file.")
        return False
    return True


def get_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for compare_paramfiles.

    Sets up argument parser for comparing two CTSM parameter files.
    Accepts two positional arguments for file paths and optional parameter filtering.

    Returns
    -------
    argparse.Namespace
        Parsed and validated command-line arguments with attributes:
        - file0 : str
            Path to the first parameter file
        - file1 : str
            Path to the second parameter file
        - params : list[str]
            List of parameter names to compare (empty list means compare all)
    """
    parser = argparse.ArgumentParser(
        description="Compare two CTSM parameter files (netCDF)",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # TODO: Add handling of PFT list option. Actually, make it more general:
    # * --dim-[dimname], --index-[dimname], --indices-[dimname], --indexes-[dimname] should all do
    #   isel along the given dim. Help text should indicate that these are 0-based.
    # * --coord-[coordname], --value-[coordname], --values-[coordname] should all do sel along the
    #   given coord.
    # * -p/--pft (i.e., PFT_FLAGS) as shortcut for --dim-pft and --coord-pftname. Decide whether
    #   user meant dimension (index) or coordinate (value) based on whether they gave integers or
    #   names, respectively.

    parser.add_argument(
        "file0",
        help="Path to first paramfile",
    )

    parser.add_argument(
        "file1",
        help="Path to second paramfile",
    )

    parser.add_argument(
        "--params",
        "--param",
        "--parameters",
        help="Comma-separated list of parameters to compare",
        type=comma_separated_list,
        default=[],
    )

    args = parser.parse_args()

    # TODO: Add option to ignore attribute differences
    # TODO: Add option to ignore raw differences if same after mask/scale
    # TODO: Add threshold option

    return args


def _read_paramfile(paramfile_path: str) -> tuple[xr.Dataset, xr.Dataset]:
    """
    Read a parameter file with and without mask/scale applied.

    Opens the same parameter file twice: once with mask_and_scale=False to preserve
    raw values, and once with mask_and_scale=True to apply _FillValue, missing_value,
    scale_factor, and add_offset transformations.

    Parameters
    ----------
    paramfile_path : str
        Path to the parameter file to read.

    Returns
    -------
    tuple of xr.Dataset
        A tuple containing:
        - paramfile_ds : xr.Dataset
            Dataset with raw values (mask_and_scale=False)
        - paramfile_ms_ds : xr.Dataset
            Dataset with masked and scaled values (mask_and_scale=True)
    """
    paramfile_ds = open_paramfile(paramfile_path, mask_and_scale=False)
    paramfile_ms_ds = open_paramfile(paramfile_path, mask_and_scale=True)
    return paramfile_ds, paramfile_ms_ds


def _get_variables_in_only_one_ds(ds_a: xr.Dataset, ds_b: xr.Dataset) -> list[str]:
    """
    Returns a sorted list of variables that exist only in the first Dataset.

    Parameters
    ----------
    ds_a : xr.Dataset
        First dataset to compare.
    ds_b : xr.Dataset
        Second dataset to compare.

    Returns
    -------
    list of str
        Sorted list of variable names present in ds_a but not in ds_b.
    """
    in_a_not_b = []
    for var in ds_a.variables:
        if var not in ds_b.variables:
            in_a_not_b.append(var)
    in_a_not_b.sort()
    return in_a_not_b


def _print_variables_in_only_one_ds(
    header: str, vars_in_only_one: list[str], any_diffs: bool, args_params: list[str]
) -> bool:
    """
    Print variables that exist in only one dataset, filtered by parameter list.

    Prints a header followed by variable names, but only if there are variables to report
    and they match the parameter filter (if provided).

    Parameters
    ----------
    header : str
        Header text to print before the variable list.
    vars_in_only_one : list[str]
        List of variable names that exist in only one dataset.
    any_diffs : bool
        Whether any differences have been found so far.
    args_params : list[str]
        List of parameters to filter by. If empty, all variables are included.

    Returns
    -------
    bool
        Updated any_diffs flag - True if any variables were printed, otherwise unchanged.
    """
    if vars_in_only_one:
        for var in vars_in_only_one:
            if args_params and var not in args_params:
                continue
            if not any_diffs:
                any_diffs = True
                print(header)
            print(INDENT + var)
        print("")
    return any_diffs


def _get_variables_in_both_ds(ds_a: xr.Dataset, ds_b: xr.Dataset) -> list[str]:
    """
    Returns a sorted list of variables that exist in both Datasets.

    Parameters
    ----------
    ds_a : xr.Dataset
        First dataset to compare.
    ds_b : xr.Dataset
        Second dataset to compare.

    Returns
    -------
    list of str
        Sorted list of variable names present in both ds_a and ds_b.
    """
    in_a_and_b = []
    for var in ds_a.variables:
        if var in ds_b.variables:
            in_a_and_b.append(var)
    in_a_and_b.sort()
    return in_a_and_b


def _get_attributes_in_only_one_da(da_a: xr.DataArray, da_b: xr.DataArray) -> list[str]:
    """
    Returns a sorted list of attributes that exist only in the first DataArray.

    Parameters
    ----------
    da_a : xr.DataArray
        First data array to compare.
    da_b : xr.DataArray
        Second data array to compare.

    Returns
    -------
    list of str
        Sorted list of attribute names present in da_a but not in da_b.
    """
    in_a_not_b = []
    for attr in list(da_a.attrs.keys()):
        if attr not in list(da_b.attrs.keys()):
            in_a_not_b.append(attr)
    in_a_not_b.sort()
    return in_a_not_b


def _get_attributes_in_both_da(da_a: xr.DataArray, da_b: xr.DataArray) -> list[str]:
    """
    Returns a sorted list of attributes that exist in both DataArrays.

    Parameters
    ----------
    da_a : xr.DataArray
        First data array to compare.
    da_b : xr.DataArray
        Second data array to compare.

    Returns
    -------
    list of str
        Sorted list of attribute names present in both da_a and da_b.
    """
    in_a_and_b = []
    for attr in list(da_a.attrs.keys()):
        if attr in list(da_b.attrs.keys()):
            in_a_and_b.append(attr)
    in_a_and_b.sort()
    return in_a_and_b


def _compare_attrs(da0: xr.DataArray, da1: xr.DataArray, msg: str) -> str:
    """
    Compare attributes for two DataArrays and append differences to message.

    Identifies attributes that are present in only one DataArray and attributes
    with different values between the two DataArrays.

    Parameters
    ----------
    da0 : xr.DataArray
        First data array to compare.
    da1 : xr.DataArray
        Second data array to compare.
    msg : str
        Existing message string to which comparison results will be appended.

    Returns
    -------
    str
        Updated message string with attribute comparison results appended.
    """
    attrs_in_0_not_1 = _get_attributes_in_only_one_da(da0, da1)
    if attrs_in_0_not_1:
        msg += INDENT + "Attribute(s) present in File 0 but not File 1:\n"
        for attr in attrs_in_0_not_1:
            msg += f"{2*INDENT} {attr}: {da0.attrs[attr]}\n"
    attrs_in_1_not_0 = _get_attributes_in_only_one_da(da1, da0)
    if attrs_in_1_not_0:
        msg += INDENT + "Attribute(s) present in File 1 but not File 0:\n"
        for attr in attrs_in_1_not_0:
            msg += f"{2*INDENT} {attr}: {da1.attrs[attr]}\n"
    any_attrs_differ = False
    for attr in _get_attributes_in_both_da(da0, da1):
        a0 = da0.attrs[attr]
        a1 = da1.attrs[attr]
        shapes_differ = np.array(a0).shape != np.array(a1).shape
        if shapes_differ:
            values_differ = True
        else:
            try:
                values_differ = not np.all((a0 == a1) | (np.isnan(a0) & np.isnan(a1)))
            except TypeError:
                values_differ = not np.all(a0 == a1)
        if values_differ:
            if not any_attrs_differ:
                msg += INDENT + "Attribute(s) with different values:\n"
                any_attrs_differ = True
            msg += INDENT * 2 + f"{attr}, file 0: {a0}\n"
            msg += INDENT * 2 + f"{attr}, file 1: {a1}\n"
    return msg


def _compare_da_values(
    da0_ms: xr.DataArray, da1_ms: xr.DataArray, da0: xr.DataArray, da1: xr.DataArray, msg: str
) -> str:
    """
    Compare values between two DataArrays and append differences to message.

    Compares both raw values and masked/scaled values, identifying where they differ
    and providing detailed information about the nature of the differences.

    Parameters
    ----------
    da0_ms : xr.DataArray
        First data array with mask/scale applied.
    da1_ms : xr.DataArray
        Second data array with mask/scale applied.
    da0 : xr.DataArray
        First data array with raw values (no mask/scale).
    da1 : xr.DataArray
        Second data array with raw values (no mask/scale).
    msg : str
        Existing message string to which comparison results will be appended.

    Returns
    -------
    str
        Updated message string with value comparison results appended.
    """
    # Ensure at least one-dimensional to enable comparing singleton arrays with scalars
    np0 = np.atleast_1d(da0)
    np1 = np.atleast_1d(da1)
    np0_ms = np.atleast_1d(da0_ms)
    np1_ms = np.atleast_1d(da1_ms)

    # Do values match?
    values_match = np.array_equal(np0, np1)
    try:
        values_match = values_match and np.array_equal(np0_ms, np1_ms, equal_nan=True)
    except TypeError:
        # Some data types will TypeError at np.array_equal(..., equal_nan=True)
        values_match = values_match and np.array_equal(np0_ms, np1_ms)

    # If not, loop through mismatches and add them to message
    if not values_match:
        # Find where they're unequal
        where_unequal = np.where(np0 != np1)
        if where_unequal:
            msg += INDENT + "Values differ:\n"

        # Get dimension names if they match between the two DataArrays
        dimnames = list(da0.dims) if da0.dims == da1.dims else None

        # Get pftnames if they match between the two DataArrays
        pftnames = None
        if dimnames:
            try:
                pftnames0 = get_pft_names(da0)
                pftnames1 = get_pft_names(da1)
                if pftnames0 == pftnames1:
                    pftnames = pftnames0
            except KeyError:
                # pftname coordinate not found on one or both DataArrays
                pass

        for indices in zip(*where_unequal):
            msg = _one_unequal_value_msg(
                np0=np0,
                np1=np1,
                np0_ms=np0_ms,
                np1_ms=np1_ms,
                indices=indices,
                dimnames=dimnames,
                pftnames=pftnames,
                msg=msg,
            )
    return msg


def _one_unequal_value_msg(
    *,
    np0: np.ndarray,
    np1: np.ndarray,
    np0_ms: np.ndarray,
    np1_ms: np.ndarray,
    indices: tuple,
    dimnames: list[str] | None,
    pftnames: list[str] | None,
    msg: str,
) -> str:
    """
    Generate a message describing a single unequal value between two arrays.

    Creates a detailed message showing how values differ at a specific index,
    including both raw and masked/scaled values when they differ.

    Parameters
    ----------
    np0 : np.ndarray
        First array with raw values.
    np1 : np.ndarray
        Second array with raw values.
    np0_ms : np.ndarray
        First array with mask/scale applied.
    np1_ms : np.ndarray
        Second array with mask/scale applied.
    indices : tuple
        Indices where the values differ.
    dimnames : list[str] | None
        List of dimension names if they match between arrays, None otherwise.
    pftnames : list[str] | None
        List of PFT names if they match between arrays, None otherwise.
    msg : str
        Existing message string to which the difference description will be appended.

    Returns
    -------
    str
        Updated message string with the value difference description appended.

    Raises
    ------
    RuntimeError
        If the function is called when values are actually equal (should not occur).
    """
    v0 = np0[indices]
    v1 = np1[indices]
    v0_ms = np0_ms[indices]
    v1_ms = np1_ms[indices]

    # List of indices for each unequal element in the array, or blank if just one element/scalar.
    # Will be put at the front of the message line for the element.
    if np.atleast_1d(np0).size == 1:  # We've already checked that the arrays are the same size.
        indices_list = ""
    else:
        # Format indices with dimension names if available
        if dimnames:
            indices_parts = []
            for i, (dimname, idx) in enumerate(zip(dimnames, indices)):
                idx_int = int(idx)
                # If this dimension is pft/pftname and we have pftnames, include the PFT name
                if pftnames and dimname.lower() in ["pft", "pftname"]:
                    pft_name = pftnames[idx_int]
                    indices_parts.append(f"{dimname} {idx_int} ({pft_name})")
                else:
                    indices_parts.append(f"{dimname} {idx_int}")
            indices_list = "[" + ", ".join(indices_parts) + "]"
        else:
            indices_list = str([int(i) for i in indices])

    raw_equal = v0 == v1

    # Check if masked/scaled values are equal, handling NaN comparison
    # Some data types (like strings) can't use np.isnan, so handle TypeError
    ms_equal = v0_ms == v1_ms
    try:
        ms_equal = ms_equal or (np.isnan(v0_ms) and np.isnan(v1_ms))
    except TypeError:
        pass

    # Raw values differ, but they're the same after masking/scaling
    if ms_equal and not raw_equal:
        msg += (
            2 * INDENT
            + f"{indices_list} raw: {v0} → {v1} (raw; both {v0_ms} after masking/scaling)\n"
        )

    # Raw values are the same, but they differ after masking/scaling
    elif raw_equal and not ms_equal:
        msg += 2 * INDENT + f"{indices_list} {v0_ms} → {v1_ms} (masked/scaled; raw both {v0})\n"

    elif not raw_equal and not ms_equal:
        # Files differ in their raw and masked/scaled values, but raw difference and m/s difference
        # are the same
        if v0 == v0_ms and v1 == v1_ms:
            msg += 2 * INDENT + f"{indices_list} {v0} → {v1}\n"

        # Files differ in their raw and masked/scaled values, AND the raw values vs. m/s values are
        # differently different.
        else:
            msg += 2 * INDENT + f"{indices_list} {v0} → {v1} (raw)\n"
            msg += (
                2 * INDENT + " " * len(str(indices_list)) + f" {v0_ms} → {v1_ms} (masked/scaled)\n"
            )

    # This shouldn't be reachable, because it means neither raw nor m/s values differ
    else:
        raise RuntimeError("How?")

    return msg


def compare_da_shapes(da0: xr.DataArray, da1: xr.DataArray, msg: str) -> tuple:
    """
    Compare shapes of two DataArrays and append differences to message.

    Parameters
    ----------
    da0 : xr.DataArray
        First data array to compare.
    da1 : xr.DataArray
        Second data array to compare.
    msg : str
        Existing message string to which comparison results will be appended.

    Returns
    -------
    tuple of (str, bool)
        A tuple containing:
        - msg : str
            Updated message string with shape comparison results appended.
        - shapes_match : bool
            True if shapes match, False otherwise.
    """
    shape0 = np.atleast_1d(da0).shape
    shape1 = np.atleast_1d(da1).shape
    shapes_match = shape0 == shape1
    if not shapes_match:
        msg += INDENT + f"Shapes differ: File 0: {list(shape0)}\n"
        msg += INDENT + f"               File 1: {list(shape1)}\n"
    return msg, shapes_match


def _compare_da_dtypes(da0: xr.DataArray, da1: xr.DataArray, msg: str) -> str:
    """
    Compare data types of two DataArrays and append differences to message.

    Parameters
    ----------
    da0 : xr.DataArray
        First data array to compare.
    da1 : xr.DataArray
        Second data array to compare.
    msg : str
        Existing message string to which comparison results will be appended.

    Returns
    -------
    str
        Updated message string with data type comparison results appended.
    """
    dtype0 = da0.values.dtype
    dtype1 = da1.values.dtype
    if dtype0 != dtype1:
        msg += INDENT + f"Data types differ: File 0: {dtype0}"
        msg += INDENT + f"                   File 1: {dtype1}"
    return msg


def _compare_da_dimnames(da0: xr.DataArray, da1: xr.DataArray, msg: str) -> str:
    """
    Compare dimension names of two DataArrays and append differences to message.

    Parameters
    ----------
    da0 : xr.DataArray
        First data array to compare.
    da1 : xr.DataArray
        Second data array to compare.
    msg : str
        Existing message string to which comparison results will be appended.

    Returns
    -------
    str
        Updated message string with dimension name comparison results appended.
    """
    if da0.dims != da1.dims:
        msg += INDENT + f"Dimension names differ: File 0: {list(da0.dims)}\n"
        msg += INDENT + f"                        File 1: {list(da1.dims)}\n"
    return msg


def main():
    """
    Main entry point for compare_paramfiles.

    Compares two CTSM parameter files and prints detailed differences including:
    - File type (netCDF format) differences
    - Variables present in only one file
    - Attribute differences for shared variables
    - Dimension name differences
    - Data type differences
    - Shape differences
    - Value differences (both raw and masked/scaled)

    The comparison is performed on both raw values and values with netCDF
    mask/scale transformations applied (_FillValue, missing_value, scale_factor,
    add_offset).

    Returns early if the files are the same (based on realpath comparison).
    Prints "Files are identical" if no differences are found.

    Command-line arguments can filter which parameters to compare using --params.
    """
    args = get_arguments()
    if not check_arguments(args):
        return
    any_diffs = False

    # Print info
    print(f"File 0: {args.file0}")
    print(f"File 1: {args.file1}\n")

    # Check whether file types differ
    file0_type = get_netcdf_format(args.file0)
    file1_type = get_netcdf_format(args.file1)
    file_types_match = file0_type == file1_type
    if not file_types_match:
        any_diffs = True
        print(f"File types differ: {file0_type} → {file1_type}\n")

    # Open files
    file0_ds, file0_ms_ds = _read_paramfile(args.file0)
    file1_ds, file1_ms_ds = _read_paramfile(args.file1)

    # TODO: Check global attributes

    # Check for variables only present in one dataset or the other
    any_diffs = _print_variables_in_only_one_ds(
        "Variable(s) present in File 0 but not File 1:",
        _get_variables_in_only_one_ds(file0_ds, file1_ds),
        any_diffs,
        args.params,
    )
    any_diffs = _print_variables_in_only_one_ds(
        "Variable(s) present in File 1 but not File 0:",
        _get_variables_in_only_one_ds(file1_ds, file0_ds),
        any_diffs,
        args.params,
    )

    # Get shared parameters, restricting to those in --params if given
    params_to_check = _get_variables_in_both_ds(file0_ds, file1_ds)
    if args.params:
        # Check that all requested parameters are in at least one file
        vars_in_neither = []
        for var in args.params:
            if var not in file0_ds and var not in file1_ds:
                vars_in_neither.append(var)
        if vars_in_neither:
            raise KeyError(f"Requested parameter(s) in neither file: {', '.join(vars_in_neither)}")
        params_to_check = [var for var in params_to_check if var in args.params]

    # Loop through parameters
    for var in params_to_check:
        msg = ""

        # Masked and scaled version: _FillValue/missing_value, scale_factor, and add_offset applied
        da0_ms = file0_ms_ds[var]
        da1_ms = file1_ms_ds[var]

        # _FillValue/missing_value, scale_factor, and add_offset NOT applied
        da0 = file0_ds[var]
        da1 = file1_ds[var]

        # Check attributes. Use NON-masked-and-scaled so as to check missing value etc. too.
        msg = _compare_attrs(da0, da1, msg)

        # Check that dimension names are the same
        msg = _compare_da_dimnames(da0, da1, msg)

        # Check that data types are the same
        msg = _compare_da_dtypes(da0, da1, msg)

        # Check that shapes are the same
        msg, shapes_match = compare_da_shapes(da0, da1, msg)

        # Check values (only if shapes match)
        if shapes_match:
            msg = _compare_da_values(da0_ms, da1_ms, da0, da1, msg)

        if msg:
            any_diffs = True
            msg = f"{var}:\n" + msg
            print(msg)
        elif file_types_match:
            # In case we missed something here. Since this is sort of just a bonus, don't do this
            # check if file types don't match, because are_paramfile_dataarrays_identical() will
            # say the DataArrays also don't match no matter what.
            identical = are_paramfile_dataarrays_identical(da0, da1)
            identical_ms = are_paramfile_dataarrays_identical(da0_ms, da1_ms)
            if not (identical and identical_ms):
                warnings.warn(
                    f"compare_paramfiles thinks {var} DataArrays are identical but"
                    " are_paramfile_dataarrays_identical() disagrees."
                )

    if not any_diffs:
        print("Files are identical.")


if __name__ == "__main__":
    main()
