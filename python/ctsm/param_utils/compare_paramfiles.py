"""
Tool for comparing two CTSM paramfiles.
"""

import os
import sys
import argparse
import numpy as np
import xarray as xr

from ctsm.param_utils.paramfile_shared import open_paramfile

INDENT = "   "


def check_arguments(args) -> None:
    """
    Validate command-line arguments for compare_paramfiles.

    Checks that both input files exist and are not the same file.
    Prints a message and exits if both input paths point to the same file.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

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
        sys.exit()


def get_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for compare_paramfiles.

    Sets up argument parser for comparing two CTSM parameter files.
    Currently accepts two positional arguments for the file paths.

    Returns
    -------
    argparse.Namespace
        Parsed and validated command-line arguments with attributes:
        - file0 : str
            Path to the first parameter file
        - file1 : str
            Path to the second parameter file
    """
    parser = argparse.ArgumentParser(
        description="Compare two CTSM parameter files (netCDF)",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # TODO: Add handling of PFT list option
    # # TODO: Add mutually-exclusive --exclude-pfts argument for PFTs you DON'T want to include
    # parser.add_argument(
    #     *PFT_FLAGS,
    #     help="Comma-separated list of PFTs to include (only applies to PFT-specific variables)",
    #     type=comma_separated_list,
    #     default=None,
    # )

    # TODO: Add option to compare only a certain list of parameters

    parser.add_argument(
        "file0",
        help="Path to first paramfile",
    )

    parser.add_argument(
        "file1",
        help="Path to second paramfile",
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
            # TODO: Print the attribute value
            msg += 2 * INDENT + attr + "\n"
    attrs_in_1_not_0 = _get_attributes_in_only_one_da(da1, da0)
    if attrs_in_1_not_0:
        msg += INDENT + "Attribute(s) present in File 1 but not File 0:\n"
        for attr in attrs_in_1_not_0:
            # TODO: Print the attribute value
            msg += 2 * INDENT + attr + "\n"
    any_attrs_differ = False
    for attr in _get_attributes_in_both_da(da0, da1):
        a0 = da0.attrs[attr]
        a1 = da1.attrs[attr]
        if not np.all(a0 == a1):
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
    # TODO: Need to handle mask/scaled values that aren't nan-capable!
    try:
        values_match = values_match and np.array_equal(np0_ms, np1_ms, equal_nan=True)
    except TypeError:
        # Some data types can't have NaNs and will thus TypeError at
        # np.array_equal(..., equal_nan=True)
        pass

    # If not, loop through mismatches and add them to message
    if not values_match:
        # Find where they're unequal
        where_unequal = np.where(np0 != np1)
        if where_unequal:
            msg += INDENT + "Values differ:\n"

        for indices in zip(*where_unequal):
            msg = _one_unequal_value_msg(
                np0=np0,
                np1=np1,
                np0_ms=np0_ms,
                np1_ms=np1_ms,
                indices=indices,
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
        # TODO: If a dimension is pftname, give "name_of_pft (i)" instead of just "i"
        indices_list = str([int(i) for i in indices]) + " "

    raw_equal = v0 == v1
    ms_equal = v0_ms == v1_ms or (np.isnan(v0_ms) and np.isnan(v1_ms))

    # Raw values differ, but they're the same after masking/scaling
    if ms_equal and not raw_equal:
        msg += (
            2 * INDENT + f"{indices_list}raw: {v0} → {v1} (but identical after masking/scaling)\n"
        )

    # Raw values are the same, but they differ after masking/scaling
    elif raw_equal and not ms_equal:
        msg += 2 * INDENT + f"{indices_list}masked/scaled (raw identical): {v0_ms} → {v1_ms}\n"

    elif not raw_equal and not ms_equal:
        # Files differ in their raw and masked/scaled values, but raw difference and m/s difference
        # are the same
        if v0 == v0_ms and v1 == v1_ms:
            msg += 2 * INDENT + f"{indices_list}raw and masked/scaled: {v0} → {v1}\n"

        # Files differ in their raw and masked/scaled values, AND the raw values vs. m/s values are
        # differently different.
        else:
            msg += 2 * INDENT + f"{indices_list}raw:           {v0} → {v1}\n"
            msg += 2 * INDENT + " " * len(str(indices_list)) + f"masked/scaled: {v0_ms} → {v1_ms}\n"

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
    - Variables present in only one file
    - Attribute differences for shared variables
    - Dimension name differences
    - Data type differences
    - Shape differences
    - Value differences (both raw and masked/scaled)

    The comparison is performed on both raw values and values with netCDF
    mask/scale transformations applied (_FillValue, missing_value, scale_factor,
    add_offset).
    """
    args = get_arguments()
    check_arguments(args)

    # Print info
    print(f"File 0: {args.file0}")
    print(f"File 1: {args.file1}\n")

    # Open files
    file0_ds, file0_ms_ds = _read_paramfile(args.file0)
    file1_ds, file1_ms_ds = _read_paramfile(args.file1)

    # Check for variables only present in one dataset or the other
    vars_in_0_not_1 = _get_variables_in_only_one_ds(file0_ds, file1_ds)
    if vars_in_0_not_1:
        print("Variable(s) present in File 0 but not File 1:")
        for var in vars_in_0_not_1:
            print(INDENT + var)
        print("")
    vars_in_1_not_0 = _get_variables_in_only_one_ds(file1_ds, file0_ds)
    if vars_in_1_not_0:
        print("Variable(s) present in File 1 but not File 0:")
        for var in vars_in_1_not_0:
            print(INDENT + var)
        print("")

    # Loop through shared variables
    for var in _get_variables_in_both_ds(file0_ds, file1_ds):
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
            msg = f"{var}:\n" + msg
            print(msg)


if __name__ == "__main__":
    main()
