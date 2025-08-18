"""
Tool for changing parameters on CTSM paramfile.

This script allows users to modify one or more parameters in a CTSM parameter file (netCDF format).
It supports selecting specific PFTs, dropping other PFTs, and changing parameter values, including
setting values to fill (missing) values. The script ensures safe file handling and provides
detailed error checking for argument validity and parameter changes.
"""

import os
import sys
from datetime import datetime
import numpy as np
import xarray as xr

from ctsm.args_utils import comma_separated_list
from ctsm.netcdf_utils import get_netcdf_format
from ctsm.param_utils.paramfile_shared import paramfile_parser_setup, open_paramfile
from ctsm.param_utils.paramfile_shared import check_pfts_in_paramfile, get_selected_pft_indices
from ctsm.param_utils.paramfile_shared import PFTNAME_VAR


def check_arguments(args):
    """
    Validate command-line arguments for set_paramfile.

    Checks for existence of input file, prevents overwriting output files,
    and ensures logical consistency of PFT-related options.
    """
    if not os.path.exists(args.input):
        raise FileNotFoundError(args.input)

    # Avoid potentially overwriting canonical files
    if os.path.exists(args.output):
        raise FileExistsError(args.output)

    # --drop-other-pfts makes no sense without --pfts
    if args.drop_other_pfts and not args.pft:
        raise RuntimeError("--drop-other-pfts makes no sense without -p/--pft")


def get_arguments():
    """
    Parse command-line arguments for setting variables on a netCDF file.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes:
            - input: Path to the input netCDF file
            - output: Path to the output netCDF file
            - pft: Optional list of PFT names whose values you want to change
            - drop_other_pfts: Boolean flag to drop PFTs not specified
            - param_changes: List of parameter changes to apply
    """
    parser, pft_flags = paramfile_parser_setup(
        "Change values of one or more parameters in a CTSM paramfile."
    )

    parser.add_argument(
        "-o", "--output", required=True, help="Output netCDF file. Must not already exist."
    )

    # TODO: Add mutually-exclusive --exclude-pfts argument for PFTs you DON'T want to include
    parser.add_argument(
        *pft_flags,
        help="Comma-separated list of PFTs to include (only applies to PFT-specific variables)",
        type=comma_separated_list,
    )

    parser.add_argument(
        "--drop-other-pfts",
        help=f"Do not include PFTs other than the ones given in {'/'.join(pft_flags)}",
        action="store_true",
    )

    parser.add_argument(
        "param_changes",
        help=(
            "Parameter changes to apply. Use nan to set to the fill value. E.g.:\n"
            "   param1=new_value1 pftparam=pft1_val,nan,... param3=nan"
        ),
        nargs="*",
    )

    args = parser.parse_args()
    check_arguments(args)

    return args


def is_integer(obj):
    """
    Determine if an object is an integer or a numpy array of integer dtype.

    Parameters
    ----------
    obj : object
        Object to check.

    Returns
    -------
    bool
        True if obj is an integer or numpy array of integer dtype, False otherwise.
    """
    if isinstance(obj, np.ndarray):
        obj_type = obj.dtype
    else:
        obj_type = type(obj)
    return np.issubdtype(obj_type, np.integer)


def check_correct_ndims(da, new_value, throw_error=False):
    """
    Check that the new value for a parameter has the correct number of dimensions.

    Parameters
    ----------
    da : xarray.DataArray
        The parameter DataArray to check against.
    new_value : array-like
        The new value to assign.
    throw_error : bool, optional
        If True, raise an error on mismatch.

    Returns
    -------
    bool
        True if dimensions match, False otherwise.
    """
    expected = da.ndim
    actual = np.array(new_value).ndim
    is_ndim_correct = expected == actual
    if throw_error and not is_ndim_correct:
        raise RuntimeError(f"Incorrect N dims: Expected {expected}, got {actual}")
    return is_ndim_correct


def drop_other_pfts(selected_pfts, ds):
    """
    Drop PFTs from the dataset that are not in the selected list.

    Parameters
    ----------
    selected_pfts : list of str
        List of PFT names to retain.
    ds : xarray.Dataset
        The parameter file dataset.

    Returns
    -------
    xarray.Dataset
        Dataset containing only the selected PFTs.
    """
    pft_names = check_pfts_in_paramfile(selected_pfts, ds)
    indices = get_selected_pft_indices(selected_pfts, pft_names)
    ds = ds.isel({"pft": indices})
    return ds


def _add_cmd_to_history(ds):
    """
    Prepend the calling command and timestamp to the netCDF history attribute.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to update.

    Returns
    -------
    xarray.Dataset
        Dataset with updated history attribute.
    """
    if "history" not in ds.attrs:
        ds.attrs["history"] = ""
    cmd_items = [f"'{x}'" if " " in x else x for x in sys.argv]
    datetime_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    ds.attrs["history"] = f"{datetime_str}: {' '.join(cmd_items)}\n{ds.attrs['history']}"
    return ds


def save_paramfile(ds_out: xr.Dataset, output_path, *, nc_format="NETCDF3_CLASSIC"):
    """
    Save an xarray Dataset to a netCDF parameter file.

    Parameters
    ----------
    ds_out : xarray.Dataset
        Dataset to save.
    output_path : str
        Path to output netCDF file.
    nc_format : str, optional
        NetCDF format to use (default: "NETCDF3_CLASSIC").
    """

    # We don't want to add _FillValue to parameters that didn't already have one. This dict will
    # track such parameters and be passed to .to_netcdf(..., encoding=encoding)
    encoding = {}
    for var in ds_out:
        if "_FillValue" not in ds_out[var].encoding:
            encoding[var] = {"_FillValue": None}

    ds_out = _add_cmd_to_history(ds_out)

    ds_out.to_netcdf(output_path, format=nc_format, encoding=encoding)


def _replace_nans_with_fill(var_encoding, new_value, *, chg=None):
    """
    Replace NaNs in the new parameter value array with the fill value.

    Parameters
    ----------
    var_encoding : dict
        Encoding of an xarray DataArray.
    new_value : numpy.ndarray
        Array of new values.
    chg : str, optional (keyword-only)
        Change string from command line, used for error messages.

    Returns
    -------
    numpy.ndarray
        Array with NaNs replaced by fill value.
    """
    if any(np.isnan(np.atleast_1d(new_value))):
        # TODO: Add code to add fill value to parameters without it
        if "_FillValue" not in var_encoding:
            msg = "Can't set parameter to fill value if it doesn't already have one"
            if chg is not None:
                msg += f": {chg}"
            raise NotImplementedError(msg)
        fill_value = var_encoding["_FillValue"]
        new_value[np.isnan(new_value)] = fill_value

    return new_value


def _convert_to_output_dtype(ds_out, var, new_value, *, chg=None):
    """
    Convert new_value array to the output variable's data type.

    Parameters
    ----------
    ds_out : xarray.Dataset
        Output dataset.
    var : str
        Variable name.
    new_value : numpy.ndarray
        Array of new values.
    chg : str, optional (keyword-only)
        Change string from command line, used for error messages.

    Returns
    -------
    numpy.ndarray
        Array converted to output dtype.
    """
    try:
        new_value = new_value.astype(type(ds_out[var].dtype))
    except ValueError as e:
        msg = str(e)
        if "invalid literal for int() with base 10" in msg:
            # Throw a nicer error message including the entire requested change
            err_msg = "Invalid assignment to an integer parameter"
            if chg is not None:
                err_msg += f": {chg}"
            raise ValueError(err_msg) from e
        raise e

    return new_value


def apply_new_value_to_parameter(args, ds_out, var, new_value, var_encoding, *, chg=None):
    """
    Apply a new value to a parameter in the output dataset, handling PFT selection, dimension
    checks, fill value replacement, and assignment.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    ds_out : xarray.Dataset
        Output dataset to modify.
    var : str
        Name of the variable to change.
    new_value : numpy.ndarray
        Array of new values to assign.
    var_encoding : dict
        Encoding dictionary for the variable, used for fill value replacement.
    chg : str, optional (keyword-only)
        Change string from command line, used for error messages.

    Returns
    -------
    xarray.Dataset
        Modified output dataset with the new parameter value applied.
    """
    # Are we acting on just some PFTs? If so, we'll need some stuff.
    just_some_pfts = PFTNAME_VAR in ds_out[var].coords and args.pft
    # pylint is probably wrong with the possibly-used-before-assignment warning, but do this
    # here just to placate it. Make it an invalid index so we get an error if we try to use
    # it.
    indices = -1
    if just_some_pfts:
        pft_names = check_pfts_in_paramfile(args.pft, ds_out)
        indices = get_selected_pft_indices(args.pft, pft_names)

    # Check that correct number of dimensions were given for new values. Special handling needed
    # if we're just acting on one PFT.
    da_to_check = ds_out[var]
    if just_some_pfts and len(args.pft) == 1:
        da_to_check = da_to_check.isel(pft=indices).squeeze()
    check_correct_ndims(da_to_check, new_value, throw_error=True)

    # Handle the situation where we're only changing values for some PFTs but keeping the others
    if just_some_pfts and not args.drop_other_pfts:
        tmp = ds_out[var].values.copy()
        tmp[indices] = new_value
        new_value = tmp

    # Ensure that any NaNs are replaced with the fill value
    new_value = _replace_nans_with_fill(var_encoding, new_value, chg=chg)

    # This can be needed if, e.g., you're selecting and changing just one PFT
    if ds_out[var].values.ndim > 0 and new_value.ndim == 0:
        new_value = np.atleast_1d(new_value)

    ds_out[var].values = new_value
    return ds_out


def main():
    """
    Main entry point for set_paramfile.

    Parses arguments, opens the input netCDF file, applies requested changes,
    and saves the modified dataset to a new netCDF file.
    """
    args = get_arguments()

    ds_in_masked_scaled = open_paramfile(args.input, mask_and_scale=True)

    ds_in = open_paramfile(args.input)
    ds_out = ds_in.copy()

    # If --drop-other-pfts was given, drop PFTs not in args.pft
    if args.drop_other_pfts:
        ds_out = drop_other_pfts(args.pft, ds_out)

    # Apply parameter changes, if any
    for chg in args.param_changes:
        var, new_value = chg.split("=")

        # TODO: Add handling of multi-dimensional parameters
        if ds_out[var].ndim > 1:
            raise NotImplementedError("Can't yet change multi-dimensional parameters")

        # Split at commas, if any, and convert to numpy array
        new_value = np.array(new_value.split(",")).squeeze()

        # TODO: Add code to set integer variables to their missing value. This is harder than it
        # sounds.
        if np.any(np.char.lower(new_value) == "nan") and is_integer(ds_out[var].values):
            raise NotImplementedError(f"Can't set integer parameter to fill value: {chg}")

        # Convert to the output data type
        new_value = _convert_to_output_dtype(ds_out, var, new_value, chg=chg)

        # Extract some information
        var_encoding = ds_in_masked_scaled[var].encoding

        # Apply new value to parameter
        ds_out = apply_new_value_to_parameter(args, ds_out, var, new_value, var_encoding, chg=chg)

    save_paramfile(ds_out, args.output, nc_format=get_netcdf_format(args.input))


if __name__ == "__main__":
    main()
