"""
Tool for changing parameters on CTSM paramfile
"""

import os
import numpy as np

from ctsm.args_utils import comma_separated_list
from ctsm.netcdf_utils import get_netcdf_format
from ctsm.param_utils.paramfile_shared import paramfile_parser_setup, open_paramfile
from ctsm.param_utils.paramfile_shared import check_pfts_in_paramfile, get_selected_pft_indices
from ctsm.param_utils.paramfile_shared import PFTNAME_VAR


def check_arguments(args):
    """
    Check arguments to set_paramfile
    """
    if not os.path.exists(args.input):
        raise FileNotFoundError(args.input)

    # Avoid potentially overwriting canonical files
    if os.path.exists(args.output):
        raise FileExistsError(args.output)


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
    """
    parser, pft_flags = paramfile_parser_setup(
        "Change values of one or more parameters in a CTSM paramfile."
    )

    parser.add_argument("-o", "--output", required=True, help="Output netCDF file")

    # TODO: Add --exclude-pfts argument for PFTs you DON'T want to include
    parser.add_argument(
        *pft_flags,
        help="Comma-separated list of PFTs to include (only applies to PFT-specific variables)",
        type=comma_separated_list,
    )

    # TODO: Add -x/--exclude argument for variables you DON'T want to extract
    parser.add_argument(
        "-v",
        "--variables",
        help="Comma-separated list of variables to extract",
        type=comma_separated_list,
    )

    parser.add_argument(
        "param_changes",
        help="Parameter changes to apply. E.g.: param1=new_value1 pftparam=pft1_val,pft2_val,...",
        nargs="*",
    )

    args = parser.parse_args()
    check_arguments(args)

    return args


def is_integer(obj):
    """
    Given an object, return True if it's (a) any type of integer or (b) a numpy array with an
    integer dtype. Note that this will return False for integer types per se.
    """
    if isinstance(obj, np.ndarray):
        obj_type = obj.dtype
    else:
        obj_type = type(obj)
    return np.issubdtype(obj_type, np.integer)


def check_correct_ndims(da, new_value, throw_error=False):
    """
    Check that the new value given for a parameter has the right number of dimensions
    """
    expected = da.ndim
    actual = np.array(new_value).ndim
    is_ndim_correct = expected == actual
    if throw_error and not is_ndim_correct:
        raise RuntimeError(f"Incorrect N dims: Expected {expected}, got {actual}")
    return is_ndim_correct


def main():
    """
    Main entry point for the script.
    Parses arguments, opens a netCDF file, makes changes, and saves a new netCDF file.
    """
    args = get_arguments()

    ds_in = open_paramfile(args.input)
    ds_out = ds_in.copy()

    # If any PFTs were specified, drop others
    if args.pft:
        pft_names = check_pfts_in_paramfile(args.pft, ds_out)
        indices = get_selected_pft_indices(args.pft, pft_names)
        ds_out = ds_out.isel({"pft": indices})

    # If any variables specified, drop others
    if args.variables:
        vars_to_drop = []
        for var in ds_in:
            if var not in args.variables:
                vars_to_drop.append(var)
        ds_out = ds_out.drop_vars(vars_to_drop)

    # Apply parameter changes, if any
    for chg in args.param_changes:
        var, new_value = chg.split("=")

        # TODO: Add handling of multi-dimensional parameters
        if ds_out[var].ndim > 1:
            raise NotImplementedError("Can't yet change multi-dimensional parameters")

        if "," in new_value:
            new_value_list = new_value.split(",")
            new_value = np.array(new_value_list, dtype=type(ds_out[var].dtype))

        check_correct_ndims(ds_out[var], new_value, throw_error=True)

        if is_integer(ds_in[var].values):
            new_value = int(new_value)

        ds_out[var].values = new_value

    ds_out.to_netcdf(args.output, format=get_netcdf_format(args.input))


if __name__ == "__main__":
    main()
