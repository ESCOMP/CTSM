"""
Query parameters in a CTSM paramfile.

This script allows users to print the values of one or more parameters from a CTSM parameter file
(netCDF format). Users can specify variables and optionally filter output by PFTs. The script
handles multi-dimensional variables and provides formatted output for PFT-specific parameters.
"""

from ctsm.args_utils import comma_separated_list
from ctsm.param_utils.paramfile_shared import paramfile_parser_setup
from ctsm.param_utils.paramfile_shared import open_paramfile
from ctsm.param_utils.paramfile_shared import PFTNAME_VAR, check_pfts_in_paramfile
from ctsm.param_utils.paramfile_shared import get_selected_pft_indices, get_pft_names


def get_arguments():
    """
    Parse command-line arguments for querying variables from a netCDF file.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes:
            - input: Path to the netCDF file
            - variables: List of variable names to query
            - pft: Optional list of PFT names to print
    """
    parser, pft_flags = paramfile_parser_setup(
        "Print values of one or more parameters from a CTSM paramfile."
    )
    parser.add_argument(
        "variables",
        help="Names of variables to query (space-separated)",
        nargs="*",
    )
    parser.add_argument(
        *pft_flags,
        help="Comma-separated list of PFT names to print (only applies to PFT-specific variables)",
        type=comma_separated_list,
    )
    args = parser.parse_args()
    return args


def print_values(ds, var, selected_pfts, pft_names):
    """
    Print the values of a variable from the dataset, optionally filtered by PFTs.

    Parameters
    ----------
    ds : xarray.Dataset
        The opened netCDF dataset.
    var : str
        Variable name to print.
    selected_pfts : list or None
        List of selected PFT names to print, or None to print all.
    pft_names : list or None
        List of all PFT names in the file.
    """
    data = ds[var].values
    if list(ds[var].dims) == ["pft"]:
        print(var + ":")
        indices = range(len(pft_names))
        if selected_pfts is not None:
            indices = get_selected_pft_indices(selected_pfts, pft_names)
            max_name_len = max(len(pft_names[i]) for i in indices) if indices else 0
        else:
            max_name_len = max(len(name) for name in pft_names)
        for p in indices:
            print(f"   {pft_names[p]:<{max_name_len}}: {data[p]}")
    elif ds[var].ndim > 1:
        print(f"{var}:")
        print(data)
    else:
        print(f"{var}: {data}")


def main():
    """
    Main entry point for query_paramfile.

    Parses arguments, opens the netCDF file, and prints requested variable values,
    optionally filtered by PFTs.
    """
    args = get_arguments()

    ds = open_paramfile(args.input, mask_and_scale=True)

    # If user didn't specify variables, print all
    if not args.variables:
        args.variables = ds.variables

    selected_pfts = args.pft
    pft_names = None
    if selected_pfts:
        pft_names = check_pfts_in_paramfile(selected_pfts, ds)
    elif PFTNAME_VAR in ds.coords:
        pft_names = get_pft_names(ds)

    for var in args.variables:
        if var in ds.variables:
            print_values(ds, var, selected_pfts, pft_names)
        else:
            print(f"Variable '{var}' not found in {args.input}")

    ds.close()


if __name__ == "__main__":
    main()
