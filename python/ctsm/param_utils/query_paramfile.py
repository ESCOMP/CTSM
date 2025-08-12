import xarray as xr

from ctsm.args_utils import comma_separated_list
from ctsm.param_utils.paramfile_shared import paramfile_parser_setup

PFTNAME_VAR = "pftname"


def get_arguments():
    """
    Parse command-line arguments for querying variables from a netCDF file.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes:
            - input: Path to the netCDF file
            - variables: Comma-separated list of variable names to extract
            - pft: Optional list of PFT names to print
    """
    parser, pft_flags = paramfile_parser_setup(
        "Print values of one or more variables from a netCDF file."
    )
    parser.add_argument(
        "variables",
        help="Comma-separated list of variable names to extract",
        type=comma_separated_list,
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
    if PFTNAME_VAR in ds[var].coords:
        print(var + ":")
        indices = range(len(pft_names))
        if selected_pfts is not None:
            indices = [i for i, name in enumerate(pft_names) if name in selected_pfts]
            max_name_len = max(len(pft_names[i]) for i in indices) if indices else 0
        else:
            max_name_len = max(len(name) for name in pft_names)
        for p in indices:
            print(f"   {pft_names[p]:<{max_name_len}}: {data[p]}")
    else:
        print(f"{var}: {data}")


def main():
    """
    Main entry point for the script.
    Parses arguments, opens the netCDF file, and prints requested variable values.
    """
    args = get_arguments()

    ds = xr.open_dataset(args.input, decode_timedelta=False)

    selected_pfts = args.pft
    pft_names = None
    if PFTNAME_VAR in ds:
        pft_names = [pft.decode().strip() for pft in ds[PFTNAME_VAR].values]
    if selected_pfts:
        pfts_not_in_file = []
        for pft in selected_pfts:
            if pft not in pft_names:
                pfts_not_in_file += [pft]
        if pfts_not_in_file:
            raise KeyError(f"PFT(s) not found in parameter file: {', '.join(pfts_not_in_file)}")

    for var in args.variables:
        if var in ds.variables:
            print_values(ds, var, selected_pfts, pft_names)
        else:
            print(f"Variable '{var}' not found in {args.input}")

    ds.close()


if __name__ == "__main__":
    main()
