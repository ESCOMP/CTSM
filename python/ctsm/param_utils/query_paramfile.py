import argparse
import xarray as xr

PFTNAME_VAR = "pftname"


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Print values of one or more variables from a netCDF file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input netCDF file")
    parser.add_argument("variables", help="Comma-separated list of variable names to extract")
    parser.add_argument(
        "-p",
        "--pft",
        help="Comma-separated list of PFT names to print (only applies to PFT-specific variables)",
    )
    args = parser.parse_args()
    return args


def print_values(ds, var, selected_pfts, pft_names):
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
    args = get_arguments()

    variable_names = [v.strip() for v in args.variables.split(",")]

    ds = xr.open_dataset(args.input, decode_timedelta=False)

    selected_pfts = None
    pft_names = None
    if PFTNAME_VAR in ds:
        pft_names = [pft.decode().strip() for pft in ds[PFTNAME_VAR].values]
    if args.pft:
        selected_pfts = [p.strip() for p in args.pft.split(",")]
        pfts_not_in_file = []
        for pft in selected_pfts:
            if pft not in pft_names:
                pfts_not_in_file += [pft]
        if pfts_not_in_file:
            raise KeyError(f"PFT(s) not found in parameter file: {', '.join(pfts_not_in_file)}")

    for var in variable_names:
        if var in ds.variables:
            print_values(ds, var, selected_pfts, pft_names)
        else:
            print(f"Variable '{var}' not found in {args.input}")

    ds.close()


if __name__ == "__main__":
    main()
