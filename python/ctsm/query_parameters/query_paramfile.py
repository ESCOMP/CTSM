import argparse
import xarray as xr

PFTNAME_VAR = "pftname"


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Print values of one or more variables from a netCDF file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input netCDF file")
    parser.add_argument("variables", help="Comma-separated list of variable names to extract")
    args = parser.parse_args()
    return args


def print_values(ds, var):
    data = ds[var].values
    if PFTNAME_VAR in ds[var].coords:
        print(var + ":")
        pft_names = [pft.decode().strip() for pft in ds[PFTNAME_VAR].values]
        max_name_len = max(len(name) for name in pft_names)
        for p, pft_name in enumerate(pft_names):
            print(f"   {pft_name:<{max_name_len}}: {data[p]}")
    else:
        print(f"{var}: {data}")


def main():
    args = get_arguments()

    variable_names = [v.strip() for v in args.variables.split(",")]

    ds = xr.open_dataset(args.input, decode_timedelta=False)

    for var in variable_names:
        if var in ds.variables:
            print_values(ds, var)
        else:
            print(f"Variable '{var}' not found in {args.input}")

    ds.close()


if __name__ == "__main__":
    main()
