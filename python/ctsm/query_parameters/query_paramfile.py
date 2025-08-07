import argparse
import xarray as xr


def get_arguments():
    parser = argparse.ArgumentParser(description="Print values of a variable from a netCDF file.")
    parser.add_argument("-i", "--input", required=True, help="Input netCDF file")
    parser.add_argument("variable", help="Name of variable to extract")
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    ds = xr.open_dataset(args.input, decode_timedelta=False)
    data = ds[args.variable].values
    print(f"{args.variable}: {data}")
    ds.close()


if __name__ == "__main__":
    main()
