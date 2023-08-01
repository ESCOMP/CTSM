import numpy as np
import xarray as xr
import sys
import argparse
import shutil
import os


def main(file_in, dir_out):

    # Checks and setup
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    if not os.path.isdir(args.output_dir):
        raise RuntimeError("-o/--output-dir needs to be a directory")
    path, ext = os.path.splitext(file_in)
    dir_in, filename_in_noext = os.path.split(path)
    file_out = os.path.join(dir_out, f"{filename_in_noext}.all_crops_everywhere{ext}")

    # Import
    print(f"Importing {file_in}...")
    in_ds = xr.open_dataset(file_in)

    pct_crop_da = in_ds["PCT_CROP"]
    pct_natveg_da = in_ds["PCT_NATVEG"]
    pct_cft_da_in = in_ds["PCT_CFT"]
    in_ds.close()

    # Move all natural land into crop
    pct_crop_da += pct_natveg_da
    pct_natveg_da -= pct_natveg_da

    # Put some of every crop in every gridcell
    pct_cft = np.full_like(pct_cft_da_in.values, 100 / in_ds.dims["cft"])
    pct_cft_da = xr.DataArray(
        data=pct_cft,
        attrs=pct_cft_da_in.attrs,
        dims=pct_cft_da_in.dims,
        name="PCT_CFT",
    )

    # Save
    print(f"Saving to {file_out}")
    shutil.copyfile(file_in, file_out)
    format = "NETCDF3_64BIT"
    mode = "a"  # Use existing file but overwrite existing variables
    pct_crop_da.to_netcdf(file_out, format=format, mode=mode)
    pct_natveg_da.to_netcdf(file_out, format=format, mode=mode)
    pct_cft_da.to_netcdf(file_out, format=format, mode=mode)


if __name__ == "__main__":
    ###############################
    ### Process input arguments ###
    ###############################
    parser = argparse.ArgumentParser(
        description="Creates a surface dataset with all natural land moved into crops, and with every crop present in every gridcell."
    )

    # Required
    parser.add_argument(
        "-i",
        "--input-file",
        help="Surface dataset (fsurdat) file to process",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Directory in which to save the new surface dataset file",
        required=True,
    )

    # Get arguments
    args = parser.parse_args(sys.argv[1:])

    # Check arguments

    # Process
    main(args.input_file, args.output_dir)
