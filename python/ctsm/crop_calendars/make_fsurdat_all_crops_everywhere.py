
import numpy as np
import xarray as xr
import sys
import argparse


def main(file_in, file_out):
    # Import
    
    file_in = "/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/release-clm5.0.18/surfdata_10x15_hist_78pfts_CMIP6_simyr2000_c190214.nc"
    in_ds = xr.open_dataset(file_in)
    
    out_ds = in_ds.copy()
    in_ds.close()
    
    # Move all natural land into crop
    out_ds['PCT_CROP'] += in_ds['PCT_NATVEG']
    out_ds['PCT_NATVEG'] -= in_ds['PCT_NATVEG']
    
    # Put some of every crop in every gridcell
    pct_cft =  np.full_like(in_ds['PCT_CFT'].values, 100 / in_ds.dims['cft'])
    out_ds['PCT_CFT'] = xr.DataArray(data=pct_cft,
                                     attrs=in_ds['PCT_CFT'].attrs,
                                     dims=in_ds['PCT_CFT'].dims)
    
    # Save
    out_ds.to_netcdf(file_out, format="NETCDF3_64BIT")


if __name__ == "__main__":
    ###############################
    ### Process input arguments ###
    ###############################
    parser = argparse.ArgumentParser(description="Creates a surface dataset with all natural land moved into crops, and with every crop present in every gridcell.")

    # Required
    parser.add_argument(
        "-i",
        "--input-file",
        help="Surface dataset (fsurdat) file to process",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Where to save the new surface dataset file",
        required=True,
    )

    # Get arguments
    args = parser.parse_args(sys.argv[1:])

    # Process
    main(args.input_file, args.output_file)

