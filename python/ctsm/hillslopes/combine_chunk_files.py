"""
Combine chunk files into single surface data file
"""
import argparse
import sys
import os
import datetime as dt
import numpy as np

# The below "pylint: disable" is because pylint complains that netCDF4 has no member Dataset, even
# though it does.
from netCDF4 import Dataset  # pylint: disable=no-name-in-module
from ctsm.hillslopes.hillslope_utils import add_variable_nc, add_longxy_latixy_nc, HillslopeVars, NETCDF_FORMAT


def parse_arguments(argv):
    """
    Parse arguments to script
    """
    parser = argparse.ArgumentParser(description="Specify a synthetic hillslope profile")

    # Input and output file settings
    parser.add_argument(
        "-i",
        "--input-file",
        help="Input surface dataset",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--input-dir",
        help="Directory containing combined-chunk files (outputs of combine_gridcell_files)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Output file",
        required=True,
    )
    parser.add_argument("--overwrite", help="overwrite", action="store_true", default=False)

    dem_source_default = "MERIT"
    parser.add_argument(
        "--dem-source",
        help=f"DEM to use (default: {dem_source_default})",
        type=str,
        default=dem_source_default,
    )

    default_n_chunks = 36
    parser.add_argument(
        "--n-chunks",
        help=f"Number of chunks (default: {default_n_chunks})",
        type=int,
        default=default_n_chunks,
    )
    default_n_bins = 4
    parser.add_argument(
        "--n-bins",
        type=int,
        default=default_n_bins,
        help=f"Number of elevation bins (default: {default_n_bins}). "
             + "Used to generate input filename template.",
    )

    default_hillslope_form = "Trapezoidal"
    parser.add_argument(
        "--hillslope-form",
        help=f"Hillslope form (default: {default_hillslope_form}). "
             + "Used to generate input filename template.",
        type=str,
        default=default_hillslope_form,
    )

    parser.add_argument("-v", "--verbose", help="print info", action="store_true", default=False)

    args = parser.parse_args(argv)

    # Check arguments
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")
    if not os.path.exists(args.input_dir):
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")
    if os.path.exists(args.output_file) and not args.overwrite:
        raise FileExistsError(f"Output file already exists: {args.output_file}")

    return args


def finish_saving(args):
    """
    Save some extra stuff to the netCDF file
    """
    with Dataset(args.input_file, "r") as ds_in:
        lon2d = ds_in.variables["LONGXY"][:]
        lat2d = ds_in.variables["LATIXY"][:]
        area = ds_in.variables["AREA"][:]
    with Dataset(args.output_file, "a", format=NETCDF_FORMAT) as ds_out:
        add_longxy_latixy_nc(lon2d, lat2d, ds_out)
        add_variable_nc(
            name="AREA",
            units="km^2",
            long_name="area",
            data = area,
            dataset=ds_out,
            dims=["lsmlat", "lsmlon"],
        )
        ds_out.setncattr("input_file", args.input_file)
        now = dt.datetime.now()
        datestr = now.strftime("%Y-%m-%d")
        ds_out.setncattr("creation_date", datestr)


def main():
    """
    See module description
    """
    args = parse_arguments(sys.argv[1:])

    # Choose data files to combine and append
    cfile0 = os.path.join(
        args.input_dir,
        f"combined_chunk_ChunkIndex_HAND_{args.n_bins}_col_hillslope_geo_params"
        + f"_{args.hillslope_form}_{args.dem_source}.nc",
    )

    surface_ds = Dataset(args.input_file, "r")
    landmask = surface_ds.variables["PFTDATA_MASK"][
        :,
    ]
    surface_ds.close()

    arrays_uninitialized = True
    for cndx in 1 + np.arange(args.n_chunks):
        print(f"{cndx} / {args.n_chunks}")
        cstr = "{:02d}".format(cndx)
        chunk_file = cfile0.replace("ChunkIndex", cstr)
        file_exists = os.path.exists(chunk_file)

        if arrays_uninitialized and file_exists:
            chunk_ds = Dataset(chunk_file, "r")

            ncolumns_per_gridcell = len(chunk_ds.dimensions["nmaxhillcol"])
            nhillslope = len(chunk_ds.dimensions["nhillslope"])
            n_lat = len(chunk_ds.dimensions["lsmlat"])
            n_lon = len(chunk_ds.dimensions["lsmlon"])

            add_bedrock = "hillslope_bedrock_depth" in chunk_ds.variables.keys()
            add_stream = "hillslope_stream_depth" in chunk_ds.variables.keys()

            chunk_ds.close()

            hillslope_vars = HillslopeVars(ncolumns_per_gridcell, nhillslope, n_lat, n_lon)
            arrays_uninitialized = False

        if not file_exists:
            if args.verbose:
                print(f"Skipping; chunk file not found: {chunk_file}")
            continue

        # Read hillslope variables from one chunk file
        hillslope_vars.read(chunk_file, add_bedrock, add_stream)

        for i in range(n_lon):
            for j in range(n_lat):
                hillslope_vars.update(i, j, add_bedrock, add_stream, landmask=landmask)

    # -- Write data to file ------------------
    hillslope_vars.save(args.input_file, args.output_file, ncolumns_per_gridcell, nhillslope, add_bedrock, add_stream)
    finish_saving(args)

    print(args.output_file + " created")
