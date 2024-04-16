"""
Combine gridcell files into single file
"""
import sys
import os
import argparse
import glob
import datetime
import numpy as np

# The below "pylint: disable" is because pylint complains that netCDF4 has no
# member Dataset, even though it does.
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm.hillslopes.hillslope_utils import HillslopeVars


def parse_arguments(argv):
    """
    Parse arguments to script
    """
    parser = argparse.ArgumentParser(description="Combine gridcell files into single file")

    parser.add_argument(
        "-i",
        "--input-file",
        help="Input surface dataset with grid information",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--input-dir",
        help="Directory containing chunk files",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Directory where output file should be saved (default: current dir)",
        default=os.getcwd(),
    )
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
        nargs=1,
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
    parser.add_argument(
        "--cndx",
        help=(
            "Chunk(s) to process. If excluded will process all chunks (see --n-chunks). To "
            + "include, specify either a single chunk or a comma-separated list of chunks."
        ),
        nargs=1,
        type=str,
        default=None,
    )
    parser.add_argument("--overwrite", help="overwrite", action="store_true", default=False)
    parser.add_argument("-v", "--verbose", help="print info", action="store_true", default=False)

    args = parser.parse_args(argv)

    # Check arguments
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")
    if not os.path.exists(args.input_dir):
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    return args


def main():
    """
    See module description
    """

    args = parse_arguments(sys.argv[1:])
    verbose = args.verbose

    if args.cndx is None:
        chunks_to_process = 1 + np.arange(args.n_chunks)
    else:
        chunks_to_process = [int(cndx) for cndx in args.cndx[0].split(",")]
        for cndx in chunks_to_process:
            if cndx < 1 or cndx > args.n_chunks:
                raise RuntimeError("All cndx must be 1-{:d}".format(args.n_chunks))

    nhillslope = None
    nmaxhillcol = None

    for cndx in chunks_to_process:

        # Gridcell file directory
        cfile = os.path.join(
            args.input_dir,
            "chunk_{:02d}_HAND_{:d}_col_hillslope_geo_params_{}_{}.nc".format(
                cndx, args.n_bins, args.hillslope_form, args.dem_source
            ),
        )

        # Output file
        outfile_path = os.path.join(
            args.output_dir, os.path.split(cfile)[-1].replace("chunk_", "combined_chunk_")
        )

        # Read output file coordinates
        fsurdat = Dataset(args.input_file, "r")
        n_lat = len(fsurdat.dimensions["lsmlat"])
        n_lon = len(fsurdat.dimensions["lsmlon"])
        fsurdat.close()

        # Check for output file existence
        if os.path.exists(outfile_path):
            if args.overwrite:
                if verbose:
                    print(outfile_path, " exists; overwriting")
            else:
                print(outfile_path, " exists; skipping")
                continue

        # Locate gridcell files
        gfile = cfile.replace(".nc", "*.nc")
        gfiles = glob.glob(gfile)
        gfiles.sort()
        if len(gfiles) == 0:
            print(f"Chunk {cndx}: Skipping; no files found matching {gfile}")
            continue
        print(f"Chunk {cndx}: Combining {len(gfiles)} files...")

        # Read hillslope data dimensions/settings, if not done yet
        if nhillslope is None:
            with Dataset(gfiles[0], "r") as first_gridcell_file:
                nhillslope = len(first_gridcell_file.dimensions["nhillslope"])
                nmaxhillcol = len(first_gridcell_file.dimensions["nmaxhillcol"])
                add_bedrock = "hillslope_bedrock_depth" in first_gridcell_file.variables.keys()
                do_add_stream_channel_vars = (
                    "hillslope_stream_depth" in first_gridcell_file.variables.keys()
                )

        write_to_file(
            outfile_path,
            n_lat,
            n_lon,
            gfiles,
            nhillslope,
            nmaxhillcol,
            add_bedrock,
            do_add_stream_channel_vars,
        )


def write_to_file(
    outfile_path,
    n_lat,
    n_lon,
    gfiles,
    nhillslope,
    nmaxhillcol,
    add_bedrock,
    do_add_stream_channel_vars,
):
    """
    Write to file
    """
    # pylint: disable=too-many-statements
    outfile = Dataset(outfile_path, "w")
    outfile.creation_date = datetime.date.today().strftime("%y%m%d")

    outfile.createDimension("lsmlon", n_lon)
    outfile.createDimension("lsmlat", n_lat)
    outfile.createDimension("nhillslope", nhillslope)
    outfile.createDimension("nmaxhillcol", nmaxhillcol)
    outfile.close()

    hillslope_vars = HillslopeVars(nmaxhillcol, nhillslope, n_lat, n_lon, incl_latlon=True)

    # loop over gridcell files
    for gfile in gfiles:
        y1, x1 = gfile.index("j_"), gfile.index("i_")  # pylint: disable=invalid-name
        j, i = int(gfile[y1 + 2 : y1 + 5]), int(gfile[x1 + 2 : x1 + 5])

        hillslope_vars.read(gfile, add_bedrock, do_add_stream_channel_vars, incl_latlon=True)

        hillslope_vars.update(
            i,
            j,
            add_bedrock,
            do_add_stream_channel_vars,
            incl_latlon=True,
            incl_chunkmask=True,
            this_chunk_1d=True,
            remove_if_too_few_aspects=False,
        )

    hillslope_vars.save(
        None,
        outfile_path,
        nmaxhillcol,
        nhillslope,
        add_bedrock,
        do_add_stream_channel_vars,
        n_lon=n_lon,
        n_lat=n_lat,
        incl_latlon=True,
        incl_chunkmask=True,
    )
    print(outfile_path + " created")
