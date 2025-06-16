"""
For each chunk in HAND analysis, combine all gridcell files into file.
"""

import sys
import os
import argparse
import glob
import datetime

# The below "pylint: disable" is because pylint complains that netCDF4 has no
# member Dataset, even though it does.
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm import ctsm_logging
from ctsm.hillslopes.hillslope_utils import HillslopeVars, get_chunks_to_process


def parse_arguments(argv):
    """
    Parse arguments to script
    """
    parser = argparse.ArgumentParser(description="Combine gridcell files into single file")

    # Use these groups to organize --help output. Otherwise, required but named (i.e., non-
    # positional) arguments (e.g., --input-file) get shown as optional.
    required_named = parser.add_argument_group("Required named arguments")
    optional_named = parser.add_argument_group("Optional named arguments")

    required_named.add_argument(
        "-i",
        "--input-file",
        help="Input surface dataset with grid information",
        required=True,
    )
    required_named.add_argument(
        "-d",
        "--input-dir",
        help="Directory containing chunk files",
        required=True,
    )
    optional_named.add_argument(
        "-o",
        "--output-dir",
        help="Directory where output file should be saved (default: current dir)",
        default=os.getcwd(),
    )
    dem_source_default = "MERIT"
    optional_named.add_argument(
        "--dem-source",
        help=f"DEM to use (default: {dem_source_default})",
        type=str,
        default=dem_source_default,
    )
    default_n_bins = 4
    optional_named.add_argument(
        "--n-bins",
        type=int,
        default=default_n_bins,
        help=f"Number of elevation bins (default: {default_n_bins}). "
        + "Used to generate input filename template.",
    )
    default_hillslope_form = "Trapezoidal"
    optional_named.add_argument(
        "--hillslope-form",
        help=f"Hillslope form (default: {default_hillslope_form}). "
        + "Used to generate input filename template.",
        type=str,
        default=default_hillslope_form,
    )
    optional_named.add_argument(
        "--cndx",
        help=(
            "Chunk(s) to process. If excluded will process all chunks (see --n-chunks). To "
            + "include, specify either a single chunk or a comma-separated list of chunks."
        ),
        nargs=1,
        type=str,
        default=None,
    )
    optional_named.add_argument(
        "--overwrite",
        help="Overwrite existing output files? If not given, will skip existing files.",
        action="store_true",
        default=False,
    )

    ctsm_logging.add_logging_args(parser)

    args = parser.parse_args(argv)
    ctsm_logging.process_logging_args(args)

    # Check arguments
    if not os.path.exists(args.input_file):
        msg = f"Input file not found: {args.input_file}"
        ctsm_logging.logger.error(msg)
        raise FileNotFoundError(msg)
    if not os.path.exists(args.input_dir):
        msg = f"Input directory not found: {args.input_dir}"
        ctsm_logging.logger.error(msg)
        raise FileNotFoundError(msg)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    return args


def main():
    """
    See module description
    """

    ctsm_logging.setup_logging_pre_config()
    args = parse_arguments(sys.argv[1:])

    chunks_to_process = get_chunks_to_process(args, "chunk")

    nhillslope = None
    add_bedrock = None
    do_add_stream_channel_vars = None
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
                ctsm_logging.logger.warning("%s exists; overwriting", outfile_path)
            else:
                ctsm_logging.logger.warning("%s exists; skipping", outfile_path)
                continue

        # Locate gridcell files
        gfile = cfile.replace(".nc", "*.nc")
        gfiles = glob.glob(gfile)
        gfiles.sort()
        if len(gfiles) == 0:
            ctsm_logging.logger.info("Chunk %d: Skipping; no files found matching %s", cndx, gfile)
            continue
        ctsm_logging.logger.info("Chunk %d: Combining %d files...", cndx, len(gfiles))

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
            outfile_path=outfile_path,
            n_lat=n_lat,
            n_lon=n_lon,
            gfiles=gfiles,
            nhillslope=nhillslope,
            nmaxhillcol=nmaxhillcol,
            add_bedrock=add_bedrock,
            do_add_stream_channel_vars=do_add_stream_channel_vars,
        )


def write_to_file(
    *,
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
        input_file=None,
        output_file=outfile_path,
        ncolumns_per_gridcell=nmaxhillcol,
        nhillslope=nhillslope,
        add_bedrock=add_bedrock,
        add_stream=do_add_stream_channel_vars,
        logger=ctsm_logging.logger,
        n_lon=n_lon,
        n_lat=n_lat,
        incl_latlon=True,
        incl_chunkmask=True,
    )
    ctsm_logging.logger.info("%s created", outfile_path)
