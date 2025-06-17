"""
Combine chunk files into a file for use in CTSM
"""

import argparse
import sys
import os
import datetime as dt
import numpy as np

# The below "pylint: disable" is because pylint complains that netCDF4 has no member Dataset, even
# though it does.
from netCDF4 import Dataset  # pylint: disable=no-name-in-module
from ctsm import ctsm_logging
from ctsm.hillslopes.hillslope_utils import (
    add_variable_nc,
    add_longxy_latixy_nc,
    NETCDF_FORMAT,
    get_chunks_to_process,
)
from ctsm.hillslopes.hillslope_vars import HillslopeVars


def parse_arguments(argv):
    """
    Parse arguments to script
    """
    parser = argparse.ArgumentParser(
        description="Combine files for each chunk into a file for use in CTSM"
    )

    # Use these groups to organize --help output. Otherwise, required but named (i.e., non-
    # positional) arguments (e.g., --input-file) get shown as optional.
    required_named = parser.add_argument_group("Required named arguments")
    optional_named = parser.add_argument_group("Optional named arguments")
    ctsm_logging.add_logging_args(parser)

    # Input and output file settings
    required_named.add_argument(
        "-i",
        "--input-file",
        help="Input surface dataset",
        required=True,
    )
    required_named.add_argument(
        "-d",
        "--input-dir",
        help="Directory containing combined-chunk files (outputs of combine_gridcell_files)",
        required=True,
    )
    required_named.add_argument(
        "-o",
        "--output-file",
        help="Output file",
        required=True,
    )
    optional_named.add_argument("--overwrite", help="overwrite", action="store_true", default=False)

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
    if os.path.exists(args.output_file) and not args.overwrite:
        msg = f"Output file already exists: {args.output_file}"
        ctsm_logging.logger.error(msg)
        raise FileExistsError(msg)

    return args


def finish_saving(args):
    """
    Save some extra stuff to the netCDF file
    """
    with Dataset(args.input_file, "r") as ds_in:
        lon2d = ds_in.variables["LONGXY"][:]
        lat2d = ds_in.variables["LATIXY"][:]
        has_area = "AREA" in ds_in.variables
        if has_area:
            area = ds_in.variables["AREA"][:]
    with Dataset(args.output_file, "a", format=NETCDF_FORMAT) as ds_out:
        add_longxy_latixy_nc(lon2d, lat2d, ds_out)
        if has_area:
            add_variable_nc(
                name="AREA",
                units="km^2",
                long_name="area",
                data=area,
                dataset=ds_out,
                dims=["lsmlat", "lsmlon"],
            )
        ds_out.setncattr("input_file", args.input_file)
        now = dt.datetime.now()
        datestr = now.strftime("%Y-%m-%d")
        ds_out.setncattr("creation_date", datestr)


def get_mask_var(surface_ds):
    """
    Get the variable to be used as a mask
    """
    mask_var = None
    mask_var_options = ["PFTDATA_MASK", "LANDFRAC_PFT"]
    for mask_var_option in mask_var_options:
        if mask_var_option in surface_ds.variables.keys():
            mask_var = mask_var_option
    if mask_var is None:
        msg = f"No variable found in sfcfile that looks like a mask ({mask_var_options})"
        ctsm_logging.logger.error(msg)
        raise KeyError(msg)

    landmask = np.asarray(surface_ds.variables[mask_var][:,])
    return mask_var, landmask


def main():
    """
    See module description
    """
    ctsm_logging.setup_logging_pre_config()
    args = parse_arguments(sys.argv[1:])

    # Choose data files to combine and append
    cfile0 = os.path.join(
        args.input_dir,
        f"combined_chunk_ChunkIndex_HAND_{args.n_bins}_col_hillslope_geo_params"
        + f"_{args.hillslope_form}_{args.dem_source}.nc",
    )

    surface_ds = Dataset(args.input_file, "r")
    mask_var, landmask = get_mask_var(surface_ds)
    surface_ds.close()
    if mask_var == "LANDFRAC_PFT":
        landmask[np.where(landmask > 0)] = 1

    hillslope_vars = None
    add_bedrock = None
    add_stream = None
    file_prefix = "combined_chunk"
    chunks_to_process = get_chunks_to_process(args, file_prefix)
    for cndx in chunks_to_process:
        ctsm_logging.logger.info("Chunk %d...", cndx)
        cstr = "{:02d}".format(cndx)
        chunk_file = cfile0.replace("ChunkIndex", cstr)
        file_exists = os.path.exists(chunk_file)

        if hillslope_vars is None and file_exists:
            chunk_ds = Dataset(chunk_file, "r")

            ncolumns_per_gridcell = len(chunk_ds.dimensions["nmaxhillcol"])
            nhillslope = len(chunk_ds.dimensions["nhillslope"])
            n_lat = len(chunk_ds.dimensions["lsmlat"])
            n_lon = len(chunk_ds.dimensions["lsmlon"])

            add_bedrock = "hillslope_bedrock_depth" in chunk_ds.variables.keys()
            add_stream = "hillslope_stream_depth" in chunk_ds.variables.keys()

            chunk_ds.close()

            hillslope_vars = HillslopeVars(ncolumns_per_gridcell, nhillslope, n_lat, n_lon)

        if not file_exists:
            ctsm_logging.logger.info("Skipping; chunk file not found: %s", chunk_file)
            continue

        # Read hillslope variables from one chunk file
        if hillslope_vars is None:
            msg = f"No chunk files found in '{args.input_dir}' starting with '{file_prefix}'"
            ctsm_logging.logger.error(msg)
            raise FileNotFoundError(msg)
        hillslope_vars.read(chunk_file, add_bedrock, add_stream)

        for i in range(n_lon):
            for j in range(n_lat):
                hillslope_vars.update(i, j, add_bedrock, add_stream, landmask=landmask)

    if hillslope_vars is None:
        msg = f"No files found in '{args.input_dir}'"
        ctsm_logging.logger.error(msg)
        raise FileNotFoundError(msg)

    # -- Write data to file ------------------
    hillslope_vars.save(
        input_file=args.input_file,
        output_file=args.output_file,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
        nhillslope=nhillslope,
        add_bedrock=add_bedrock,
        add_stream=add_stream,
        logger=ctsm_logging.logger,
    )
    finish_saving(args)

    ctsm_logging.logger.info("%s created", args.output_file)
