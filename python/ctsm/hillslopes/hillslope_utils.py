"""
Utilities for hillslope scripts to share
"""

import os
import glob
import re
import numpy as np

NETCDF_FORMAT = "NETCDF3_64BIT_OFFSET"

def add_longxy_latixy_nc(lon2d, lat2d, ds_out):
    """
    Add LONGXY and LATIXY to a netCDF file
    """
    add_variable_nc(
        name="LONGXY",
        units="degrees east",
        long_name="longitude",
        data=lon2d,
        dataset=ds_out,
        dims=["lsmlat", "lsmlon"],
    )
    add_variable_nc(
        name="LATIXY",
        units="degrees north",
        long_name="latitude",
        data=lat2d,
        dataset=ds_out,
        dims=["lsmlat", "lsmlon"],
    )


def create_variables(outfile):
    """
    Create variables
    """
    ohand = create_variable(
        outfile,
        "hillslope_elevation",
        "m",
        "hillslope elevation",
    )

    odtnd = create_variable(
        outfile,
        "hillslope_distance",
        "m",
        "hillslope distance",
    )

    owidth = create_variable(
        outfile,
        "hillslope_width",
        "m",
        "hillslope width",
    )

    oarea = create_variable(
        outfile,
        "hillslope_area",
        "m2",
        "hillslope area",
    )

    oslop = create_variable(
        outfile,
        "hillslope_slope",
        "m/m",
        "hillslope slope",
    )

    oasp = create_variable(
        outfile,
        "hillslope_aspect",
        "radians",
        "hillslope aspect (clockwise from North)",
    )

    onhill = create_variable(
        outfile,
        "nhillcolumns",
        "unitless",
        "number of columns per landunit",
        dims=(
            "lsmlat",
            "lsmlon",
        ),
        data_type=np.int32,
    )

    opcthill = create_variable(
        outfile,
        "pct_hillslope",
        "per cent",
        "percent hillslope of landunit",
        dims=(
            "nhillslope",
            "lsmlat",
            "lsmlon",
        ),
    )

    ohillndx = create_variable(
        outfile,
        "hillslope_index",
        "unitless",
        "hillslope_index",
        data_type=np.int32,
    )

    ocolndx = create_variable(
        outfile,
        "column_index",
        "unitless",
        "column index",
        data_type=np.int32,
    )

    odcolndx = create_variable(
        outfile,
        "downhill_column_index",
        "unitless",
        "downhill column index",
        data_type=np.int32,
    )

    obed = create_variable(
        outfile,
        "hillslope_bedrock_depth",
        "meters",
        "hillslope bedrock depth",
    )

    return (
        ohand,
        odtnd,
        owidth,
        oarea,
        oslop,
        oasp,
        onhill,
        opcthill,
        ohillndx,
        ocolndx,
        odcolndx,
        obed,
    )


def add_stream_channel_vars(outfile):
    """
    Add stream channel variables
    """
    dims = ("lsmlat", "lsmlon")

    wdepth = create_variable(
        outfile, "hillslope_stream_depth", "m", "stream channel bankfull depth", dims=dims
    )
    wwidth = create_variable(
        outfile, "hillslope_stream_width", "m", "stream channel bankfull width", dims=dims
    )
    wslope = create_variable(
        outfile, "hillslope_stream_slope", "m/m", "stream channel slope", dims=dims
    )

    return wdepth, wwidth, wslope


def create_variable(
    outfile,
    name,
    units,
    long_name,
    *,
    dims=(
        "nmaxhillcol",
        "lsmlat",
        "lsmlon",
    ),
    data_type=np.float64,
):
    """
    Convenient function to use for making hillslope variables
    """

    # Variable in netCDF output file
    nc_var = outfile.createVariable(
        name,
        data_type,
        dims,
    )
    nc_var.units = units
    nc_var.long_name = long_name

    return nc_var


def add_variable_nc(
    *,
    name,
    units,
    long_name,
    data,
    dataset,
    data_type=np.float64,
    dims=("nmaxhillcol", "lsmlat", "lsmlon"),
):
    """
    Convenient function to use for adding hillslope variables to a Dataset with netcdf
    """

    if isinstance(dims, list):
        dims = tuple(dims)

    # Make variable
    nc_var = create_variable(
        dataset,
        name,
        units,
        long_name,
        dims=dims,
        data_type=data_type,
    )

    # Fill with data
    nc_var[:] = data


def get_chunks_to_process(args, prefix):
    """
    Get list of chunks to process
    """
    if not hasattr(args, "cndx") or args.cndx is None:
        # List of gridcell files
        file_list = glob.glob(os.path.join(args.input_dir, prefix + "_[0-9]*nc"))
        if not file_list:
            raise FileNotFoundError(f"No files found in '{args.input_dir}'")
        # Extract the chunk number from the file names
        chunk_list = [re.search(r"chunk_\d+", x).group() for x in file_list]
        chunk_list = [x.replace("chunk_", "") for x in chunk_list]
        # Get the list of unique chunk numbers
        chunks_to_process = [int(x) for x in list(set(chunk_list))]
        chunks_to_process.sort()
        if not chunks_to_process:
            raise FileNotFoundError(f"No MATCHING chunk files found in '{args.input_dir}'")
    else:
        chunks_to_process = [int(cndx) for cndx in args.cndx[0].split(",")]
        for cndx in chunks_to_process:
            if cndx < 1 or cndx > args.n_chunks:
                raise RuntimeError("All cndx must be 1-{:d}".format(args.n_chunks))
    return chunks_to_process
