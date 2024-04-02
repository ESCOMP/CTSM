"""
Utilities for hillslope scripts to share
"""
import numpy as np


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
        )
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
        outfile,
        "hillslope_stream_depth",
        "m",
        "stream channel bankfull depth",
        dims=dims
    )
    wwidth = create_variable(
        outfile,
        "hillslope_stream_width",
        "m",
        "stream channel bankfull width",
        dims=dims
    )
    wslope = create_variable(
        outfile,
        "hillslope_stream_slope",
        "m/m",
        "stream channel slope",
        dims=dims
    )

    return wdepth, wwidth, wslope


def create_variable(
        outfile, name, units, long_name,
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
