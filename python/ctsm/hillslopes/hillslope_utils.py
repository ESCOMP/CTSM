"""
Utilities for hillslope scripts to share
"""
import numpy as np


def create_variables(outfile):
    """
    Create variables
    """
    ohand = outfile.createVariable(
        "hillslope_elevation",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    ohand.units = "m"
    ohand.long_name = "hillslope elevation"

    odtnd = outfile.createVariable(
        "hillslope_distance",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    odtnd.units = "m"
    odtnd.long_name = "hillslope distance"

    owidth = outfile.createVariable(
        "hillslope_width",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    owidth.units = "m"
    owidth.long_name = "hillslope width"

    oarea = outfile.createVariable(
        "hillslope_area",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    oarea.units = "m2"
    oarea.long_name = "hillslope area"

    oslop = outfile.createVariable(
        "hillslope_slope",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    oslop.units = "m/m"
    oslop.long_name = "hillslope slope"

    oasp = outfile.createVariable(
        "hillslope_aspect",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    oasp.units = "radians"
    oasp.long_name = "hillslope aspect (clockwise from North)"

    onhill = outfile.createVariable(
        "nhillcolumns",
        np.int32,
        (
            "lsmlat",
            "lsmlon",
        ),
    )
    onhill.units = "unitless"
    onhill.long_name = "number of columns per landunit"

    opcthill = outfile.createVariable(
        "pct_hillslope",
        np.float64,
        (
            "nhillslope",
            "lsmlat",
            "lsmlon",
        ),
    )
    opcthill.units = "per cent"
    opcthill.long_name = "percent hillslope of landunit"

    ohillndx = outfile.createVariable(
        "hillslope_index",
        np.int32,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    ohillndx.units = "unitless"
    ohillndx.long_name = "hillslope_index"

    ocolndx = outfile.createVariable(
        "column_index",
        np.int32,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    ocolndx.units = "unitless"
    ocolndx.long_name = "column index"

    odcolndx = outfile.createVariable(
        "downhill_column_index",
        np.int32,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    odcolndx.units = "unitless"
    odcolndx.long_name = "downhill column index"

    obed = outfile.createVariable(
        "hillslope_bedrock_depth",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    obed.units = "meters"
    obed.long_name = "hillslope bedrock depth"

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
