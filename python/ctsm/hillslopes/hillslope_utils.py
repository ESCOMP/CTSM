"""
Utilities for hillslope scripts to share
"""
import numpy as np
import xarray as xr


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

def add_variable_nc(
        name, units, long_name,
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
        dataset, name, units, long_name,
        dims=dims,
        data_type=data_type,
        )

    # Fill with data
    nc_var[:] = data

def add_variable_xr(
        name, units, long_name,
        data,
        dataset,
        lsmlat, lsmlon,
        dims=["nmaxhillcol", "lsmlat", "lsmlon"],
        ncolumns_per_gridcell=None,
        nhillslope=None,
        ):  # pylint: disable=dangerous-default-value
    # pylint disable above: pylint thinks the dims default is empty list []
    """
    Convenient function to use for adding hillslope variables to a Dataset with xarray
    """

    if isinstance(dims, tuple):
        dims = list(dims)

    coords = {}
    if dims == ["nmaxhillcol", "lsmlat", "lsmlon"]:
        if ncolumns_per_gridcell is None:
            raise RuntimeError("If nmaxhillcol is in dims, provide ncolumns_per_gridcell argument")
        coords = {
            "nmaxhillcol": np.arange(ncolumns_per_gridcell),
        }
    elif dims == ["nhillslope", "lsmlat", "lsmlon"]:
        if nhillslope is None:
            raise RuntimeError("If nhillslope is in dims, provide nhillslope argument")
        coords = {
            "nhillslope": np.arange(nhillslope),
        }
    elif dims != ["lsmlat", "lsmlon"]:
        raise RuntimeError(f"Unhandled dim list: {dims}")
    coords["lsmlat"] = lsmlat
    coords["lsmlon"] = lsmlon

    data_array = xr.DataArray(
        data=data,
        dims=dims,
        coords=coords,
        attrs={
            "units": units,
            "long_name": long_name,
        }
    )

    dataset[name] = data_array

    return dataset
