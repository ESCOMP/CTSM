"""
Combine chunk files into single surface data file
"""
import argparse
import sys
import os
import shutil
import datetime
import numpy as np
import xarray as xr

# The below "pylint: disable" is because pylint complains that netCDF4 has no member Dataset, even
# though it does.
from netCDF4 import Dataset  # pylint: disable=no-name-in-module
from ctsm.hillslopes.hillslope_utils import add_variable_xr



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
        nargs=1,
        type=int,
        default=default_n_chunks,
    )
    default_n_bins = 4
    parser.add_argument(
        "--n-bins",
        type=int,
        default=default_n_bins,
        help=f"Number of elevation bins (default: {default_n_bins})",
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


def main():
    """
    See module description
    """
    args = parse_arguments(sys.argv[1:])

    # Choose data files to combine and append
    cfile0 = os.path.join(
        args.input_dir,
        f"combined_chunk_ChunkIndex_HAND_{args.n_bins}_col_hillslope_geo_params"
        + f"_trapezoid_{args.dem_source}.nc",
    )

    f = Dataset(args.input_file, "r")
    landmask = f.variables["PFTDATA_MASK"][
        :,
    ]
    f.close()

    arrays_uninitialized = True
    for cndx in 1 + np.arange(args.n_chunks):
        print(f"{cndx} / {args.n_chunks}")
        cstr = "{:02d}".format(cndx)
        cfile = cfile0.replace("ChunkIndex", cstr)
        file_exists = os.path.exists(cfile)

        if arrays_uninitialized and file_exists:
            f = Dataset(cfile, "r")

            ncolumns_per_gridcell = len(f.dimensions["nmaxhillcol"])
            nhillslope = len(f.dimensions["nhillslope"])
            sjm = len(f.dimensions["lsmlat"])
            sim = len(f.dimensions["lsmlon"])

            # initialize new fields to be added to surface data file
            h_elev = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_dist = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_area = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_slope = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_aspect = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_width = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_pftndx = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_bedrock = np.zeros((ncolumns_per_gridcell, sjm, sim))
            h_stream_depth = np.zeros((sjm, sim))
            h_stream_width = np.zeros((sjm, sim))
            h_stream_slope = np.zeros((sjm, sim))

            nhillcolumns = np.zeros((sjm, sim), dtype=int)
            pct_hillslope = np.zeros((nhillslope, sjm, sim))
            hillslope_index = np.zeros((ncolumns_per_gridcell, sjm, sim), dtype=int)
            column_index = np.zeros((ncolumns_per_gridcell, sjm, sim), dtype=int)
            downhill_column_index = np.zeros((ncolumns_per_gridcell, sjm, sim), dtype=int)

            add_bedrock = "hillslope_bedrock_depth" in f.variables.keys()
            add_stream = "hillslope_stream_depth" in f.variables.keys()

            arrays_uninitialized = False
            f.close()

        if not file_exists:
            if args.verbose:
                print(f"Skipping; chunk file not found: {cfile}")
            continue

        f = Dataset(cfile, "r")
        nhillslope = len(f.dimensions["nhillslope"])
        chunk_mask = f.variables["chunk_mask"][:]
        try:
            h_elev0 = f.variables["hillslope_elevation"][:]
        except KeyError:
            h_elev0 = f.variables["h_height"][:]
        try:
            h_dist0 = f.variables["hillslope_distance"][:]
        except KeyError:
            h_dist0 = f.variables["h_length"][:]
        try:
            h_width0 = f.variables["hillslope_width"][:]
        except KeyError:
            h_width0 = f.variables["h_width"][:]
        try:
            h_area0 = f.variables["hillslope_area"][:]
        except KeyError:
            h_area0 = f.variables["h_area"][:]
        try:
            h_slope0 = f.variables["hillslope_slope"][:]
        except KeyError:
            h_slope0 = f.variables["h_slope"][:]
        try:
            h_aspect0 = f.variables["hillslope_aspect"][:]
        except KeyError:
            h_aspect0 = f.variables["h_aspect"][:]
        if add_bedrock:
            try:
                h_bedrock0 = f.variables["hillslope_bedrock_depth"][:]
            except KeyError:
                h_bedrock0 = f.variables["h_bedrock"][:]
        if add_stream:
            try:
                h_stream_depth0 = f.variables["h_stream_depth"][:]
            except KeyError:
                h_stream_depth0 = f.variables["hillslope_stream_depth"][:]
            try:
                h_stream_width0 = f.variables["hillslope_stream_width"][:]
            except KeyError:
                h_stream_width0 = f.variables["h_stream_width"][:]
            try:
                h_stream_slope0 = f.variables["hillslope_stream_slope"][:]
            except KeyError:
                h_stream_slope = f.variables["hillslope_stream_slope"][:]

        nhillcolumns0 = f.variables["nhillcolumns"][
            :,
        ].astype(int)
        pct_hillslope0 = f.variables["pct_hillslope"][
            :,
        ]
        hillslope_index0 = f.variables["hillslope_index"][
            :,
        ].astype(int)
        column_index0 = f.variables["column_index"][
            :,
        ].astype(int)
        downhill_column_index0 = f.variables["downhill_column_index"][
            :,
        ].astype(int)
        f.close()

        for i in range(sim):
            for j in range(sjm):
                if not(chunk_mask[j, i] > 0 and landmask[j, i] > 0):
                    continue

                h_elev[:, j, i] = h_elev0[:, j, i]
                h_dist[:, j, i] = h_dist0[:, j, i]
                h_width[:, j, i] = h_width0[:, j, i]
                h_area[:, j, i] = h_area0[:, j, i]
                h_slope[:, j, i] = h_slope0[:, j, i]
                h_aspect[:, j, i] = h_aspect0[:, j, i]
                if add_bedrock:
                    h_bedrock[:, j, i] = h_bedrock0[:, j, i]
                if add_stream:
                    h_stream_depth[j, i] = h_stream_depth0[j, i]
                    h_stream_width[j, i] = h_stream_width0[j, i]
                    h_stream_slope[j, i] = h_stream_slope0[j, i]

                nhillcolumns[j, i] = nhillcolumns0[j, i]
                pct_hillslope[:, j, i] = pct_hillslope0[:, j, i]
                hillslope_index[:, j, i] = hillslope_index0[:, j, i]
                column_index[:, j, i] = column_index0[:, j, i]
                downhill_column_index[:, j, i] = downhill_column_index0[:, j, i]

                # if 2 or less valid aspects, remove all hillslope data
                if nhillcolumns0[j, i] > 0:
                    # check number of hillslopes
                    h_ndx = hillslope_index0[:, j, i]
                    nactual_hillslopes = np.unique(h_ndx[h_ndx > 0]).size
                    if nactual_hillslopes < 3:
                        h_elev[:, j, i] = 0
                        h_dist[:, j, i] = 0
                        h_width[:, j, i] = 0
                        h_area[:, j, i] = 0
                        h_slope[:, j, i] = 0
                        h_aspect[:, j, i] = 0
                        h_bedrock[:, j, i] = 0
                        h_stream_depth[j, i] = 0
                        h_stream_width[j, i] = 0
                        h_stream_slope[j, i] = 0

                        nhillcolumns[j, i] = 0
                        pct_hillslope[:, j, i] = 0
                        hillslope_index[:, j, i] = 0
                        column_index[:, j, i] = 0
                        downhill_column_index[:, j, i] = 0

    # -- Write data to file ------------------

    print("saving")

    with xr.open_dataset(args.input_file) as ds_in:
        lsmlat = ds_in["lsmlat"].load()
        lsmlon = ds_in["lsmlon"].load()

    ds_out = xr.Dataset()

    ds_out = add_variable_xr(
        name="hillslope_elevation",
        units="m",
        long_name="hillslope elevation above channel",
        data=h_elev,
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )

    ds_out = add_variable_xr(
        name="hillslope_distance",
        units="m",
        long_name="hillslope distance from channel",
        data=h_dist,
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )

    ds_out = add_variable_xr(
        name="hillslope_width",
        units="m",
        long_name="hillslope width",
        data=h_width,
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )

    ds_out = add_variable_xr(
        name="hillslope_area",
        units="m2",
        long_name="hillslope area",
        data=h_area,
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )

    ds_out = add_variable_xr(
        name="hillslope_slope",
        units="m/m",
        long_name="hillslope slope",
        data=h_slope,
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )

    ds_out = add_variable_xr(
        name="hillslope_aspect",
        units="radians",
        long_name="hillslope aspect (clockwise from North)",
        data=h_aspect,
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )


    if add_bedrock:
        ds_out = add_variable_xr(
            name="hillslope_bedrock_depth",
            units="meters",
            long_name="hillslope bedrock depth",
            data=h_bedrock,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

    if add_stream:
        ds_out = add_variable_xr(
            name="hillslope_stream_depth",
            units="meters",
            long_name="stream channel bankfull depth",
            data=h_stream_depth,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            dims=["lsmlat", "lsmlon"],
        )
        ds_out = add_variable_xr(
            name="hillslope_stream_width",
            units="meters",
            long_name="stream channel bankfull width",
            data=h_stream_width,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            dims=["lsmlat", "lsmlon"],
        )
        ds_out = add_variable_xr(
            name="hillslope_stream_slope",
            units="m/m",
            long_name="stream channel slope",
            data=h_stream_slope,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            dims=["lsmlat", "lsmlon"],
        )

    ds_out = add_variable_xr(
        name="nhillcolumns",
        units="unitless",
        long_name="number of columns per landunit",
        data=nhillcolumns.astype(np.int32),
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        dims=["lsmlat", "lsmlon"],
    )

    ds_out = add_variable_xr(
        name="pct_hillslope",
        units="per cent",
        long_name="percent hillslope of landunit",
        data=pct_hillslope,
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        dims=["nhillslope", "lsmlat", "lsmlon"],
        nhillslope=nhillslope,
    )

    ds_out = add_variable_xr(
        name="hillslope_index",
        units="unitless",
        long_name="hillslope_index",
        data=hillslope_index.astype(np.int32),
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )

    ds_out = add_variable_xr(
        name="column_index",
        units="unitless",
        long_name="column index",
        data=column_index.astype(np.int32),
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )

    ds_out = add_variable_xr(
        name="downhill_column_index",
        units="unitless",
        long_name="downhill column index",
        data=downhill_column_index.astype(np.int32),
        dataset=ds_out,
        lsmlat=lsmlat,
        lsmlon=lsmlon,
        ncolumns_per_gridcell=ncolumns_per_gridcell,
    )


    # Before saving, drop coordinate variables (which aren't present on fsurdat)
    ds_out = ds_out.drop_vars(
        [
            "lsmlat",
            "lsmlon",
            "nmaxhillcol",
            "nhillslope",
        ]
    )

    # Save
    ds_out.to_netcdf(args.output_file, "w")

    print(args.output_file + " created")
