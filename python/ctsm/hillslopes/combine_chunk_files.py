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

class HillslopeVars:
    """
    Fields to be added to hillslope_file
    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, ncolumns_per_gridcell, nhillslope, sjm, sim, recurse=True):

        # Variables that will actually be saved
        self.h_elev = np.zeros((ncolumns_per_gridcell, sjm, sim))
        self.h_dist = np.zeros((ncolumns_per_gridcell, sjm, sim))
        self.h_area = np.zeros((ncolumns_per_gridcell, sjm, sim))
        self.h_slope = np.zeros((ncolumns_per_gridcell, sjm, sim))
        self.h_aspect = np.zeros((ncolumns_per_gridcell, sjm, sim))
        self.h_width = np.zeros((ncolumns_per_gridcell, sjm, sim))
        self.h_bedrock = np.zeros((ncolumns_per_gridcell, sjm, sim))
        self.h_stream_depth = np.zeros((sjm, sim))
        self.h_stream_width = np.zeros((sjm, sim))
        self.h_stream_slope = np.zeros((sjm, sim))
        self.nhillcolumns = np.zeros((sjm, sim), dtype=int)
        self.pct_hillslope = np.zeros((nhillslope, sjm, sim))
        self.hillslope_index = np.zeros((ncolumns_per_gridcell, sjm, sim), dtype=int)
        self.column_index = np.zeros((ncolumns_per_gridcell, sjm, sim), dtype=int)
        self.downhill_column_index = np.zeros((ncolumns_per_gridcell, sjm, sim), dtype=int)

        # Placeholders for read-in data from each chunk
        if recurse:
            self.chunk_mask = None
            self.this_chunk = HillslopeVars(ncolumns_per_gridcell, nhillslope, sjm, sim, recurse=False)


    def read(self, chunk_file, read_bedrock, read_stream):
        """Read hillslope variables from one chunk file

        Args:
            cfile (str): Path to chunk file
            read_bedrock (logical): Whether to read bedrock variable(s)
            read_stream (logical): Whether to read stream variable(s)
        """
        # pylint: disable=too-many-statements

        chunk_ds = Dataset(chunk_file, "r")
        self.chunk_mask = chunk_ds.variables["chunk_mask"][:]
        try:
            self.this_chunk.h_elev = chunk_ds.variables["hillslope_elevation"][:]
        except KeyError:
            self.this_chunk.h_elev = chunk_ds.variables["h_height"][:]
        try:
            self.this_chunk.h_dist = chunk_ds.variables["hillslope_distance"][:]
        except KeyError:
            self.this_chunk.h_dist = chunk_ds.variables["h_length"][:]
        try:
            self.this_chunk.h_width = chunk_ds.variables["hillslope_width"][:]
        except KeyError:
            self.this_chunk.h_width = chunk_ds.variables["h_width"][:]
        try:
            self.this_chunk.h_area = chunk_ds.variables["hillslope_area"][:]
        except KeyError:
            self.this_chunk.h_area = chunk_ds.variables["h_area"][:]
        try:
            self.this_chunk.h_slope = chunk_ds.variables["hillslope_slope"][:]
        except KeyError:
            self.this_chunk.h_slope = chunk_ds.variables["h_slope"][:]
        try:
            self.this_chunk.h_aspect = chunk_ds.variables["hillslope_aspect"][:]
        except KeyError:
            self.this_chunk.h_aspect = chunk_ds.variables["h_aspect"][:]
        if read_bedrock:
            try:
                self.this_chunk.h_bedrock = chunk_ds.variables["hillslope_bedrock_depth"][:]
            except KeyError:
                self.this_chunk.h_bedrock = chunk_ds.variables["h_bedrock"][:]
        if read_stream:
            try:
                self.this_chunk.h_stream_depth = chunk_ds.variables["h_stream_depth"][:]
            except KeyError:
                self.this_chunk.h_stream_depth = chunk_ds.variables["hillslope_stream_depth"][:]
            try:
                self.this_chunk.h_stream_width = chunk_ds.variables["hillslope_stream_width"][:]
            except KeyError:
                self.this_chunk.h_stream_width = chunk_ds.variables["h_stream_width"][:]
            try:
                self.this_chunk.h_stream_slope = chunk_ds.variables["hillslope_stream_slope"][:]
            except KeyError:
                self.h_stream_slope = chunk_ds.variables["hillslope_stream_slope"][:]

        self.this_chunk.nhillcolumns = chunk_ds.variables["nhillcolumns"][
                :,
            ].astype(int)
        self.this_chunk.pct_hillslope = chunk_ds.variables["pct_hillslope"][
                :,
            ]
        self.this_chunk.hillslope_index = chunk_ds.variables["hillslope_index"][
                :,
            ].astype(int)
        self.this_chunk.column_index = chunk_ds.variables["column_index"][
                :,
            ].astype(int)
        self.this_chunk.downhill_column_index = chunk_ds.variables["downhill_column_index"][
                :,
            ].astype(int)
        chunk_ds.close()

    def update(self, i, j, add_bedrock, add_stream, landmask):
        """
        Update a gridcell in chunk
        """
        if not(self.chunk_mask[j, i] > 0 and landmask[j, i] > 0):
            return

        self.h_elev[:, j, i] = self.this_chunk.h_elev[:, j, i]
        self.h_dist[:, j, i] = self.this_chunk.h_dist[:, j, i]
        self.h_width[:, j, i] = self.this_chunk.h_width[:, j, i]
        self.h_area[:, j, i] = self.this_chunk.h_area[:, j, i]
        self.h_slope[:, j, i] = self.this_chunk.h_slope[:, j, i]
        self.h_aspect[:, j, i] = self.this_chunk.h_aspect[:, j, i]
        if add_bedrock:
            self.h_bedrock[:, j, i] = self.this_chunk.h_bedrock[:, j, i]
        if add_stream:
            self.h_stream_depth[j, i] = self.this_chunk.h_stream_depth[j, i]
            self.h_stream_width[j, i] = self.this_chunk.h_stream_width[j, i]
            self.h_stream_slope[j, i] = self.this_chunk.h_stream_slope[j, i]

        self.nhillcolumns[j, i] = self.this_chunk.nhillcolumns[j, i]
        self.pct_hillslope[:, j, i] = self.this_chunk.pct_hillslope[:, j, i]
        self.hillslope_index[:, j, i] = self.this_chunk.hillslope_index[:, j, i]
        self.column_index[:, j, i] = self.this_chunk.column_index[:, j, i]
        self.downhill_column_index[:, j, i] = self.this_chunk.downhill_column_index[:, j, i]

        # if 2 or less valid aspects, remove all hillslope data
        if self.this_chunk.nhillcolumns[j, i] > 0:
            # check number of hillslopes
            h_ndx = self.this_chunk.hillslope_index[:, j, i]
            nactual_hillslopes = np.unique(h_ndx[h_ndx > 0]).size
            if nactual_hillslopes < 3:
                self.h_elev[:, j, i] = 0
                self.h_dist[:, j, i] = 0
                self.h_width[:, j, i] = 0
                self.h_area[:, j, i] = 0
                self.h_slope[:, j, i] = 0
                self.h_aspect[:, j, i] = 0
                self.h_bedrock[:, j, i] = 0
                self.h_stream_depth[j, i] = 0
                self.h_stream_width[j, i] = 0
                self.h_stream_slope[j, i] = 0

                self.nhillcolumns[j, i] = 0
                self.pct_hillslope[:, j, i] = 0
                self.hillslope_index[:, j, i] = 0
                self.column_index[:, j, i] = 0
                self.downhill_column_index[:, j, i] = 0

    def save(self, input_file, output_file, ncolumns_per_gridcell, nhillslope, add_bedrock, add_stream):
        """
        Save to netCDF
        """
        print("saving")

        with xr.open_dataset(input_file) as ds_in:
            lsmlat = ds_in["lsmlat"].load()
            lsmlon = ds_in["lsmlon"].load()

        ds_out = xr.Dataset()

        ds_out = add_variable_xr(
            name="hillslope_elevation",
            units="m",
            long_name="hillslope elevation above channel",
            data = self.h_elev,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

        ds_out = add_variable_xr(
            name="hillslope_distance",
            units="m",
            long_name="hillslope distance from channel",
            data = self.h_dist,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

        ds_out = add_variable_xr(
            name="hillslope_width",
            units="m",
            long_name="hillslope width",
            data = self.h_width,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

        ds_out = add_variable_xr(
            name="hillslope_area",
            units="m2",
            long_name="hillslope area",
            data = self.h_area,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

        ds_out = add_variable_xr(
            name="hillslope_slope",
            units="m/m",
            long_name="hillslope slope",
            data = self.h_slope,
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

        ds_out = add_variable_xr(
            name="hillslope_aspect",
            units="radians",
            long_name="hillslope aspect (clockwise from North)",
            data = self.h_aspect,
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
                data = self.h_bedrock,
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
                data = self.h_stream_depth,
                dataset=ds_out,
                lsmlat=lsmlat,
                lsmlon=lsmlon,
                dims=["lsmlat", "lsmlon"],
            )
            ds_out = add_variable_xr(
                name="hillslope_stream_width",
                units="meters",
                long_name="stream channel bankfull width",
                data = self.h_stream_width,
                dataset=ds_out,
                lsmlat=lsmlat,
                lsmlon=lsmlon,
                dims=["lsmlat", "lsmlon"],
            )
            ds_out = add_variable_xr(
                name="hillslope_stream_slope",
                units="m/m",
                long_name="stream channel slope",
                data = self.h_stream_slope,
                dataset=ds_out,
                lsmlat=lsmlat,
                lsmlon=lsmlon,
                dims=["lsmlat", "lsmlon"],
            )

        ds_out = add_variable_xr(
            name="nhillcolumns",
            units="unitless",
            long_name="number of columns per landunit",
            data = self.nhillcolumns.astype(np.int32),
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            dims=["lsmlat", "lsmlon"],
        )

        ds_out = add_variable_xr(
            name="pct_hillslope",
            units="per cent",
            long_name="percent hillslope of landunit",
            data = self.pct_hillslope,
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
            data = self.hillslope_index.astype(np.int32),
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

        ds_out = add_variable_xr(
            name="column_index",
            units="unitless",
            long_name="column index",
            data = self.column_index.astype(np.int32),
            dataset=ds_out,
            lsmlat=lsmlat,
            lsmlon=lsmlon,
            ncolumns_per_gridcell=ncolumns_per_gridcell,
        )

        ds_out = add_variable_xr(
            name="downhill_column_index",
            units="unitless",
            long_name="downhill column index",
            data = self.downhill_column_index.astype(np.int32),
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
        ds_out.to_netcdf(output_file, "w", format="NETCDF4_CLASSIC")

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
            sjm = len(chunk_ds.dimensions["lsmlat"])
            sim = len(chunk_ds.dimensions["lsmlon"])

            add_bedrock = "hillslope_bedrock_depth" in chunk_ds.variables.keys()
            add_stream = "hillslope_stream_depth" in chunk_ds.variables.keys()

            chunk_ds.close()

            hillslope_vars = HillslopeVars(ncolumns_per_gridcell, nhillslope, sjm, sim)
            arrays_uninitialized = False

        if not file_exists:
            if args.verbose:
                print(f"Skipping; chunk file not found: {chunk_file}")
            continue

        # Read hillslope variables from one chunk file
        hillslope_vars.read(chunk_file, add_bedrock, add_stream)

        for i in range(sim):
            for j in range(sjm):
                hillslope_vars.update(i, j, add_bedrock, add_stream, landmask)

    # -- Write data to file ------------------
    hillslope_vars.save(args.input_file, args.output_file, ncolumns_per_gridcell, nhillslope, add_bedrock, add_stream)

    print(args.output_file + " created")
