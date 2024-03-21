"""
Combine gridcell files into single file
"""
import sys
import os
import argparse
import glob
import datetime
import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module
from ctsm.hillslopes.hillslope_utils import create_variables as shared_create_variables

# The above "pylint: disable" is because pylint complains that netCDF4 has no
# member Dataset, even though it does.


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

    for cndx in chunks_to_process:

        # Gridcell file directory
        cfile = os.path.join(
            args.input_dir,
            "chunk_{:02d}_HAND_4_col_hillslope_geo_params_section_quad_{}.nc".format(
                cndx, args.dem_source
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

        # Read hillslope data dimensions
        first_gridcell_file = Dataset(gfiles[0], "r")
        nhillslope = len(first_gridcell_file.dimensions["nhillslope"])
        nmaxhillcol = len(first_gridcell_file.dimensions["nmaxhillcol"])

        add_bedrock = "hillslope_bedrock_depth" in first_gridcell_file.variables.keys()
        add_stream_channel_vars = "hillslope_stream_depth" in first_gridcell_file.variables.keys()

        first_gridcell_file.close()

        write_to_file(
            outfile_path,
            n_lat,
            n_lon,
            gfiles,
            nhillslope,
            nmaxhillcol,
            add_bedrock,
            add_stream_channel_vars,
        )


def write_to_file(
    outfile_path,
    n_lat,
    n_lon,
    gfiles,
    nhillslope,
    nmaxhillcol,
    add_bedrock,
    add_stream_channel_vars,
):
    """
    Write to file
    """
    outfile = Dataset(outfile_path, "w")
    outfile.creation_date = datetime.date.today().strftime("%y%m%d")

    outfile.createDimension("lsmlon", n_lon)
    outfile.createDimension("lsmlat", n_lat)
    outfile.createDimension("nhillslope", nhillslope)
    outfile.createDimension("nmaxhillcol", nmaxhillcol)

    (
        olon,
        olat,
        olon2d,
        olat2d,
        ohand,
        odtnd,
        owidth,
        oarea,
        oslop,
        oasp,
        obed,
        onhill,
        opcthill,
        ohillndx,
        ocolndx,
        odcolndx,
        ocmask,
        osdepth,
        oswidth,
        osslope,
    ) = create_variables(add_stream_channel_vars, outfile)

    # loop over gridcell files
    for gfile in gfiles:
        y1, x1 = gfile.index("j_"), gfile.index("i_")  # pylint: disable=invalid-name
        j, i = int(gfile[y1 + 2 : y1 + 5]), int(gfile[x1 + 2 : x1 + 5])

        gfile_ds = Dataset(gfile, "r")
        lon = gfile_ds.variables["longitude"][
            :,
        ]
        lat = gfile_ds.variables["latitude"][
            :,
        ]
        lon2d = gfile_ds.variables["LONGXY"][
            :,
        ]
        lat2d = gfile_ds.variables["LATIXY"][
            :,
        ]
        chunk_mask = gfile_ds.variables["chunk_mask"][
            :,
        ]

        hillslope_elev = np.asarray(
            gfile_ds.variables["hillslope_elevation"][
                :,
            ]
        )
        hillslope_dist = np.asarray(
            gfile_ds.variables["hillslope_distance"][
                :,
            ]
        )
        hillslope_width = np.asarray(
            gfile_ds.variables["hillslope_width"][
                :,
            ]
        )
        hillslope_area = np.asarray(
            gfile_ds.variables["hillslope_area"][
                :,
            ]
        )
        hillslope_slope = np.asarray(
            gfile_ds.variables["hillslope_slope"][
                :,
            ]
        )
        hillslope_aspect = np.asarray(
            gfile_ds.variables["hillslope_aspect"][
                :,
            ]
        )
        if add_bedrock:
            hillslope_bedrock = np.asarray(
                gfile_ds.variables["hillslope_bedrock_depth"][
                    :,
                ]
            )
        if add_stream_channel_vars:
            hillslope_stream_depth = np.asarray(
                gfile_ds.variables["hillslope_stream_depth"][
                    :,
                ]
            )
            hillslope_stream_width = np.asarray(
                gfile_ds.variables["hillslope_stream_width"][
                    :,
                ]
            )
            hillslope_stream_slope = np.asarray(
                gfile_ds.variables["hillslope_stream_slope"][
                    :,
                ]
            )

        nhillcolumns = gfile_ds.variables["nhillcolumns"][
            :,
        ].astype(int)
        pct_hillslope = gfile_ds.variables["pct_hillslope"][
            :,
        ]
        hillslope_index = gfile_ds.variables["hillslope_index"][
            :,
        ].astype(int)
        column_index = gfile_ds.variables["column_index"][
            :,
        ].astype(int)
        downhill_column_index = gfile_ds.variables["downhill_column_index"][
            :,
        ].astype(int)
        gfile_ds.close()

        olon[i] = lon
        olat[j] = lat
        olon2d[j, i] = lon2d
        olat2d[j, i] = lat2d

        ohand[:, j, i] = hillslope_elev
        odtnd[:, j, i] = hillslope_dist
        oarea[:, j, i] = hillslope_area
        owidth[:, j, i] = hillslope_width
        oslop[:, j, i] = hillslope_slope
        oasp[:, j, i] = hillslope_aspect
        opcthill[:, j, i] = pct_hillslope
        onhill[j, i] = np.int32(nhillcolumns)
        ohillndx[:, j, i] = hillslope_index.astype(np.int32)
        ocolndx[:, j, i] = column_index.astype(np.int32)
        odcolndx[:, j, i] = downhill_column_index.astype(np.int32)
        ocmask[j, i] = np.int32(chunk_mask)
        if add_bedrock:
            obed[:, j, i] = hillslope_bedrock

        if add_stream_channel_vars:
            osdepth[j, i] = hillslope_stream_depth
            oswidth[j, i] = hillslope_stream_width
            osslope[j, i] = hillslope_stream_slope

    outfile.close()
    print(outfile_path + " created")


def create_variables(add_stream_channel_vars, outfile):
    """
    Create variables
    """
    olon = outfile.createVariable("longitude", float, ("lsmlon",))
    olon.units = "degrees"
    olon.long_name = "longitude"

    olat = outfile.createVariable("latitude", float, ("lsmlat",))
    olat.units = "degrees"
    olat.long_name = "latitude"

    olon2d = outfile.createVariable(
        "LONGXY",
        float,
        (
            "lsmlat",
            "lsmlon",
        ),
    )
    olon2d.units = "degrees"
    olon2d.long_name = "longitude - 2d"

    olat2d = outfile.createVariable(
        "LATIXY",
        float,
        (
            "lsmlat",
            "lsmlon",
        ),
    )
    olat2d.units = "degrees"
    olat2d.long_name = "latitude - 2d"

    (
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
    ) = shared_create_variables(outfile)

    ocmask = outfile.createVariable(
        "chunk_mask",
        np.int32,
        (
            "lsmlat",
            "lsmlon",
        ),
    )
    ocmask.units = "unitless"
    ocmask.long_name = "chunk mask"

    if add_stream_channel_vars:
        wdims = outfile["LONGXY"].dimensions
        osdepth = outfile.createVariable("hillslope_stream_depth", float, wdims)
        oswidth = outfile.createVariable("hillslope_stream_width", float, wdims)
        osslope = outfile.createVariable("hillslope_stream_slope", float, wdims)

        osdepth.long_name = "stream channel bankfull depth"
        osdepth.units = "m"

        oswidth.long_name = "stream channel bankfull width"
        oswidth.units = "m"

        osslope.long_name = "stream channel slope"
        osslope.units = "m/m"

    # Initialize arrays
    olon[
        :,
    ] = 0
    olat[
        :,
    ] = 0
    olon2d[
        :,
    ] = 0
    olat2d[
        :,
    ] = 0

    ohand[
        :,
    ] = 0
    odtnd[
        :,
    ] = 0
    oarea[
        :,
    ] = 0
    owidth[
        :,
    ] = 0
    oslop[
        :,
    ] = 0
    obed[
        :,
    ] = 0
    oasp[
        :,
    ] = 0
    opcthill[
        :,
    ] = 0
    onhill[
        :,
    ] = 0
    ohillndx[
        :,
    ] = 0
    ocolndx[
        :,
    ] = 0
    odcolndx[
        :,
    ] = 0
    ocmask[
        :,
    ] = 0

    if add_stream_channel_vars:
        osdepth[
            :,
        ] = 0
        oswidth[
            :,
        ] = 0
        osslope[
            :,
        ] = 0

    return (
        olon,
        olat,
        olon2d,
        olat2d,
        ohand,
        odtnd,
        owidth,
        oarea,
        oslop,
        oasp,
        obed,
        onhill,
        opcthill,
        ohillndx,
        ocolndx,
        odcolndx,
        ocmask,
        osdepth,
        oswidth,
        osslope,
    )
