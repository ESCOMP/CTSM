"""
Combine chunk files into single surface data file
"""
import subprocess
import argparse
import sys
import os
import netCDF4 as netcdf4
import numpy as np


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

    parser.add_argument("-v", "--verbose", help="print info", action="store_true", default=False)

    args = parser.parse_args(argv)

    # Check arguments
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")
    if not os.path.exists(args.input_dir):
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")
    if os.path.exists(args.output_file) and not args.overwrite:
        raise FileExistsError(f"Output file already exists: {args.output_dir}")

    return args


def main():
    """
    See module description
    """
    args = parse_arguments(sys.argv[1:])

    # Choose data files to combine and append
    cfile0 = os.path.join(
        args.input_dir, f"combined_chunk_ChunkIndex_HAND_4_col_hillslope_geo_params_section_quad_{args.dem_source}.nc"
    )

    f = netcdf4.Dataset(args.input_file, "r")
    landmask = f.variables["PFTDATA_MASK"][
        :,
    ]
    f.close()

    initializeArrays = True
    for n in range(args.n_chunks):
        print(n, args.n_chunks)
        ndx = n + 1
        cstr = "{:02d}".format(ndx)
        cfile = cfile0.replace("ChunkIndex", cstr)
        command = ["ls", cfile]
        file_exists = subprocess.run(command, capture_output=True).returncode

        if initializeArrays and file_exists == 0:
            f = netcdf4.Dataset(cfile, "r")

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

            if "hillslope_bedrock_depth" in f.variables.keys():
                addBedrock = True
            else:
                addBedrock = False

            if "hillslope_stream_depth" in f.variables.keys():
                addStream = True
            else:
                addStream = False

            initializeArrays = False
            f.close()

        if file_exists > 0:
            if args.verbose:
                print(f"Skipping; chunk file not found: {cfile}")
            continue
        
        f = netcdf4.Dataset(cfile, "r")
        nhillslope = len(f.dimensions["nhillslope"])
        chunk_mask = f.variables["chunk_mask"][
            :,
        ]
        h_elev0 = f.variables["hillslope_elevation"][
            :,
        ]
        h_dist0 = f.variables["hillslope_distance"][
            :,
        ]
        h_width0 = f.variables["hillslope_width"][
            :,
        ]
        h_area0 = f.variables["hillslope_area"][
            :,
        ]
        h_slope0 = f.variables["hillslope_slope"][
            :,
        ]
        h_aspect0 = f.variables["hillslope_aspect"][
            :,
        ]
        if addBedrock:
            h_bedrock0 = f.variables["hillslope_bedrock_depth"][
                :,
            ]
        if addStream:
            h_stream_depth0 = f.variables["hillslope_stream_depth"][
                :,
            ]
            h_stream_width0 = f.variables["hillslope_stream_width"][
                :,
            ]
            h_stream_slope0 = f.variables["hillslope_stream_slope"][
                :,
            ]

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
                if chunk_mask[j, i] > 0 and landmask[j, i] > 0:

                    h_elev[:, j, i] = h_elev0[:, j, i]
                    h_dist[:, j, i] = h_dist0[:, j, i]
                    h_width[:, j, i] = h_width0[:, j, i]
                    h_area[:, j, i] = h_area0[:, j, i]
                    h_slope[:, j, i] = h_slope0[:, j, i]
                    h_aspect[:, j, i] = h_aspect0[:, j, i]
                    if addBedrock:
                        h_bedrock[:, j, i] = h_bedrock0[:, j, i]
                    if addStream:
                        h_stream_depth[j, i] = h_stream_depth0[j, i]
                        h_stream_width[j, i] = h_stream_width0[j, i]
                        h_stream_slope[j, i] = h_stream_slope0[j, i]

                    nhillcolumns[j, i] = nhillcolumns0[j, i]
                    pct_hillslope[:, j, i] = pct_hillslope0[:, j, i]
                    hillslope_index[:, j, i] = hillslope_index0[:, j, i]
                    column_index[:, j, i] = column_index0[:, j, i]
                    downhill_column_index[:, j, i] = downhill_column_index0[:, j, i]

                    # if 2 or less valid aspects, remove all hillslope data
                    minAspect = True
                    if minAspect:
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

    command = 'date "+%y%m%d"'
    timetag = (
        subprocess.Popen(command, stdout=subprocess.PIPE, shell="True")
        .communicate()[0]
        .strip()
        .decode()
    )

    # copy original file
    command = ["cp", args.input_file, args.output_file]
    x = subprocess.call(command, stderr=subprocess.PIPE)

    print("appending file")
    w = netcdf4.Dataset(args.output_file, "a")
    w.creation_date = timetag

    w.createDimension("nhillslope", nhillslope)
    w.createDimension("nmaxhillcol", ncolumns_per_gridcell)

    ohand = w.createVariable(
        "hillslope_elevation",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    ohand.units = "m"
    ohand.long_name = "hillslope elevation above channel"

    odtnd = w.createVariable(
        "hillslope_distance",
        np.float64,
        (
            "nmaxhillcol",
            "lsmlat",
            "lsmlon",
        ),
    )
    odtnd.units = "m"
    odtnd.long_name = "hillslope  distance from channel"

    owidth = w.createVariable(
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

    oarea = w.createVariable(
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

    oslop = w.createVariable(
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

    oasp = w.createVariable(
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

    if addBedrock:
        obed = w.createVariable(
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

    if addStream:
        osdepth = w.createVariable(
            "hillslope_stream_depth",
            np.float64,
            (
                "lsmlat",
                "lsmlon",
            ),
        )
        osdepth.units = "meters"
        osdepth.long_name = "stream channel bankfull depth"
        oswidth = w.createVariable(
            "hillslope_stream_width",
            np.float64,
            (
                "lsmlat",
                "lsmlon",
            ),
        )
        oswidth.units = "meters"
        oswidth.long_name = "stream channel bankfull width"
        osslope = w.createVariable(
            "hillslope_stream_slope",
            np.float64,
            (
                "lsmlat",
                "lsmlon",
            ),
        )
        osslope.units = "m/m"
        osslope.long_name = "stream channel slope"

    onhill = w.createVariable(
        "nhillcolumns",
        np.int32,
        (
            "lsmlat",
            "lsmlon",
        ),
    )
    onhill.units = "unitless"
    onhill.long_name = "number of columns per landunit"

    opcthill = w.createVariable(
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

    ohillndx = w.createVariable(
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

    ocolndx = w.createVariable(
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

    odcolndx = w.createVariable(
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

    ohand[
        :,
    ] = h_elev
    odtnd[
        :,
    ] = h_dist
    oarea[
        :,
    ] = h_area
    owidth[
        :,
    ] = h_width
    oslop[
        :,
    ] = h_slope
    # aspect should be in radians on surface data file
    oasp[
        :,
    ] = h_aspect
    if addBedrock:
        obed[
            :,
        ] = h_bedrock
    if addStream:
        osdepth[
            :,
        ] = h_stream_depth
        oswidth[
            :,
        ] = h_stream_width
        osslope[
            :,
        ] = h_stream_slope
    opcthill[
        :,
    ] = pct_hillslope
    onhill[
        :,
    ] = nhillcolumns.astype(np.int32)
    ohillndx[
        :,
    ] = hillslope_index.astype(np.int32)
    ocolndx[
        :,
    ] = column_index.astype(np.int32)
    odcolndx[
        :,
    ] = downhill_column_index.astype(np.int32)

    w.close()
    print(args.output_file + " created")
