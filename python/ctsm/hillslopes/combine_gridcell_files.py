#!/usr/bin/env python
# coding: utf-8

# Import modules

import sys
import string
import subprocess
import copy
import argparse
import netCDF4 as netcdf4
import numpy as np

# ---------------------------------------------------------- #
#
#  Combine gridcell files into single file
#
# ---------------------------------------------------------- #

parser = argparse.ArgumentParser(description="Combine gridcell files")
parser.add_argument("cndx", help="chunk", nargs="?", type=int, default=0)
parser.add_argument("-o", "--overwrite", help="overwrite", action="store_true", default=False)
parser.add_argument("-v", "--verbose", help="print info", action="store_true", default=False)
args = parser.parse_args()

cndx = args.cndx
verbose = args.verbose

totalChunks = 36
if cndx < 1 or cndx > totalChunks:
    print("cndx must be 1-{:d}".format(totalChunks))
    stop

printFlush = True

# Specify data files to combine

# Input file with grid information
gridfile = "/fs/cgd/csm/inputdata/lnd/clm2/surfdata_map/surfdata_0.9x1.25_78pfts_CMIP6_simyr2000_c170824.nc"

# Gridcell file directory
dem_source = "MERIT"
cdir = "./{}_gridcell_files/".format(dem_source)
cfile = cdir + "chunk_{:02d}_HAND_4_col_hillslope_geo_params_section_quad_{}.nc".format(
    cndx, dem_source
)

# Output file
odir = "./"
outfile = odir + cfile.split("/")[-1].replace("chunk_", "combined_chunk_")

# Read output file coordinates
f = netcdf4.Dataset(gridfile, "r")
sjm = len(f.dimensions["lsmlat"])
sim = len(f.dimensions["lsmlon"])
f.close()

# Check for output file existence
command = ["ls", outfile]
file_exists = subprocess.run(command, capture_output=True).returncode
# if file_exists !=0, file not found
if file_exists == 0:
    if args.overwrite:
        if verbose:
            print(outfile, " exists; overwriting", flush=printFlush)
    else:
        print(outfile, " exists; stopping", flush=printFlush)
        stop

# Locate gridcell files
gfile = cfile.replace(".nc", "*.nc")

command = "ls " + gfile
gfiles = (
    subprocess.run(command, capture_output=True, shell="True")
    .stdout.decode("utf-8")
    .split("\n")[:-1]
)

if len(gfiles) == 0:
    print("No files found")
    stop

# Read hillslope data dimensions
f = netcdf4.Dataset(gfiles[0], "r")
nhillslope = len(f.dimensions["nhillslope"])
nmaxhillcol = len(f.dimensions["nmaxhillcol"])

if "hillslope_bedrock_depth" in f.variables.keys():
    addBedrock = True
else:
    addBedrock = False

if "hillslope_stream_depth" in f.variables.keys():
    addStreamChannelVariables = True
else:
    addStreamChannelVariables = False
f.close()

# initialize outfile
command = 'date "+%y%m%d"'
timetag = (
    subprocess.Popen(command, stdout=subprocess.PIPE, shell="True")
    .communicate()[0]
    .strip()
    .decode()
)

w = netcdf4.Dataset(outfile, "w")
w.creation_date = timetag

w.createDimension("lsmlon", sim)
w.createDimension("lsmlat", sjm)
w.createDimension("nhillslope", nhillslope)
w.createDimension("nmaxhillcol", nmaxhillcol)

olon = w.createVariable("longitude", float, ("lsmlon",))
olon.units = "degrees"
olon.long_name = "longitude"

olat = w.createVariable("latitude", float, ("lsmlat",))
olat.units = "degrees"
olat.long_name = "latitude"

olon2d = w.createVariable(
    "LONGXY",
    float,
    (
        "lsmlat",
        "lsmlon",
    ),
)
olon2d.units = "degrees"
olon2d.long_name = "longitude - 2d"

olat2d = w.createVariable(
    "LATIXY",
    float,
    (
        "lsmlat",
        "lsmlon",
    ),
)
olat2d.units = "degrees"
olat2d.long_name = "latitude - 2d"

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
odtnd.long_name = "hillslope distance from channel"

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

ocmask = w.createVariable(
    "chunk_mask",
    np.int32,
    (
        "lsmlat",
        "lsmlon",
    ),
)
ocmask.units = "unitless"
ocmask.long_name = "chunk mask"

if addStreamChannelVariables:
    wdims = w["LONGXY"].dimensions
    osdepth = w.createVariable("hillslope_stream_depth", float, wdims)
    oswidth = w.createVariable("hillslope_stream_width", float, wdims)
    osslope = w.createVariable("hillslope_stream_slope", float, wdims)

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

if addStreamChannelVariables:
    osdepth[
        :,
    ] = 0
    oswidth[
        :,
    ] = 0
    osslope[
        :,
    ] = 0

# loop over gridcell files
for gfile in gfiles:

    y1, x1 = gfile.index("j_"), gfile.index("i_")
    j, i = int(gfile[y1 + 2 : y1 + 5]), int(gfile[x1 + 2 : x1 + 5])

    f = netcdf4.Dataset(gfile, "r")
    lon = f.variables["longitude"][
        :,
    ]
    lat = f.variables["latitude"][
        :,
    ]
    lon2d = f.variables["LONGXY"][
        :,
    ]
    lat2d = f.variables["LATIXY"][
        :,
    ]
    chunk_mask = f.variables["chunk_mask"][
        :,
    ]
    # h_elev = f.variables['hillslope_elevation'][:,]
    # h_dist = f.variables['hillslope_distance'][:,]
    # h_width = f.variables['hillslope_width'][:,]
    # h_area = f.variables['hillslope_area'][:,]
    # h_slope = f.variables['hillslope_slope'][:,]
    # h_aspect = f.variables['hillslope_aspect'][:,]
    # if addBedrock:
    #     h_bedrock = f.variables['hillslope_bedrock_depth'][:,]
    # if addStreamChannelVariables:
    #     h_stream_depth = f.variables['hillslope_stream_depth'][:,]
    #     h_stream_width = f.variables['hillslope_stream_width'][:,]
    #     h_stream_slope = f.variables['hillslope_stream_slope'][:,]

    h_elev = np.asarray(
        f.variables["h_height"][
            :,
        ]
    )
    h_dist = np.asarray(
        f.variables["h_length"][
            :,
        ]
    )
    h_width = np.asarray(
        f.variables["h_width"][
            :,
        ]
    )
    h_area = np.asarray(
        f.variables["h_area"][
            :,
        ]
    )
    h_slope = np.asarray(
        f.variables["h_slope"][
            :,
        ]
    )
    h_aspect = np.asarray(
        f.variables["h_aspect"][
            :,
        ]
    )
    if addBedrock:
        h_bedrock = np.asarray(
            f.variables["h_bedrock"][
                :,
            ]
        )
    if addStreamChannelVariables:
        h_stream_depth = np.asarray(
            f.variables["h_stream_depth"][
                :,
            ]
        )
        h_stream_width = np.asarray(
            f.variables["h_stream_width"][
                :,
            ]
        )
        h_stream_slope = np.asarray(
            f.variables["h_stream_slope"][
                :,
            ]
        )

    nhillcolumns = f.variables["nhillcolumns"][
        :,
    ].astype(int)
    pct_hillslope = f.variables["pct_hillslope"][
        :,
    ]
    hillslope_index = f.variables["hillslope_index"][
        :,
    ].astype(int)
    column_index = f.variables["column_index"][
        :,
    ].astype(int)
    downhill_column_index = f.variables["downhill_column_index"][
        :,
    ].astype(int)
    f.close()

    olon[i] = lon
    olat[j] = lat
    olon2d[j, i] = lon2d
    olat2d[j, i] = lat2d

    ohand[:, j, i] = h_elev
    odtnd[:, j, i] = h_dist
    oarea[:, j, i] = h_area
    owidth[:, j, i] = h_width
    oslop[:, j, i] = h_slope
    oasp[:, j, i] = h_aspect
    opcthill[:, j, i] = pct_hillslope
    onhill[j, i] = np.int32(nhillcolumns)
    ohillndx[:, j, i] = hillslope_index.astype(np.int32)
    ocolndx[:, j, i] = column_index.astype(np.int32)
    odcolndx[:, j, i] = downhill_column_index.astype(np.int32)
    ocmask[j, i] = np.int32(chunk_mask)
    if addBedrock:
        obed[:, j, i] = h_bedrock

    if addStreamChannelVariables:
        osdepth[j, i] = h_stream_depth
        oswidth[j, i] = h_stream_width
        osslope[j, i] = h_stream_slope


w.close()
print(outfile + " created")
