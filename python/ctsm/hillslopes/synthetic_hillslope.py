"""
 specify a synthetic hillslope profile
"""

import argparse
import os
import sys
import numpy as np

# The below "pylint: disable" is because pylint complains that netCDF4 has no member Dataset, even
# though it does.
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm.hillslopes.hillslope_utils import HillslopeVars


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
        "-o",
        "--output-file",
        help=(
            "Output surface dataset (default: append .synth_hillslopes before extension of "
            + "--input-file)"
        ),
        default=None,
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing output file",
        dest="overwrite",
        action="store_true",
    )

    # Synthetic hillslope settings
    default = 1.0
    parser.add_argument(
        "--delx",
        help=(
            "increments to use in numerical integration of mean elevation (m) "
            + "(default: {default})"
        ),
        type=float,
        default=default,
    )
    default = "slope_aspect"
    parser.add_argument(
        "--hcase",
        help=f"hcase (default: {default})",
        type=str,
        default=default,
        choices=["slope_aspect"],
    )
    default = 500.0
    parser.add_argument(
        "--hillslope-distance",
        help=f"distance from channel to ridge (m) (default: {default})",
        type=float,
        default=default,
    )
    default = 16
    parser.add_argument(
        "--nmaxhillcol",
        help=f"max. number of hillslope columns (default: {default})",
        type=int,
        default=default,
    )
    default = 4
    parser.add_argument(
        "--num-hillslopes",
        help=f"number of hillslopes (default: {default})",
        type=int,
        default=default,
    )
    default = 1.0
    parser.add_argument(
        "--phill",
        help=f"shape parameter (power law exponent) (default: {default})",
        type=float,
        default=default,
    )
    default = 2.0
    parser.add_argument(
        "--thresh",
        help=f"threshold for freating specified fractional bins (default: {default})",
        type=float,
        default=default,
    )
    default = 500.0
    parser.add_argument(
        "--width-reach",
        help=f"uniform width of reach (m) (default: {default})",
        type=float,
        default=default,
    )

    args = parser.parse_args(argv)

    if args.output_file is None:
        stem, ext = os.path.splitext(args.input_file)
        args.output_file = stem + ".synth_hillslopes" + ext

    if os.path.exists(args.output_file) and not args.overwrite:
        raise FileExistsError(f"Output file already exists: {args.output_file}")

    return args


def calc_stream_geom(hillslope_vars):
    """
    Calculate stream geometry variables
    """
    uharea = np.sum(hillslope_vars.h_area, axis=0)
    adepth, bdepth = 1e-3, 0.4
    hillslope_vars.h_stream_depth[
        :,
    ] = adepth * (uharea**bdepth)
    awidth, bwidth = 1e-3, 0.6
    hillslope_vars.h_stream_width[
        :,
    ] = awidth * (uharea**bwidth)
    hillslope_vars.h_stream_slope[
        :,
    ] = 1e-2
    return hillslope_vars


"""
---------------------------------------------------
#cosine - power law hillslope
create bins of equal height
this form ensures a near-zero slope at the hill top
---------------------------------------------------
"""  # pylint: disable=pointless-string-statement


def cosp_height(x, hlen, hhgt, phill):
    """
    Get elevation
    """
    fx = 0.5 * (1.0 + np.cos(np.pi * (1.0 + (x / hlen))))
    h = hhgt * np.power(fx, phill)
    return h


def icosp_height(h, hlen, hhgt, phill):
    """
    Fill lbins
    """
    if hhgt <= 0.0:
        x = 0.0
    else:
        fh = np.arccos(2.0 * np.power((h / hhgt), (1.0 / phill)) - 1)
        # np.arccos returns [0,pi]
        # want [pi,2pi] based on cosp_height definition
        fh = 2.0 * np.pi - fh
        x = hlen * ((1.0 / np.pi) * fh - 1.0)
    return x


def calc_mean_elevation(args, hhgt, uedge, ledge):
    """
    numerically integrate to calculate mean elevation
    """
    nx = int(uedge - ledge)
    mean_elev = 0.0
    for k in range(nx):
        x1 = uedge - (k + 0.5) * args.delx
        mean_elev += cosp_height(x1, args.hillslope_distance, hhgt, args.phill)
    mean_elev = mean_elev / float(nx)
    return mean_elev


def create_bins(args, max_columns_per_hillslope, bin_fractions, hhgt):
    """
    create specified fractional bins
    """

    # create height bins
    hbins = np.zeros(max_columns_per_hillslope + 1)
    hbins[1] = args.thresh
    # array needs to be length max_columns_per_hillslope-1
    hbins[2 : max_columns_per_hillslope + 1] = hbins[1] + (hhgt - args.thresh) * bin_fractions

    # create length bins from height bins
    lbins = np.zeros(max_columns_per_hillslope + 1)
    for k in range(max_columns_per_hillslope + 1):
        if hhgt > 0.0:
            lbins[k] = icosp_height(hbins[k], args.hillslope_distance, hhgt, args.phill)
    return hbins, lbins


def main():
    """
    See module description
    """

    args = parse_arguments(sys.argv[1:])

    with Dataset(args.input_file, "r") as infile:
        n_lon = len(infile.dimensions["lsmlon"])
        n_lat = len(infile.dimensions["lsmlat"])
        std_elev = np.asarray(infile.variables["STD_ELEV"][:, :])
        lfrac = np.asarray(infile.variables["LANDFRAC_PFT"][:, :])
        lmask = lfrac > 0
        pct_natveg = np.asarray(infile.variables["PCT_NATVEG"][:, :])

    # are any points in land mask but have zero % natveg?
    print("zero natveg pts ", np.sum(np.where(np.logical_and(lmask == 1, pct_natveg == 0), 1, 0)))
    lmask = np.where(np.logical_and(lmask == 1, pct_natveg > 0), 1, 0).astype(int)

    max_columns_per_hillslope = args.nmaxhillcol // args.num_hillslopes

    if max_columns_per_hillslope == 4:
        bin_fractions = np.array((0.25, 0.75, 1.0))
    elif max_columns_per_hillslope == 5:
        bin_fractions = np.array((0.25, 0.50, 0.75, 1.0))
    elif max_columns_per_hillslope == 6:
        bin_fractions = np.array((0.20, 0.40, 0.60, 0.80, 1.0))
    else:
        raise RuntimeError(f"Unhandled max_columns_per_hillslope: {max_columns_per_hillslope}")

    max_columns_per_landunit = args.num_hillslopes * max_columns_per_hillslope

    hillslope_vars = HillslopeVars(max_columns_per_landunit, args.num_hillslopes, n_lat, n_lon, recurse=False, incl_latlon=True, incl_pftndx=True)

    cndx = 0
    for i in range(n_lon):
        for j in range(n_lat):
            if lmask[j, i] != 1:
                continue

            # slope tangent (y/x)
            beta = np.min((std_elev[j, i], 200.0)) / args.hillslope_distance

            # specify hill height from slope and length
            hhgt = beta * args.hillslope_distance
            hhgt = np.max([hhgt, 4.0])

            # create specified fractional bins
            hbins, lbins = create_bins(args, max_columns_per_hillslope, bin_fractions, hhgt)

            # loop over aspect bins
            for naspect in range(args.num_hillslopes):
                hillslope_vars.pct_hillslope[naspect, j, i] = 100 / float(args.num_hillslopes)
                # index from ridge to channel (i.e. downhill)
                for k in range(max_columns_per_hillslope):
                    ncol = k + naspect * max_columns_per_hillslope

                    cndx += 1  # start at 1 not zero (oceans are 0)
                    hillslope_vars.column_index[ncol, j, i] = cndx
                    hillslope_vars.hillslope_index[ncol, j, i] = naspect + 1

                    uedge = lbins[k + 1]
                    ledge = lbins[k]
                    #      lowland column
                    if k == 0:
                        hillslope_vars.downhill_column_index[ncol, j, i] = -999
                    else:  # upland columns
                        hillslope_vars.downhill_column_index[ncol, j, i] = hillslope_vars.column_index[ncol, j, i] - 1

                    hillslope_vars.h_dist[ncol, j, i] = 0.5 * (uedge + ledge)
                    hillslope_vars.h_area[ncol, j, i] = args.width_reach * (uedge - ledge)
                    hillslope_vars.h_width[ncol, j, i] = args.width_reach

                    # numerically integrate to calculate mean elevation
                    hillslope_vars.h_elev[ncol, j, i] = calc_mean_elevation(args, hhgt, uedge, ledge)

                    hillslope_vars.h_slope[ncol, j, i] = (hbins[k + 1] - hbins[k]) / (lbins[k + 1] - lbins[k])
                    if 0 <= naspect <= 3:
                        # 0 = north
                        # 1 = east
                        # 2 = south
                        # 3 = west
                        hillslope_vars.h_aspect[ncol, j, i] = naspect * np.pi / 2
                    else:
                        raise RuntimeError(f"Unhandled naspect: {naspect}")

    # Fill stream geometry variables
    hillslope_vars = calc_stream_geom(hillslope_vars)

    # Fill other variables
    hillslope_vars.nhillcolumns = max_columns_per_landunit * lmask
    hillslope_vars.h_bedrock = 2.0
    hillslope_vars.hillslope_pftndx[:] = 13

    # write to file  --------------------------------------------
    hillslope_vars.save(
        args.input_file,
        args.output_file,
        max_columns_per_landunit,
        args.num_hillslopes,
        add_bedrock=True,
        add_stream=True,
    )

    # Save settings as global attributes
    with Dataset(args.output_file, "a") as outfile:
        outfile.synth_hillslopes_delx = args.delx
        outfile.synth_hillslopes_hcase = args.hcase
        outfile.synth_hillslopes_hillslope_distance = args.hillslope_distance
        outfile.synth_hillslopes_nmaxhillcol = np.int32(args.nmaxhillcol)
        outfile.synth_hillslopes_num_hillslopes = np.int32(args.num_hillslopes)
        outfile.synth_hillslopes_phill = args.phill
        outfile.synth_hillslopes_thresh = args.thresh
        outfile.synth_hillslopes_width_reach = args.width_reach
