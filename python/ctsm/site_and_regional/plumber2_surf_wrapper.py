#! /usr/bin/env python3
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This script is a simple wrapper for neon sites that performs the
following:
    1) For neon sites, subset surface dataset from global dataset
        (i.e. ./subset_data.py )
    2) Download neon and update the created surface dataset
       based on the downloaded neon data.
       (i.e. modify_singlept_site_neon.py)

Instructions for running using conda python environments:

../../py_env_create
conda activate ctsm_py

"""
#  Import libraries
from __future__ import print_function

import argparse
import logging
import sys
import tqdm

# pylint:disable=wrong-import-position
from ctsm.site_and_regional.plumber2_shared import PLUMBER2_SITES_CSV, read_plumber2_sites_csv
from ctsm import subset_data
from ctsm.pft_utils import MAX_PFT_MANAGEDCROPS, is_valid_pft


def get_args():
    """
    Get arguments for this script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        "-v",
        "--verbose",
        help="Verbose mode will print more information. ",
        action="store_true",
        dest="verbose",
    )

    parser.add_argument(
        "--crop",
        help=f"Create and/or modify {MAX_PFT_MANAGEDCROPS}-PFT "
        "surface datasets (e.g. for a non-FATES run)",
        action="store_true",
        dest="use_managed_crops",
    )

    parser.add_argument(
        "--overwrite",
        help="Overwrite any existing files",
        action="store_true",
    )

    parser.add_argument(
        "--plumber2-sites-csv",
        help=f"Comma-separated value (CSV) file with Plumber2 sites. Default: {PLUMBER2_SITES_CSV}",
        default=PLUMBER2_SITES_CSV,
    )

    return parser.parse_args()


def execute(command):
    """
    Runs subset_data with given arguments.
    Args:
        command (list):
            list of args for command that we want to run.
    Raises:
        Whatever error subset_data gives, if any.
    """
    print("\n", " >>  ", *command, "\n")

    sys.argv = command
    subset_data.main()


def main():
    """
    Read plumber2_sites from csv, iterate through sites, and add dominant PFT
    """

    args = get_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    plumber2_sites = read_plumber2_sites_csv(args.plumber2_sites_csv)

    for _, row in tqdm.tqdm(plumber2_sites.iterrows()):
        lat = row["Lat"]
        lon = row["Lon"]
        site = row["Site"]

        clmsite = "1x1_PLUMBER2_" + site
        print("Now processing site :", site)

        # Set up part of subset_data command that is shared among all options
        subset_command = [
            "./subset_data",
            "point",
            "--lat",
            str(lat),
            "--lon",
            str(lon),
            "--site",
            clmsite,
            "--create-surface",
            "--uniform-snowpack",
            "--cap-saturation",
            "--lon-type",
            "180",
        ]

        # Read info for first PFT
        pft1 = row["pft1"]
        if not is_valid_pft(pft1, args.use_managed_crops):
            raise RuntimeError(f"pft1 must be a valid PFT; got {pft1}")
        pctpft1 = row["pft1-%"]
        cth1 = row["pft1-cth"]
        cbh1 = row["pft1-cbh"]

        # Read info for second PFT, if a valid one is given in the .csv file
        pft2 = row["pft2"]
        if is_valid_pft(pft2, args.use_managed_crops):
            pctpft2 = row["pft2-%"]
            cth2 = row["pft2-cth"]
            cbh2 = row["pft2-cbh"]

        # Set dominant PFT(s)
        if is_valid_pft(pft2, args.use_managed_crops):
            subset_command += [
                "--dompft",
                str(pft1),
                str(pft2),
                "--pctpft",
                str(pctpft1),
                str(pctpft2),
            ]
        else:
            subset_command += [
                "--dompft",
                str(pft1),
                "--pctpft",
                str(pctpft1),
            ]

        if not args.use_managed_crops:
            # use surface dataset with 78 pfts, but overwrite to 100% 1 dominant PFT
            # don't set crop flag
            # set canopy top and bottom heights
            if is_valid_pft(pft2, args.use_managed_crops):
                subset_command += [
                    "--cth",
                    str(cth1),
                    str(cth2),
                    "--cbh",
                    str(cbh1),
                    str(cbh2),
                ]
            else:
                subset_command += [
                    "--cth",
                    str(cth1),
                    "--cbh",
                    str(cbh1),
                ]
        else:
            # use surface dataset with 78 pfts, and overwrite to 100% 1 dominant PFT
            # NOTE: FATES will currently not run with a 78-PFT surface dataset
            # set crop flag
            subset_command += ["--crop"]
            # don't set canopy top and bottom heights

        if args.verbose:
            subset_command += ["--verbose"]
        if args.overwrite:
            subset_command += ["--overwrite"]

        execute(subset_command)


if __name__ == "__main__":
    main()
