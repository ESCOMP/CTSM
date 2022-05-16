#! /usr/bin/env python
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This script is for modifying surface dataset at neon sites
using data available from the neon server.

After creating a single point surface data file from a global 
surface data file using subset_data.py, use this script to 
overwrite some fields with site-specific data for neon sites.

This script will do the following:
- Download neon data for the specified site if it does not exist 
    in the specified directory : (i.e. ../../../neon_surffiles).
- Modify surface dataset with downloaded data.

-------------------------------------------------------------------
Instructions for running on Cheyenne/Casper:

load the following into your local environment
    module load python
    ncar_pylib

To remove NPL from your environment on Cheyenne/Casper:
    deactivate
-------------------------------------------------------------------
To see the available options:
    ./modify_singlept_site_neon.py --help
-------------------------------------------------------------------
Example:
    ./modify_singlept_site_neon.py --neon_site PUUM --debug
-------------------------------------------------------------------
"""
# TODO (NS)
# --[] If subset file not found run subset_data.py
# --[] List of valid neon sites for all scripts come from one place.
# --[] Download files only when available.

#  Import libraries
from __future__ import print_function

import os
import sys
import glob
import argparse
import requests

import logging
import numpy as np
import pandas as pd
import xarray as xr

from datetime import date
from getpass import getuser


myname = getuser()


# -- valid neon sites
valid_neon_sites = [
    "ABBY",
    "BARR",
    "BART",
    "BLAN",
    "BONA",
    "CLBJ",
    "CPER",
    "DCFS",
    "DEJU",
    "DELA",
    "DSNY",
    "GRSM",
    "GUAN",
    "HARV",
    "HEAL",
    "JERC",
    "JORN",
    "KONA",
    "KONZ",
    "LAJA",
    "LENO",
    "MLBS",
    "MOAB",
    "NIWO",
    "NOGP",
    "OAES",
    "ONAQ",
    "ORNL",
    "OSBS",
    "PUUM",
    "RMNP",
    "SCBI",
    "SERC",
    "SJER",
    "SOAP",
    "SRER",
    "STEI",
    "STER",
    "TALL",
    "TEAK",
    "TOOL",
    "TREE",
    "UKFS",
    "UNDE",
    "WOOD",
    "WREF",
    "YELL",
]


def get_parser():
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        "--neon_site",
        help="4-letter neon site code.",
        action="store",
        dest="site_name",
        choices=valid_neon_sites,
        required=True,
    )
    parser.add_argument(
        "--surf_dir",
        help="""
                Directory of single point surface dataset. 
                [default: %(default)s]
                """,
        action="store",
        dest="surf_dir",
        type=str,
        required=False,
        default="/glade/scratch/" + myname + "/single_point/",
    )
    parser.add_argument(
        "--out_dir",
        help="""
                Directory to write updated single point surface dataset.
                [default: %(default)s] 
                """,
        action="store",
        dest="out_dir",
        type=str,
        required=False,
        default="/glade/scratch/" + myname + "/single_point_neon_updated/",
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="""
                Debug mode will print more information. 
                [default: %(default)s] 
                """,
        action="store_true",
        dest="debug",
        default=False,
    )

    return parser


def get_neon(neon_dir, site_name):
    """
    Function for finding neon data files
    and download from neon server if the
    file does not exits.

    Args:
        neon_dir (str): local directory for downloading neon data.
        site_name (str): 4 letter neon site name

    Raises:
        Error if the download was not successful (exit code:404).
        In case the data does not exist in the neon server or if
        neon server is down.

    Returns:
        neon_file (str) : complete file name of the downloaded data
    """

    # -- create directory if not exists
    if not os.path.exists(neon_dir):
        os.makedirs(neon_dir)

    neon_file = os.path.join(neon_dir, site_name + "_surfaceData.csv")

    # -- Download the file if it does not exits
    if os.path.isfile(neon_file):
        print("neon file for", site_name, "already exists! ")
        print("Skipping download from neon for", site_name, "...")
    else:
        print("------------------------------------------------")
        print("Beginning download from neon server for", site_name, "...")

        url = (
            "https://s3.data.neonscience.org/neon-ncar/NEON/surf_files/v1/"
            + site_name
            + "_surfaceData.csv"
        )
        response = requests.get(url)

        with open(neon_file, "wb") as f:
            f.write(response.content)

        # -- Check if download status_code
        if response.status_code == 200:
            print("Download finished successfully for", site_name)
        elif response.status_code == 404:
            sys.exit(
                "Data for this site "
                + site_name
                + " was not available on the neon server:"
                + url
            )

        print("Download exit status code:  ", response.status_code)
        print("Downloaded file type     :  ", response.headers["content-type"])
        print("Downloaded file encoding :  ", response.encoding)
        print("------------------------------------------------")

        response.close()

    return neon_file


def find_surffile(surf_dir, site_name):
    """
    Function for finding and choosing surface file for
    a neon site.
    These files are created using ./subset_data.py script.
    In case multiple files exist for the neon site, it
    will choose the file created the latest.

    Args:
        surf_dir (str): directory of single point surface data
        site_name (str): 4 letter neon site name

    Raises:
        Error if the surface data for the site is not created

    Returns:
        surf_file (str): name of the surface dataset file
    """

    # sf_name = "surfdata_hist_16pfts_Irrig_CMIP6_simyr2000_"+site_name+"*.nc"
    sf_name = "surfdata_*hist_78pfts_CMIP6_simyr2000_" + site_name + "*.nc"
    print (os.path.join(surf_dir , sf_name))
    surf_file = sorted(glob.glob(os.path.join(surf_dir , sf_name)))

    if len(surf_file) > 1:
        print("The following files found :", *surf_file, sep="\n- ")
        print("The latest file is chosen :", surf_file[-1])
        surf_file = surf_file[-1]
    elif len(surf_file) == 1:
        print("File found : ")
        print(surf_file)
        surf_file = surf_file[0]
    else:
        sys.exit(
            "Surface data for this site " + site_name + "was not found:" + surf_file,
            ".",
            "\n",
            "Please run ./subset_data.py for this site.",
        )
    return surf_file


def find_soil_structure(surf_file):
    """
    Function for finding surface dataset soil
    strucutre using surface data metadata.

    In CLM surface data, soil layer information
    is in a file from surface data metadata
    under "Soil_texture_raw_data_file_name".
    This function finds this file for the surface
    dataset, read it, and find soil layers.

    Args:
        surf_file (str): single point surface data filename

    Raises:
        error if the soil layer strucutre file does not exist

    Returns:
        soil_bot : array of soil layers top depths
        soil_top : array of soil layers bottom depths
    """
    # TODO: What if not cheyenne? Self-contained depth info.

    print("------------")
    print("surf_file : ", surf_file)
    f1 = xr.open_dataset(surf_file)
    print("------------")
    # print (f1.attrs["Soil_texture_raw_data_file_name"])

    clm_input_dir = "/glade/p/cesmdata/cseg/inputdata/lnd/clm2/rawdata/"
    surf_soildepth_file = os.path.join(
        clm_input_dir, f1.attrs["Soil_texture_raw_data_file_name"]
    )

    if os.path.exists(surf_soildepth_file):
        print(
            "\n\n Reading",
            surf_soildepth_file,
            "for surface data soil structure information:",
        )
        f1_soildepth = xr.open_dataset(surf_soildepth_file)
        print(f1_soildepth["DZSOI"])
        soil_bot = f1_soildepth["DZSOI"].values

        # -- soil layer top
        soil_top = soil_bot[:-1]
        soil_top = np.insert(soil_top, 0, 0)

    else:
        sys.exit(
            "Cannot find soil structure file : "
            + surf_soildepth_file
            + "for the surface dataset."
        )

    return soil_bot, soil_top


def update_metadata(nc, surf_file, neon_file, zb_flag):
    """
    Function for updating modified surface dataset
    metadat for neon sites.

    Args:
        nc (xr Dataset): netcdf file including updated neon surface data
        surf_file (str): single point surface data filename
        neon_file (str): filename of neon downloaded surface dataset

    Returns:
        nc (xr Dataset): netcdf file including updated neon surface data
    """
    today = date.today()
    today_string = today.strftime("%Y-%m-%d")

    nc.attrs["Updated_on"] = today_string
    nc.attrs["Updated_by"] = myname
    nc.attrs["Updated_with"] = os.path.abspath(__file__)
    nc.attrs["Updated_from"] = surf_file
    nc.attrs["Updated_using"] = neon_file
    if zb_flag:
        nc.attrs["Updated_fields"] = "PCT_CLAY, PCT_SAND, ORGANIC, zbedrock"
    else:
        nc.attrs["Updated_fields"] = "PCT_CLAY, PCT_SAND, ORGANIC"

    return nc


def update_time_tag(fname_in):
    """
    Function for updating time tag on surface dataset
    files.
    Expects file to end with [._]cYYMMDD.nc or [._]YYMMDD.nc
    Add the tag to just before that ending part.

    Args:
        fname_in (str) : file name with the old time tag

    Raises:
        error if the file does not end with with
         [._]cYYMMDD.nc or [._]YYMMDD.nc

    Returns:
        fname_out (str) : file name with the updated time tag
    """
    today = date.today()
    today_string = today.strftime("%y%m%d")

    basename = os.path.basename(fname_in)
    cend = -10
    if basename[cend] == "c":
        cend = cend - 1
    if (basename[cend] != ".") and (basename[cend] != "_"):
        sys.exit("Trouble figuring out where to add tag to filename:" + fname_in)

    fname_out = basename[:cend] + "_" + "c" + today_string + ".nc"
    return fname_out


def sort_print_soil_layers(obs_bot, soil_bot):
    """
    Function for pretty printing soil structure of
    original surface dataset and neon dataset.

    Args:
        obs_bot  : array of neon soil layers bottom depths
        soil_bot : array of soil layers bottom depths
    """

    obs_bot_df = pd.DataFrame({"depth": obs_bot, "type": "obs"})
    soil_bot_df = pd.DataFrame({"depth": soil_bot, "type": "sfc"})
    depth_df = pd.concat([obs_bot_df, soil_bot_df])

    depth_df = depth_df.sort_values("depth")

    space = " "
    print("================================", "================================")

    print("  Neon data soil structure:     ", "  Surface data soil structure:  ")

    print("================================", "================================")

    for index, row in depth_df.iterrows():
        if row["type"] == "obs":
            print("-------------", "{0:.3f}".format(row["depth"]), "------------")
        else:
            print(
                33 * space + "-------------",
                "{0:.3f}".format(row["depth"]),
                "-----------",
            )

    print("--------------------------------" + "--------------------------------")


def check_neon_time():
    """
    A function to download and parse neon listing file.

    Returns:
        dict_out (str) :
            dictionary of *_surfaceData.csv files with the last modified
    """
    listing_file = "listing.csv"
    url = "https://storage.neonscience.org/neon-ncar/listing.csv"

    download_file(url, listing_file)

    df = pd.read_csv(listing_file)
    df = df[df["object"].str.contains("_surfaceData.csv")]
    # df=df.join(df['object'].str.split("/", expand=True))
    dict_out = dict(zip(df["object"], df["last_modified"]))
    print(dict_out)
    return dict_out


def download_file(url, fname):
    """
    Function to download a file.
    Args:
        url (str):
            url of the file for downloading
        fname (str) :
            file name to save the downloaded file.
    """
    try:
        response = requests.get(url)

        with open(fname, "wb") as f:
            f.write(response.content)

        # -- Check if download status_code
        if response.status_code == 200:
            print("Download finished successfully for", fname, ".")
        elif response.status_code == 404:
            print("File " + fname + "was not available on the neon server:" + url)
    except Exception as err:
        print ('The server could not fulfill the request.')
        print ('Something went wrong in downloading', fname)
        print ('Error code:', err.code)


def fill_interpolate(f2, var, method):
    """
    Function to interpolate a variable in a
    xarray dataset a specific method
    """
    print("=====================================")
    print("Filling in ", var, "with interpolation (method =" + method + ").")

    print("Variable before filling : ")
    print(f2[var])

    tmp_df = pd.DataFrame(f2[var].values.ravel())

    tmp_df = tmp_df.interpolate(method=method, limit_direction="both")
    # tmp_df = tmp_df.interpolate(method ='spline',order = 2,   limit_direction ='both')
    # tmp_df = tmp_df.interpolate(method="pad", limit=5, limit_direction = 'forward')

    tmp = tmp_df.to_numpy()

    soil_levels = f2[var].size
    for soil_lev in range(soil_levels):
        f2[var][soil_lev] = tmp[soil_lev].reshape(1, 1)

    print("Variable after filling : ")
    print(f2[var])
    print("=====================================")


def main():

    args = get_parser().parse_args()

    # -- debugging option
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    file_time = check_neon_time()

    # --  specify site from which to extract data
    site_name = args.site_name

    # --  Look for surface data
    surf_dir = args.surf_dir
    surf_file = find_surffile(surf_dir, site_name)

    # --  directory structure
    current_dir = os.getcwd()
    parent_dir = os.path.dirname(current_dir)
    clone_dir = os.path.abspath(os.path.join(__file__, "../../.."))
    neon_dir = os.path.join(clone_dir, "neon_surffiles")

    print("Present Directory", current_dir)

    # --  download neon data if needed
    neon_file = get_neon(neon_dir, site_name)

    # -- Read neon data
    df = pd.read_csv(neon_file)

    # -- Read surface dataset files
    print("surf_file:", surf_file)
    f1 = xr.open_dataset(surf_file)

    # -- Find surface dataset soil depth information
    soil_bot, soil_top = find_soil_structure(surf_file)

    # -- Find surface dataset soil levels
    # TODO: how? NS uses metadata on file to find
    # soil strucure
    # better suggestion by WW to write dzsoi to neon surface dataset
    # This todo needs to go to the subset_data

    # TODO Will: if I sum them up , are they 3.5? (m) YES
    print("soil_top:", soil_top)
    print("soil_bot:", soil_bot)
    print("Sum of soil top depths    :", sum(soil_top))
    print("Sum of soil bottom depths :", sum(soil_bot))

    soil_top = np.cumsum(soil_top)
    soil_bot = np.cumsum(soil_bot)
    soil_mid = 0.5 * (soil_bot - soil_top) + soil_top
    # print ("Cumulative sum of soil bottom depths :", sum(soil_bot))

    obs_top = df["biogeoTopDepth"] / 100
    obs_bot = df["biogeoBottomDepth"] / 100

    # -- Mapping surface dataset and neon soil levels
    bins = df["biogeoTopDepth"] / 100
    bin_index = np.digitize(soil_mid, bins) - 1

    """
    print ("================================")
    print ("  Neon data soil structure:     ")
    print ("================================")

    print ("------------","ground","------------")
    for i in range(len(obs_bot)):
        print ("layer",i)
        print ("-------------",
                "{0:.2f}".format(obs_bot[i]),
                "-------------")

    print ("================================")
    print ("Surface data soil structure:    ")
    print ("================================")

    print ("------------","ground","------------")
    for b in range(len(bin_index)):
        print ("layer",b)
        print ("-------------",
                "{0:.2f}".format(soil_bot[b]),
                "-------------")
    """

    # -- update fields with neon
    f2 = f1
    soil_levels = f2["PCT_CLAY"].size
    for soil_lev in range(soil_levels):
        print("--------------------------")
        print("soil_lev:", soil_lev)
        print(df["clayTotal"][bin_index[soil_lev]])
        f2["PCT_CLAY"][soil_lev] = df["clayTotal"][bin_index[soil_lev]]
        f2["PCT_SAND"][soil_lev] = df["sandTotal"][bin_index[soil_lev]]

        bulk_den = df["bulkDensExclCoarseFrag"][bin_index[soil_lev]]
        carbon_tot = df["carbonTot"][bin_index[soil_lev]]
        estimated_oc = df["estimatedOC"][bin_index[soil_lev]]

        # -- estimated_oc in neon data is rounded to the nearest integer.
        # -- Check to make sure the rounded oc is not higher than carbon_tot.
        # -- Use carbon_tot if estimated_oc is bigger than carbon_tot.

        if estimated_oc > carbon_tot:
            estimated_oc = carbon_tot

        layer_depth = (
            df["biogeoBottomDepth"][bin_index[soil_lev]]
            - df["biogeoTopDepth"][bin_index[soil_lev]]
        )

        # f2["ORGANIC"][soil_lev] = estimated_oc * bulk_den / 0.58

        # -- after adding caco3 by NEON:
        # -- if caco3 exists:
        # -- inorganic = caco3/100.0869*12.0107
        # -- organic = carbon_tot - inorganic
        # -- else:
        # -- oranigc = estimated_oc * bulk_den /0.58

        caco3 = df["caco3Conc"][bin_index[soil_lev]]
        inorganic = caco3 / 100.0869 * 12.0107
        print("inorganic:", inorganic)

        if not np.isnan(inorganic):
            actual_oc = carbon_tot - inorganic
        else:
            actual_oc = estimated_oc

        f2["ORGANIC"][soil_lev] = actual_oc * bulk_den / 0.58

        print("~~~~~~~~~~~~~~~~~~~~~~~~")
        print("inorganic:")
        print("~~~~~~~~~~~~~~~~~~~~~~~~")
        print(inorganic)
        print("~~~~~~~~~~~~~~~~~~~~~~~~")

        print("bin_index    : ", bin_index[soil_lev])
        print("layer_depth  : ", layer_depth)
        print("carbon_tot   : ", carbon_tot)
        print("estimated_oc : ", estimated_oc)
        print("bulk_den     : ", bulk_den)
        print("organic      :", f2["ORGANIC"][soil_lev].values)
        print("--------------------------")

    # -- Interpolate missing values
    method = "linear"
    fill_interpolate(f2, "PCT_CLAY", method)
    fill_interpolate(f2, "PCT_SAND", method)
    fill_interpolate(f2, "ORGANIC", method)

    # -- Update zbedrock if neon observation does not make it down to 2m depth
    rock_thresh = 2

    zb_flag = False

    if obs_bot.iloc[-1] < rock_thresh:
        print("zbedrock is updated.")
        f2["zbedrock"].values[:, :] = obs_bot.iloc[-1]
        zb_flag = True

    sort_print_soil_layers(obs_bot, soil_bot)

    # -- updates for ag sites : KONA and STER
    ag_sites = ["KONA", "STER"]
    if site_name in ag_sites:
        print("Updating PCT_NATVEG")
        print("Original : ", f2.PCT_NATVEG.values)
        f2.PCT_NATVEG.values = [[0.0]]
        print("Updated  : ", f2.PCT_NATVEG.values)

        print("Updating PCT_CROP")
        print("Original : ", f2.PCT_CROP.values)
        f2.PCT_CROP.values = [[100.0]]
        print("Updated  : ", f2.PCT_CROP.values)

        print("Updating PCT_NAT_PFT")
        #print (f2.PCT_NAT_PFT)
        print(f2.PCT_NAT_PFT.values[0])
        f2.PCT_NAT_PFT.values[0] = [[100.0]]
        print(f2.PCT_NAT_PFT[0].values)

    out_dir = args.out_dir

    # -- make out_dir if it does not exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # -- update time tag for the output file
    wfile = out_dir + update_time_tag(surf_file)

    # -- update netcdf metadata
    f2 = update_metadata(f2, surf_file, neon_file, zb_flag)

    print(f2.attrs)
    f2.to_netcdf(path=wfile, mode="w", format="NETCDF3_64BIT")

    print(
        "Successfully updated surface data file for neon site("
        + site_name
        + "):\n - "
        + wfile
    )


if __name__ == "__main__":
    main()
