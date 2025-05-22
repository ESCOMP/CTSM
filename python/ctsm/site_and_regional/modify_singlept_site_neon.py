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
    in the specified directory : (i.e. ../../../neon_surf_files).
- Modify surface dataset with downloaded data.

-------------------------------------------------------------------
Instructions for running using conda python environments:

../../py_env_create
conda activate ctsm_py

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
# --[] Download files only when available.

#  Import libraries
from __future__ import print_function

import argparse
from datetime import date
from getpass import getuser
import glob
import logging
import os
import sys
import requests

import numpy as np
import pandas as pd
import xarray as xr
from packaging import version

from ctsm.path_utils import path_to_ctsm_root

myname = getuser()

# Seconds to wait before requests.get() times out
TIMEOUT = 60


# -- valid neon sites
valid = glob.glob(
    os.path.join(path_to_ctsm_root(), "cime_config", "usermods_dirs", "clm", "NEON", "[!d]*")
)
valid_neon_sites = [x[-4:] for x in valid]  # last 4 letters in each string


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
        default="/glade/derecho/scratch/" + myname + "/single_point/",
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
        default="/glade/derecho/scratch/" + myname + "/single_point_neon_updated/",
    )
    parser.add_argument(
        "--inputdata-dir",
        help="""
                Directory containing standard input files from CESM input data such as the surf_soildepth_file.
                [default: %(default)s]
                """,
        action="store",
        dest="inputdatadir",
        type=str,
        required=False,
        default="/glade/campaign/cesm/cesmdata/cseg/inputdata",
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

    parser.add_argument(
        "--16pft",
        help="Modify 16-pft surface data files (e.g. for a FATES run)",
        action="store_true",
        dest="pft_16",
        default=False,
    )

    return parser


def get_neon(neon_dir, site_name):
    """
    Function for finding neon data files
    and download from neon server if the
    file does not exist.

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

    # -- Download the file if it does not exist
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
        response = requests.get(url, timeout=TIMEOUT)

        with open(neon_file, "wb") as a_file:
            a_file.write(response.content)

        # -- Check if download status_code
        if response.status_code == 200:
            print("Download finished successfully for", site_name)
        elif response.status_code == 404:
            sys.exit(
                "Data for this site " + site_name + " was not available on the neon server:" + url
            )

        print("Download exit status code:  ", response.status_code)
        print("Downloaded file type     :  ", response.headers["content-type"])
        print("Downloaded file encoding :  ", response.encoding)
        print("------------------------------------------------")

        response.close()

    return neon_file


def find_surffile(surf_dir, site_name, pft_16):
    """
    Function for finding and choosing surface file for
    a neon site.
    These files are created using ./subset_data.py script.
    In case multiple files exist for the neon site, it
    will choose the file created the latest.

    Args:
        surf_dir (str): directory of single point surface data
        site_name (str): 4 letter neon site name
        pft_16 (bool):    if true, use 16-PFT version of surface data file

    Raises:
        Error if the surface data for the site is not created

    Returns:
        surf_file (str): name of the surface dataset file
    """

    if pft_16:
        sf_name = "surfdata_1x1_NEON_" + site_name + "*hist_2000_16pfts*.nc"
    else:
        sf_name = "surfdata_1x1_NEON_" + site_name + "*hist_2000_78pfts*.nc"

    print(os.path.join(surf_dir, sf_name))
    surf_file = sorted(glob.glob(os.path.join(surf_dir, sf_name)))

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
            "Surface data for this site "
            + str(site_name)
            + " was not found:"
            + str(surf_dir)
            + str(sf_name)
            + "."
            + "\n"
            + "Please run ./subset_data.py for this site."
        )
    return surf_file


def find_soil_structure(args, surf_file):
    """
    Function for finding surface dataset soil
    structure using surface data metadata.

    In CLM surface data, soil layer information
    is in a file from surface data metadata
    under "Soil_texture_raw_data_file_name".
    This function finds this file for the surface
    dataset, read it, and find soil layers.

    args:
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
    f_1 = xr.open_dataset(surf_file)
    print("------------")
    # print (f_1.attrs["Soil_texture_raw_data_file_name"])

    clm_input_dir = os.path.join(args.inputdatadir, "lnd/clm2/rawdata/")
    surf_soildepth_file = os.path.join(
        clm_input_dir, f_1.attrs["soil_texture_lookup_raw_data_file_name"]
    )

    if os.path.exists(surf_soildepth_file):
        print(
            "\n\n Reading",
            surf_soildepth_file,
            "for surface data soil structure information:",
        )
        f_1_soildepth = xr.open_dataset(surf_soildepth_file)
        print(f_1_soildepth["DZSOI"])
        soil_bot = f_1_soildepth["DZSOI"].values

        # -- soil layer top
        soil_top = soil_bot[:-1]
        soil_top = np.insert(soil_top, 0, 0)

    else:
        sys.exit(
            "Cannot find soil structure file : " + surf_soildepth_file + "for the surface dataset."
        )

    return soil_bot, soil_top


def update_metadata(nc_file, surf_file, neon_file, zb_flag):
    """
    Function for updating modified surface dataset
    metadata for neon sites.

    Args:
        nc_file (xr Dataset): netcdf file including updated neon surface data
        surf_file (str): single point surface data filename
        neon_file (str): filename of neon downloaded surface dataset
        zb_flag (bool): update bedrock

    Returns:
        nc_file (xr Dataset): netcdf file including updated neon surface data
    """
    today = date.today()
    today_string = today.strftime("%Y-%m-%d")

    nc_file.attrs["Updated_on"] = today_string
    nc_file.attrs["Updated_by"] = myname
    nc_file.attrs["Updated_with"] = os.path.abspath(__file__)
    nc_file.attrs["Updated_from"] = surf_file
    nc_file.attrs["Updated_using"] = neon_file
    if zb_flag:
        nc_file.attrs["Updated_fields"] = "PCT_CLAY, PCT_SAND, ORGANIC, zbedrock"
    else:
        nc_file.attrs["Updated_fields"] = "PCT_CLAY, PCT_SAND, ORGANIC"

    return nc_file


def update_time_tag(fname_in):
    """
    Function for updating time tag on surface dataset
    files.
    Expects file to end with [._]cYYMMDD.nc or [._]YYMMDD.nc
    Add the tag to just before that ending part.

    Args:
        fname_in (str) : file name with the old time tag

    Raises:
        error if the file does not end with
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

    for _, row in depth_df.iterrows():
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

    d_f = pd.read_csv(listing_file)
    d_f = d_f[d_f["object"].str.contains("_surfaceData.csv")]
    dict_out = dict(zip(d_f["object"], d_f["last_modified"]))
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
        response = requests.get(url, timeout=TIMEOUT)

        with open(fname, "wb") as a_file:
            a_file.write(response.content)

        # -- Check if download status_code
        if response.status_code == 200:
            print("Download finished successfully for", fname, ".")
        elif response.status_code == 404:
            print("File " + fname + "was not available on the neon server:" + url)
    except Exception as err:
        print("The server could not fulfill the request.")
        print("Something went wrong in downloading", fname)
        raise err


def fill_interpolate(f_2, var, method):
    """
    Function to interpolate a variable in a
    xarray dataset a specific method
    """
    print("=====================================")
    print("Filling in ", var, "with interpolation (method =" + method + ").")

    print("Variable before filling : ")
    print(f_2[var])

    tmp_df = pd.DataFrame(f_2[var].values.ravel())

    tmp_df = tmp_df.interpolate(method=method, limit_direction="both")

    tmp = tmp_df.to_numpy()

    soil_levels = f_2[var].size
    for soil_lev in range(soil_levels):
        f_2[var][soil_lev] = tmp[soil_lev].reshape(1, 1)

    print("Variable after filling : ")
    print(f_2[var])
    print("=====================================")


def print_neon_data_soil_structure(obs_bot, soil_bot, bin_index):
    """
    Print info about NEON data soil structure
    """
    print("================================")
    print("  Neon data soil structure:     ")
    print("================================")

    print("------------", "ground", "------------")
    for i, this_obs_bot in enumerate(obs_bot):
        print("layer", i)
        print("-------------", "{0:.2f}".format(this_obs_bot), "-------------")

    print("================================")
    print("Surface data soil structure:    ")
    print("================================")

    print("------------", "ground", "------------")
    for this_bin in range(len(bin_index)):
        print("layer", this_bin)
        print("-------------", "{0:.2f}".format(soil_bot[this_bin]), "-------------")


def print_soil_quality(
    *, inorganic, bin_index, soil_lev, layer_depth, carbon_tot, estimated_oc, bulk_den, f_2
):
    """
    Prints information about soil quality
    """
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
    print("organic      :", f_2["ORGANIC"][soil_lev].values)
    print("--------------------------")


def update_agri_site_info(site_name, f_2):
    """
    Updates agricultural sites
    """
    ag_sites = ["KONA", "STER"]
    if site_name not in ag_sites:
        return f_2

    print("Updating PCT_NATVEG")
    print("Original : ", f_2.PCT_NATVEG.values)
    f_2.PCT_NATVEG.values = [[0.0]]
    print("Updated  : ", f_2.PCT_NATVEG.values)

    print("Updating PCT_CROP")
    print("Original : ", f_2.PCT_CROP.values)
    f_2.PCT_CROP.values = [[100.0]]
    print("Updated  : ", f_2.PCT_CROP.values)

    print("Updating PCT_NAT_PFT")
    print(f_2.PCT_NAT_PFT.values[0])
    print(f_2.PCT_NAT_PFT[0].values)

    return f_2


def update_fields_with_neon(f_1, d_f, bin_index):
    """
    update fields with neon
    """
    f_2 = f_1
    soil_levels = f_2["PCT_CLAY"].size
    for soil_lev in range(soil_levels):
        print("--------------------------")
        print("soil_lev:", soil_lev)
        print(d_f["clayTotal"][bin_index[soil_lev]])
        f_2["PCT_CLAY"][soil_lev] = d_f["clayTotal"][bin_index[soil_lev]]
        f_2["PCT_SAND"][soil_lev] = d_f["sandTotal"][bin_index[soil_lev]]

        bulk_den = d_f["bulkDensExclCoarseFrag"][bin_index[soil_lev]]
        carbon_tot = d_f["carbonTot"][bin_index[soil_lev]]
        estimated_oc = d_f["estimatedOC"][bin_index[soil_lev]]

        # -- estimated_oc in neon data is rounded to the nearest integer.
        # -- Check to make sure the rounded oc is not higher than carbon_tot.
        # -- Use carbon_tot if estimated_oc is bigger than carbon_tot.

        estimated_oc = min(estimated_oc, carbon_tot)

        layer_depth = (
            d_f["biogeoBottomDepth"][bin_index[soil_lev]]
            - d_f["biogeoTopDepth"][bin_index[soil_lev]]
        )

        # f_2["ORGANIC"][soil_lev] = estimated_oc * bulk_den / 0.58

        # -- after adding caco3 by NEON:
        # -- if caco3 exists:
        # -- inorganic = caco3/100.0869*12.0107
        # -- organic = carbon_tot - inorganic
        # -- else:
        # -- organic = estimated_oc * bulk_den /0.58

        caco3 = d_f["caco3Conc"][bin_index[soil_lev]]
        inorganic = caco3 / 100.0869 * 12.0107
        print("inorganic:", inorganic)

        if not np.isnan(inorganic):
            actual_oc = carbon_tot - inorganic
        else:
            actual_oc = estimated_oc

        f_2["ORGANIC"][soil_lev] = actual_oc * bulk_den / 0.58

        print_soil_quality(
            inorganic=inorganic,
            bin_index=bin_index,
            soil_lev=soil_lev,
            layer_depth=layer_depth,
            carbon_tot=carbon_tot,
            estimated_oc=estimated_oc,
            bulk_den=bulk_den,
            f_2=f_2,
        )
    return f_2


def main():
    """modify_singlept_site_neon main function"""
    args = get_parser().parse_args()

    # -- debugging option
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    # Check if pandas is a recent enough version
    pdvers = pd.__version__
    if version.parse(pdvers) < version.parse("1.1.0"):
        sys.exit(
            """The pandas version in your python environment is too old,
            update to a newer version of pandas (>=1.1.0): version=%s""",
            pdvers,
        )

    # file_time = check_neon_time()

    # --  specify site from which to extract data
    site_name = args.site_name

    # --  Look for surface data
    surf_dir = args.surf_dir
    surf_file = find_surffile(surf_dir, site_name, args.pft_16)

    # --  directory structure
    clone_dir = os.path.abspath(os.path.join(__file__, "../../../.."))
    neon_dir = os.path.join(clone_dir, "neon_surffiles")

    # --  download neon data if needed
    neon_file = get_neon(neon_dir, site_name)

    # -- Read neon data
    d_f = pd.read_csv(neon_file)

    # -- Read surface dataset files
    print("surf_file:", surf_file)
    f_1 = xr.open_dataset(surf_file)

    # -- Find surface dataset soil depth information
    soil_bot, soil_top = find_soil_structure(args, surf_file)

    # -- Find surface dataset soil levels
    # TODO: how? NS uses metadata on file to find
    # soil strucure
    # better suggestion by WW to write dzsoi to neon surface dataset
    # This todo needs to go to the subset_data

    soil_top = np.cumsum(soil_top)
    soil_bot = np.cumsum(soil_bot)
    soil_mid = 0.5 * (soil_bot - soil_top) + soil_top
    # print ("Cumulative sum of soil bottom depths :", sum(soil_bot))

    obs_bot = d_f["biogeoBottomDepth"] / 100

    # -- Mapping surface dataset and neon soil levels
    bins = d_f["biogeoTopDepth"] / 100
    bin_index = np.digitize(soil_mid, bins) - 1

    print_neon_data_soil_structure(obs_bot, soil_bot, bin_index)

    # -- update fields with neon
    f_2 = update_fields_with_neon(f_1, d_f, bin_index)

    # -- Interpolate missing values
    method = "linear"
    fill_interpolate(f_2, "PCT_CLAY", method)
    fill_interpolate(f_2, "PCT_SAND", method)
    fill_interpolate(f_2, "ORGANIC", method)

    # -- Update zbedrock if neon observation does not make it down to 2m depth
    rock_thresh = 2

    zb_flag = False

    if obs_bot.iloc[-1] < rock_thresh:
        print("zbedrock is updated.")
        f_2["zbedrock"].values[:, :] = obs_bot.iloc[-1]
        zb_flag = True

    sort_print_soil_layers(obs_bot, soil_bot)

    # -- updates for ag sites
    update_agri_site_info(site_name, f_2)

    out_dir = args.out_dir

    # -- make out_dir if it does not exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # -- update time tag for the output file
    wfile = out_dir + update_time_tag(surf_file)

    # -- update netcdf metadata
    f_2 = update_metadata(f_2, surf_file, neon_file, zb_flag)

    print(f_2.attrs)
    f_2.to_netcdf(path=wfile, mode="w", format="NETCDF3_64BIT")

    print("Successfully updated surface data file for neon site(" + site_name + "):\n - " + wfile)
