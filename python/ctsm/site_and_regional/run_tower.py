#! /usr/bin/env python3

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This is a wrapper script for running CTSM simulation for one or more
tower (neon or plumber) sites.

This script is only for supported tower sites and we will develop a
more general code later.

This script first creates and builds a generic base case.
Next, it will clone the base_case for different tower sites and run
types to reduce the need to build ctsm everytime.

This script will do the following:
    1) Create a generic base case for cloning.
    2) Make the case for the specific neon or plumber site(s).
    3) Make changes to the case, for:
        a. AD spinup
        b. post-AD spinup
        c. transient
        #---------------
        d. SASU or Matrix spinup
    4) Build and submit the case.

-------------------------------------------------------------------
Instructions for running using conda python environments:

../../py_env_create
conda activate ctsm_py

-------------------------------------------------------------------
To see the available options:
    ./run_tower.py --help
-------------------------------------------------------------------
"""
# TODO (NS)
# - [ ] Case dependency and the ability to check case status
# - [ ] If Case dependency works we don't need finidat given explicilty for post-ad and transient.

# - [ ] "./bin/git-fleximod update" instead of using env variable
# - [ ] wget the fields available and run for those available

# - [ ] Matrix spin-up if (SASU) Eric merged it in
# - [ ] Make sure both AD and SASU are not on at the same time

# - [ ] Make sure CIME and other dependencies are checked out.


# Import libraries
import glob
import logging
import os
import sys
import pandas as pd

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm.path_utils import path_to_ctsm_root
from ctsm.download_utils import download_file
from ctsm.site_and_regional.tower_arg_parse import get_parser
from ctsm.site_and_regional.neon_site import NeonSite
from ctsm.site_and_regional.plumber_site import Plumber2Site

# pylint: disable=import-error, wildcard-import, wrong-import-order
from standard_script_setup import *

logger = logging.getLogger(__name__)


def check_neon_listing(valid_neon_sites):
    """
    A function to download and parse neon listing file.
    """
    listing_file = "listing.csv"
    url = "https://storage.neonscience.org/neon-ncar/listing.csv"

    download_file(url, listing_file)
    available_list = parse_neon_listing(listing_file, valid_neon_sites)
    return available_list


def parse_neon_listing(listing_file, valid_neon_sites):
    """
    A function to parse neon listing file
    and find neon sites with the dates
    where data is available.

    Args:
        listing_file (str): downloaded listing file

    Returns:
        available_list :
            list of neon_site objects that is found
            on the downloaded listing file.
    """

    # pd.set_option("display.max_rows", None, "display.max_columns", None)

    available_list = []

    listing_df = pd.read_csv(listing_file)

    # check for finidat files for transient run
    finidatlist = listing_df[listing_df["object"].str.contains("lnd/ctsm")]

    # -- filter lines with atm/cdep
    listing_df = listing_df[listing_df["object"].str.contains("atm/cdeps/")]

    # -- split the object str to extract site name
    listing_df = listing_df["object"].str.split("/", expand=True)

    # -- groupby site name
    grouped_df = listing_df.groupby(8)
    for key, _ in grouped_df:
        # -- check if it is a valid neon site
        if any(key in x for x in valid_neon_sites):
            site_name = key
            tmp_df = grouped_df.get_group(key)

            # -- filter files only ending with YYYY-MM.nc
            tmp_df = tmp_df[tmp_df[9].str.contains(r"\d\d\d\d-\d\d.nc")]

            # -- find all the data versions
            # versions = tmp_df[7].unique()
            # print ("all versions available for ", site_name,":", *versions)
            latest_version = tmp_df[7].iloc[-1]
            # print ("latest version available for ", site_name,":", latest_version)

            tmp_df = tmp_df[tmp_df[7].str.contains(latest_version)]
            # -- remove .nc from the file names
            tmp_df[9] = tmp_df[9].str.replace(".nc", "", regex=False)

            tmp_df2 = tmp_df[9].str.split("-", expand=True)

            # ignore any prefix in file name and just get year
            tmp_df2[0] = tmp_df2[0].str.slice(-4)

            # -- figure out start_year and end_year
            start_year = tmp_df2[0].iloc[0]
            end_year = tmp_df2[0].iloc[-1]

            # -- figure out start_month and end_month
            start_month = tmp_df2[1].iloc[0]
            end_month = tmp_df2[1].iloc[-1]

            logger.debug("Valid neon site %s found!", site_name)
            logger.debug("File version %s", latest_version)
            logger.debug("start_year=%s", start_year)
            logger.debug("end_year=%s", end_year)
            logger.debug("start_month=%s", start_month)
            logger.debug("end_month=%s", end_month)
            finidat = None
            for line in finidatlist["object"]:
                if site_name in line:
                    finidat = line.split(",")[0].split("/")[-1]

            neon_site = NeonSite(
                name=site_name,
                start_year=start_year,
                end_year=end_year,
                start_month=start_month,
                end_month=end_month,
                finidat=finidat,
            )
            logger.debug(neon_site)
            available_list.append(neon_site)

    return available_list


def setup_plumber_data(valid_plumber_sites):
    """
    A function to set up plumber site objects
    with dummy start and end years and months.
    This allows us to use the list of plumber site objects.

    Returns:
        available_list :
            list of plumber_site objects that is found
    """

    available_list = []

    for site_name in valid_plumber_sites:

        # start_year and end_year are set in shell commands, so these get overwritten
        start_year = 8888
        end_year = 9999
        start_month = 1
        end_month = 12

        logger.debug("Valid plumber site %s found!", site_name)
        finidat = None

        plumber_site = Plumber2Site(
            name=site_name,
            start_year=start_year,
            end_year=end_year,
            start_month=start_month,
            end_month=end_month,
            finidat=finidat,
        )
        available_list.append(plumber_site)

    return available_list


def main(description):
    """
    Determine valid tower sites. Make an output directory if it does not exist.
    Loop through requested sites and run CTSM at that site.
    """
    cesmroot = path_to_ctsm_root()
    # Get the list of supported neon sites from usermods
    # The [!Fd]* portion means that we won't retrieve cases that start with:
    # F (FATES) or d (default). We should be aware of adding cases that start with these.
    valid_neon_sites = glob.glob(
        os.path.join(cesmroot, "cime_config", "usermods_dirs", "clm", "NEON", "[!Fd]*")
    )
    valid_neon_sites = sorted([v.split("/")[-1] for v in valid_neon_sites])

    # Get the list of supported plumber sites from usermods
    valid_plumber_sites = glob.glob(
        os.path.join(cesmroot, "cime_config", "usermods_dirs", "clm", "PLUMBER2", "[!d]*")
    )

    valid_plumber_sites = sorted([v.split("/")[-1] for v in valid_plumber_sites])

    (
        neon_site_list,
        plumber_site_list,
        output_root,
        run_type,
        experiment,
        prism,
        overwrite,
        base_case_root,
        run_from_postad,
        setup_only,
        no_batch,
        rerun,
        no_input_data_check,
        user_version,
        xmlchange,
    ) = get_parser(sys.argv, description, valid_neon_sites, valid_plumber_sites)

    if output_root:
        logger.debug("output_root : %s", output_root)
        if not os.path.exists(output_root):
            os.makedirs(output_root)

    # -- check neon listing file for available data:
    available_list = check_neon_listing(valid_neon_sites)

    # =================================
    # -- all tower sites can be cloned from one generic case
    # -- so no need to define a base_case for every site.

    res = "CLM_USRDAT"
    if run_type == "transient":
        compset = "IHist1PtClm60Bgc"
    else:
        compset = "I1PtClm60Bgc"

    # --  Looping over neon sites
    if neon_site_list:
        for neon_site in available_list:
            if neon_site.name in neon_site_list:
                if run_from_postad:
                    neon_site.finidat = None
                if not base_case_root:
                    neon_site.set_default_user_mods_dirs()
                    base_case_root = neon_site.build_base_case(
                        cesmroot,
                        output_root,
                        res,
                        compset,
                        overwrite=overwrite,
                        setup_only=setup_only,
                    )
                logger.info("-----------------------------------")
                logger.info("Running CTSM for neon site : %s", neon_site.name)

                neon_site.run_case(
                    base_case_root=base_case_root,
                    run_type=run_type,
                    prism=prism,
                    user_version=user_version,
                    overwrite=overwrite,
                    setup_only=setup_only,
                    no_batch=no_batch,
                    rerun=rerun,
                    experiment=experiment,
                    no_input_data_check=no_input_data_check,
                    xmlchange=xmlchange,
                )

    # -- check for available plumber data:
    available_plumber_list = setup_plumber_data(valid_plumber_sites)

    # --  Looping over plumber sites
    if plumber_site_list:
        for plumber_site in available_plumber_list:
            if plumber_site.name in plumber_site_list:
                if run_from_postad:
                    plumber_site.finidat = None
                if not base_case_root:
                    plumber_site.set_default_user_mods_dirs()
                    base_case_root = plumber_site.build_base_case(
                        cesmroot,
                        output_root,
                        res,
                        compset,
                        overwrite=overwrite,
                        setup_only=setup_only,
                    )
                logger.info("-----------------------------------")
                logger.info("Running CTSM for plumber site : %s", plumber_site.name)
                plumber_site.run_case(
                    base_case_root=base_case_root,
                    run_type=run_type,
                    prism=prism,
                    user_version=user_version,
                    overwrite=overwrite,
                    setup_only=setup_only,
                    no_batch=no_batch,
                    rerun=rerun,
                    experiment=experiment,
                    no_input_data_check=no_input_data_check,
                    xmlchange=xmlchange,
                )
