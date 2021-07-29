#! /usr/bin/env python3
 
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This is a wrapper script for running CTSM simulation for one or more
neon sites. 

This script is only for neon site and we will develop a more general
code later.

This script first creates and builds a generic base case. 
Next, it will clone the base_case for different neon sites and run
types to reduce the need to build ctsm everytime. 

This script will do the following:
    1) Create a generic base case for cloning. 
    2) Make the case for the specific neon site(s).
    3) Make changes to the case, for:
        a. AD spinup
	    b. post-AD spinup
        c. transient
    	#---------------
    	d. SASU or Matrix spinup
    4) Build and submit the case.
 
-------------------------------------------------------------------
Instructions for running on Cheyenne/Casper:
 
load the following into your local environment
    module load python
    ncar_pylib
 
To remove NPL from your environment on Cheyenne/Casper:
    deactivate
 
-------------------------------------------------------------------
To see the available options:
    ./run_neon.py --help
-------------------------------------------------------------------
"""
#TODO (NS)
#- [ ]
#- [ ] Case dependency and the ability to check case status
#- [ ] If Case dependency works we don't need finidat given explicilty for post-ad and transient.

#- [ ] checkout_externals instead of using env varaiable
#- [ ] wget the fields available and run for those available
 
#- [ ] Matrix spin-up if (SASU) Eric merged it in 
#- [ ] Make sure both AD and SASU are not on at the same time

#- [ ] Make sure CIME and other dependencies is checked out.

 
#Import libraries

import os
import sys
import time 
import shutil
import logging
import requests
import argparse
import subprocess
import pandas as pd
import glob 
from datetime import date
from getpass import getuser
 
#-- get the environment variable
# Then something like this:
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..","..","..",'python'))
print("{}".format(_CTSM_PYTHON))
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root

import CIME.build as build
from standard_script_setup import *
from CIME.case             import Case
from CIME.utils            import safe_copy
from argparse              import RawTextHelpFormatter
from CIME.locked_files     import lock_file, unlock_file

logger = logging.getLogger(__name__) 
 
def get_parser(args, description, valid_neon_sites):                                                                                                                                           
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(description=description,
                           formatter_class=argparse.RawDescriptionHelpFormatter)

    CIME.utils.setup_standard_logging_options(parser)
    
    parser.print_usage = parser.print_help

    parser.add_argument('--neon-sites',
                help='4-letter neon site code.', 
                action="store",
                required=False,
                choices=valid_neon_sites + ['all'],
                dest="neon_sites",
                default=["OSBS"],
                nargs='+')

# not used
#    parser.add_argument('--surf-dir',
#                help='''
#                Directory of single point surface dataset.
#                [default: %(default)s]
#                ''', 
#                action="store", 
#                dest="surf_dir",
#                type =str,
#                required=False,
#                default="/glade/scratch/"+myname+"/single_point/")

    parser.add_argument('--case-root',
                help='''
                Root Directory of cases
                [default: %(default)s] 
                ''',
                action="store", 
                dest="case_root",
                type =str,
                required=False,
                default="CIME_OUTPUT_ROOT as defined in cime")

    subparsers = parser.add_subparsers (
                        dest='run_type',
                        help='Four different ways to run this script.')

    ad_parser = subparsers.add_parser ('ad',
                help=''' AD spin-up options ''') 

    pad_parser = subparsers.add_parser ('postad',
                help=''' Post-AD spin-up options ''')

    tr_parser = subparsers.add_parser ('transient',
                help=''' Transient spin-up options ''')

    sasu_parser = subparsers.add_parser ('sasu',
                help=''' Sasu spin-up options --not in CTSM yet''')

    ad_parser.add_argument ('--ad-length',
                help='''
                How many years to run AD spin-up
                [default: %(default)s]
                ''', 
                required = False,
                type = int, 
                default = 200)

    pad_parser.add_argument ('--postad-length',
                help='''
                How many years to run in post-AD mode
                [default: %(default)s]
                ''', 
                required = False,
                type = int, 
                default = 100)

    pad_parser.add_argument('--finidat',
                help='''
                finidat file location from spinup step to start from.
                [default: %(default)s]
                ''', 
                action="store",
                dest="finidat_postad",
                required = True,
                type = str)

    tr_parser.add_argument('--start-year',
                help='''           
                Start year for running CTSM simulation.
                [default: %(default)s]
                ''',               
                action="store", 
                dest="start_year",                                                                                                                                                 
                required = False,
                type = int,        
                default = 2018)

    tr_parser.add_argument('--end-year',
                help='''
                End year for running CTSM simulation.
                [default: %(default)s]
                ''', 
                action="store",
                dest="end_year",
                required = False,
                type = int,
                default = 2020)

    tr_parser.add_argument('--finidat',
                help='''
                finidat file location from spinup step to start from.
                [default: %(default)s]
                ''', 
                action="store",
                dest="finidat_transient",
                required = True,
                type = str)

    parser.add_argument('--overwrite',
                help='''
                overwrite existing case directories
                [default: %(default)s]
                ''', 
                action="store_true",
                dest="overwrite",
                required = False,
                default = False)

    #parser.add_argument('--spinup',
    #            help='''
    #            AD spin-up
    #            [default: %(default)s]
    #            ''', 
    #            action="store_true",
    #            dest="ad_flag",
    #            required = False,
    #            default = True)
    #parser.add_argument('--postad',
    #            help='''
    #            Post-AD spin-up
    #            [default: %(default)s]
    #            ''', 
    #            action="store_true",
    #            dest="postad_flag",
    #            required = False,
    #            default = True)
    #parser.add_argument('--transient',
    #            help='''
    #            Transient
    #            [default: %(default)s]
    #            ''', 
    #            action="store_true",
    #            dest="transient_flag",
    #            required = False,
    #            default = True)

    #parser.add_argument('--sasu','--matrix',
    #            help='''
    #            Matrix (SASU) spin-up
    #            [default: %(default)s]
    #            ''', 
    #            action="store_true",
    #            dest="sasu_flag",
    #            required = False,
    #            default = False)

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args, parser)

    if 'all' in args.neon_sites:
        neon_sites = valid_neon_sites
    else:
        neon_sites = args.neon_sites
        for site in neon_sites:
            if site not in valid_neon_sites:
                raise ValueError("Invalid site name {}".format(site))

    if "CIME_OUTPUT_ROOT" in args.case_root:
        args.case_root = None
            
    return neon_sites, args.case_root, args.run_type, args.overwrite

class NeonSite :
    """
    A class for encapsulating neon sites.
       
    ...
       
    Attributes
    ----------
       
    Methods
    -------
    """
    def __init__(self, name, start_year, end_year, start_month, end_month):
        self.name = name
        self.start_year= start_year
        self.end_year = end_year
        self.start_month = start_month
        self.end_month = end_month
        
    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' 
                    for item in (self.__dict__)))

    def build_base_case(self, cesmroot, case_root, res, compset, overwrite):
        """
        Function for building a base_case to clone.
        To spend less time on building ctsm for the neon cases,
        all the other cases are cloned from this case
        
        Args:
        self: 
            The NeonSite object
        base_root (str): 
            root of the base_case CIME 
        res (str):
            base_case resolution or gridname
        compset (str):
            base case compset
        overwrite (bool) : 
            Flag to overwrite the case if exists
        """
        logger.info("---- building a base case -------")
        
        user_mods_dirs = [os.path.join(cesmroot,"cime_config","usermods_dirs","NEON",self.name)]
        if case_root:
            case_path = os.path.join(case_root,self.name)
            logger.info ('case_root      : {}'.format(case_root))
        else:
            case_path = self.name
            
        logger.info ('base_case_name : {}'.format(self.name))
        logger.info ('user_mods_dir  : {}'.format(user_mods_dirs[0]))

        if overwrite and os.path.isdir(case_path):
            logger.info ("Removing the existing case at: ", case_path)
            shutil.rmtree(case_path)
            
        with Case(case_path, read_only=False) as case:
            if not os.path.isdir(case_path):
                logger.info("---- creating a base case -------")

                case.create(case_path, cesmroot, compset, res, mpilib="mpi-serial",
                            run_unsupported=True, answer="r",walltime="04:00:00",
                            user_mods_dirs = user_mods_dirs, driver="nuopc")

                logger.info("---- base case created ------")

                #--change any config for base_case:
                #case.set_value("RUN_TYPE","startup")


                logger.info("---- base case setup ------")
                case.case_setup()
            else:
                case.case_setup(reset=True)

            logger.info("---- base case build ------")
            # always walk through the build process to make sure it's up to date. 
            t0 = time.time()
            build.case_build(case_path, case=case)
            t1 = time.time()
            total = t1-t0
            logger.info ("Time required to building the base case: {} s.".format(total))
            # update case_path to be the full path to the base case
            case_path = case.get_value("CASEROOT")
        return case_path


def check_neon_listing(valid_neon_sites):
    """
    A function to download and parse neon listing file.
    """
    listing_file = 'listing.csv'
    url = 'https://neon-ncar.s3.data.neonscience.org/listing.csv'
 
    download_file(url, listing_file)
    available_list= parse_neon_listing(listing_file, valid_neon_sites)
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
 
    #pd.set_option("display.max_rows", None, "display.max_columns", None)
 
    available_list = []
 
    df = pd.read_csv(listing_file)
 
    #-- TODO: do we want to check for v2 in future?
 
    #-- filter lines with atm/cdep/v1
    df = df[df['object'].str.contains("atm/cdeps/v1")]
 
    #-- split the object str to extract site name
    df=df['object'].str.split("/", expand=True)
 
    #-- groupby site name
    grouped_df = df.groupby(7)
 
    for key, item in grouped_df:
        #-- check if it is a valid neon site
        if any(key in x for x in valid_neon_sites):
            site_name = key
            logger.info ("Valid neon site " + site_name+" found!")

            tmp_df = grouped_df.get_group(key)
 
            #-- filter files only ending with .nc
            tmp_df = tmp_df[tmp_df[8].str.contains('.nc')]
 
            #-- remove .nc from the file names
            tmp_df[8] = tmp_df[8].str.replace('.nc','')
 
            tmp_df2 = tmp_df[8].str.split("-",  expand=True)

            #-- figure out start_year and end_year
            start_year = tmp_df2[0].iloc[0]
            end_year = tmp_df2[0].iloc[-1]
 
            #-- figure out start_month and end_month
            start_month = tmp_df2[1].iloc[0]  
            end_month = tmp_df2[1].iloc[-1] 

            logger.debug ('start_year='+start_year)
            logger.debug ('start_year=',end_year)
            logger.debug ('start_month=',start_month)
            logger.debug ('end_month=',end_month)
 
            neon_site = NeonSite(site_name, start_year, end_year, start_month, end_month)
            logger.debug (neon_site)
            available_list.append(neon_site)
 
    return available_list

def download_file(url, fname):
    """
    Function to download a file.
 
    Args:
        url (str): 
            url of the file for downloading
 
        fname (str) : 
            file name to save the downloaded file.
    """
    response = requests.get(url)
 
    with open(fname, 'wb') as f:
        f.write(response.content)
 
    #-- Check if download status_code
    if response.status_code == 200:
        print('Download finished successfully for', fname,'.')
    elif response.status_code == 404:
        print('File '+fname+'was not available on the neon server:'+ url)



#-----------------------------------
# NS thinks these should be inside NeonSite object.
#-----------------------------------
def run_spinup_ad(orig_root, case_root, user_mods_dir, overwrite, ad_length, neon_site):
    """
    Function to run accelerated decomposition over specific length.
    Accelerated Decomposition is denoted as AD in the code.
    For creating the ad case we clone the base case
 
    Args:
        orig_root (str): 
            root of the base_case
            that we want to clone from
        case_root (str):
            Name of the AD case root 
        user_mods_dir (list):
            list of user_mods_dir
            Note, cime does not accept str only
            and it should be a list of strs.
        overwrite (bool) : 
            Flag to overwrite the case if exists
        ad_length (int) :
            Length of the AD spinup in years
        neon_site (NeonSite):
            NeonSite object that we are running.
    """

    if not os.path.isdir(orig_root):
        sys.exit('Base case does not exist in', orig_root)

    if overwrite and os.path.isdir(case_root):
        shutil.rmtree(case_root)
        logger.info("---- removing the existing case -------")

    if not os.path.isdir(case_root):
        with Case(orig_root, read_only=False) as clone:
            logger.info("---- cloning the base case -------")
            clone.create_clone(case_root, keepexe=True)

    with Case(case_root, read_only=False) as case:
        case.set_value("DATM_YR_ALIGN",neon_site.start_year)
        case.set_value("DATM_YR_START",neon_site.start_year)
        case.set_value("DATM_YR_END",neon_site.end_year)
        case.set_value("STOP_OPTION", "nyears")

        case.set_value("CLM_FORCE_COLDSTART","on")
        case.set_value("CLM_ACCELERATED_SPINUP","on")

        case.set_value("STOP_N", ad_length.__str__())
        case.set_value("REST_N", "100")

        case.set_value("RUN_REFDATE", "0018-01-01")
        case.set_value("RUN_STARTDATE", "0018-01-01")

        case.set_value("CONTINUE_RUN", "FALSE")
        case.set_value("RESUBMIT", "1")

        user_nl_fname = os.path.join(case_root, "user_nl_clm")

        user_nl_file  = open(user_nl_fname, "a")
        user_nl_lines = [
                        "hist_mfilt = 20",
                        "hist_nhtfrq = -8760",
                        "hist_empty_htapes = .true.",
                        "hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'"]

        for line in user_nl_lines:
            user_nl_file.write("%s\n" % line)

        user_nl_file.close()

        case.create_namelists()
        case.submit()


def run_postad(orig_root, case_root, user_mods_dir, overwrite, postad_length, neon_site):
    """
    Function to run post-ad simultion for a proper spin-up over specific length.
    For creating the ad case we clone the base case
 
    Args:
        orig_root (str): 
            root of the base_case
            that we want to clone from
        case_root (str):
            Name of the AD case root 
        user_mods_dir (list):
            list of user_mods_dir
            Note, cime does not accept str only
            and it should be a list of strs.
        overwrite (bool) : 
            Flag to overwrite the case if exists
        postad_length (int) :
            Length of the AD spinup in years
        neon_site (NeonSite):
            NeonSite object that we are running.
    """

    if not os.path.isdir(orig_root):
        sys.exit('Base case does not exist in', orig_root)

    if overwrite and os.path.isdir(case_root):
        shutil.rmtree(case_root)
        logger.info("---- removing the existing case -------")

    if not os.path.isdir(case_root):
        with Case(orig_root, read_only=False) as clone:
            logger.info("---- cloning the base case -------")
            clone.create_clone(case_root, keepexe=True)

    with Case(case_root, read_only=False) as case:
        case.set_value("DATM_YR_ALIGN",neon_site.start_year)
        case.set_value("DATM_YR_START",neon_site.start_year)
        case.set_value("DATM_YR_END",neon_site.end_year)
        case.set_value("STOP_OPTION", "nyears")

        case.set_value("CLM_FORCE_COLDSTART","off")
        case.set_value("CLM_ACCELERATED_SPINUP","off")

        case.set_value("STOP_N", postad_length.__str__())
        case.set_value("REST_N", "100")

        #TODO: change this instead of hardcoding it to 218
        case.set_value("RUN_REFDATE", "0218-01-01")
        case.set_value("RUN_STARTDATE", "0218-01-01")

        case.set_value("CONTINUE_RUN", "FALSE")
        case.set_value("RESUBMIT","0")

        user_nl_fname = os.path.join(case_root, "user_nl_clm")

        user_nl_file  = open(user_nl_fname, "a")
        user_nl_lines = [
                        "hist_mfilt = 20",
                        "hist_nhtfrq = -8760",
                        "hist_empty_htapes = .true.",
                        "hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'"]

        #-- TODO: If we add 
        #-- point to the correct finidat
        fini_line = ("finidat="+post_finidat) 

        for line in user_nl_lines:
            user_nl_file.write("%s\n" % line)

        user_nl_file.write("%s\n" % fini_line)

        user_nl_file.close()

        case.create_namelists()
        case.submit()

def run_transient(orig_root, case_root, user_mods_dir, overwrite, start_year, end_year):
    """
    Function to run transient simultion for a single case.
    For creating the transient case we clone the base case
 
    Args:
        orig_root (str): 
            root of the base_case
            that we want to clone from
        case_root (str):
            Name of the AD case root 
        user_mods_dir (list):
            list of user_mods_dir
            Note, cime does not accept str only
            and it should be a list of strs.
        overwrite (bool) : 
            Flag to overwrite the case if exists
        postad_length (int) :
            Length of the AD spinup in years
        neon_site (NeonSite):
            NeonSite object that we are running.
        start_year (int):
            start year that the user specifies
        end_year (int):
            end year that the user specifies.
        """
    if not os.path.isdir(orig_root):
        sys.exit('Base case does not exist in', orig_root)

    if overwrite and os.path.isdir(case_root):
        shutil.rmtree(case_root)
        logger.info("---- removing the existing case -------")

    if not os.path.isdir(case_root):
        with Case(orig_root, read_only=False) as clone:
            logger.info("---- cloning the base case -------")
            clone.create_clone(case_root, keepexe=True)

    with Case(case_root, read_only=False) as case:
        case.set_value("DATM_YR_ALIGN",start_year)
        case.set_value("DATM_YR_START",start_year)
        case.set_value("DATM_YR_END",end_year)
        case.set_value("STOP_OPTION", "nyears")

        case.set_value("CLM_FORCE_COLDSTART","off")
        case.set_value("CLM_ACCELERATED_SPINUP","off")

        case.set_value("STOP_N", "12")
        case.set_value("REST_N", "12")

        case.set_value("RUN_REFDATE", "2018-01-01")
        case.set_value("RUN_STARTDATE", "2018-01-01")

        case.set_value("CONTINUE_RUN", "FALSE")
        case.set_value("RESUBMIT","0")

        user_nl_fname = os.path.join(case_root, "user_nl_clm")

        # No need to make any changes to this ?
        #Confirm this

        user_nl_file  = open(user_nl_fname, "a")
        #user_nl_lines = [
        #                "hist_mfilt = 48",
        #                "hist_nhtfrq = -1",
        #                "hist_empty_htapes = .true.",
        #                "hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'"]
        # point to the correct finidat
        #for line in user_nl_lines:
        #    user_nl_file.write("%s\n" % line)
        fini_line = ("finidat="+post_finidat) 

        user_nl_file.close()

        case.create_namelists()
        case.submit()


def main(description):
    cesmroot = path_to_ctsm_root()
    # Get the list of supported neon sites from usermods
    valid_neon_sites = glob.glob(os.path.join(cesmroot,"cime_config","usermods_dirs","NEON","[!d]*"))
    valid_neon_sites = [v.split('/')[-1] for v in valid_neon_sites]

    site_list, case_root, run_type, overwrite = get_parser(sys.argv, description, valid_neon_sites)

    logger.debug ("case_root : "+ case_root)

    if not os.path.exists(case_root):
        os.makedirs(case_root)

    #-- check neon listing file for available data: 
    available_list = check_neon_listing(valid_neon_sites)

    #=================================
    #-- all neon sites can be cloned from one generic case
    #-- so no need to define a base_case for every site.

    res = "CLM_USRDAT"
    compset = "I1PtClm51Bgc"

    #--  Looping over neon sites
    orig_root = None
    for neon_site in available_list: 
        if neon_site.name in site_list:
            if not orig_root:

                orig_root = neon_site.build_base_case(cesmroot, case_root, res,
                                                      compset, overwrite)
            logger.info ("-----------------------------------")
            logger.info ("Running CTSM for neon site : {}".format(neon_site.name))


            if (run_type=="ad"):
                ad_length = args.ad_length
                print ("Running Accelerated Decomposition Spinup for: ", args.ad_length)
                print ("-----------------------------------")

                ad_case_name = "spinup_AD_"+neon_site.name
                ad_case_root = os.path.join(out_dir,"neon_sims",ad_case_name) 
                user_mods_dir=["NEON/"+neon_site.name]
                overwrite = True
                run_spinup_ad(orig_root, ad_case_root, user_mods_dir, overwrite, ad_length, neon_site)

            #if args.postad_flag:
            if (run_type=="postad"):
                postad_length = args.postad_length
                print ("Running in post-AD mode for: ", args.ad_length)
                print ("-----------------------------------")
                postad_case_name = "spinup_postAD_"+neon_site.name
                postad_case_root = os.path.join(out_dir,"neon_sims",postad_case_name) 
                user_mods_dir=["NEON/"+neon_site.name]

                run_postad(orig_root, postad_case_root, user_mods_dir, overwrite, postad_length, neon_site)

            #if args.transient_flag:
            if (run_type=="transient"):
                start_year = args.start_year
                end_year = args.start_year
                print ("Running in transient mode for: ", start_year, "to", end_year)
                print ("-----------------------------------")
                transient_case_name = "transient_"+neon_site.name
                transient_case_root = os.path.join(out_dir,"neon_sims",transient_case_name) 
                user_mods_dir=["NEON/"+neon_site.name]

                run_transient(orig_root, transient_case_root, user_mods_dir, overwrite, start_year, end_year)



if __name__ == "__main__":                                                                                                                                  
        main(__doc__) 


