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
import re
import subprocess
import pandas as pd
import glob 
import datetime
from getpass import getuser
 
# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..","..",'python'))
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root

import CIME.build as build
from standard_script_setup import *
from CIME.case             import Case
from CIME.utils            import safe_copy, expect, symlink_force
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

    parser.add_argument('--base-case',
                help='''
                Root Directory of base case build
                [default: %(default)s] 
                ''',
                action="store", 
                dest="base_case_root",
                type =str,
                required=False,
                default=None)

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

    parser.add_argument('--overwrite',
                help='''
                overwrite existing case directories
                [default: %(default)s]
                ''', 
                action="store_true",
                dest="overwrite",
                required = False,
                default = False)

    parser.add_argument('--run-type',
                        help='''
                        Type of run to do
                        [default: %(default)s]
                        ''',
                        choices = ["ad", "postad", "transient", "sasu"],
                        default = "transient")

    parser.add_argument ('--run-length',
                help='''
                How long to run (modified ISO 8601 duration)
                [default: %(default)s]
                ''', 
                required = False,
                type = str, 
                default = '200Y')

    parser.add_argument('--start-date',
                help='''           
                Start date for running CTSM simulation in ISO format.
                [default: %(default)s]
                ''',               
                action="store", 
                dest="start_date",                                                                                     
                required = False,
                type = datetime.date.fromisoformat,
                default = datetime.datetime.strptime("2018-01-01",'%Y-%m-%d'))

    parser.add_argument('--end-date',
                help='''
                End date for running CTSM simulation in ISO format.
                [default: %(default)s]
                ''', 
                action="store",
                dest="end_date",
                required = False,
                type = datetime.date.fromisoformat,
                default = datetime.datetime.strptime("2020-01-01",'%Y-%m-%d'))

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


    run_length = parse_isoduration(args.run_length)
    
    return neon_sites, args.case_root, args.run_type, args.overwrite, run_length, args.base_case_root

def get_isosplit(s, split):
    if split in s:
        n, s = s.split(split)
    else:
        n = 0
    return n, s

def parse_isoduration(s):
    '''
    simple ISO 8601 duration parser, does not account for leap years and assumes 30 day months
    '''
    # Remove prefix
    s = s.split('P')[-1]
    
    # Step through letter dividers
    years, s = get_isosplit(s, 'Y')
    months, s = get_isosplit(s, 'M')
    days, s = get_isosplit(s, 'D')
    
    # Convert all to timedelta
    dt = datetime.timedelta(days=int(days)+365*int(years)+30*int(months))
    return int(dt.total_seconds()/86400)

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
        self.start_year= int(start_year)
        self.end_year = int(end_year)
        self.start_month = int(start_month)
        self.end_month = int(end_month)
        self.cesmroot = path_to_ctsm_root()
        
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
        self.base_case_root = case_root
        user_mods_dirs = [os.path.join(cesmroot,"cime_config","usermods_dirs","NEON",self.name)]
        if not case_root:
            case_root = os.getcwd()
        case_path = os.path.join(case_root,self.name)
            
        logger.info ('base_case_name : {}'.format(self.name))
        logger.info ('user_mods_dir  : {}'.format(user_mods_dirs[0]))

        if overwrite and os.path.isdir(case_path):
            logger.info ("Removing the existing case at: {}".format(case_path))
            shutil.rmtree(case_path)
            
        with Case(case_path, read_only=False) as case:
            if not os.path.isdir(case_path):
                logger.info("---- creating a base case -------")

                case.create(case_path, cesmroot, compset, res, mpilib="mpi-serial",
                            run_unsupported=True, answer="r",output_root=case_root,
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

    def diff_month(self):
        d1 = datetime(self.end_year,self.end_month, 1)
        d2 = datetime(self.start_year, self.start_month, 1)
        return (d1.year - d2.year) * 12 + d1.month - d2.month
    

    
    def run_case(self, base_case_root, run_type, run_length, overwrite=False):
        user_mods_dirs = [os.path.join(self.cesmroot,"cime_config","usermods_dirs","NEON",self.name)]
        expect(os.path.isdir(base_case_root), "Error base case does not exist in {}".format(base_case_root))
        case_root = os.path.abspath(os.path.join(base_case_root,"..", self.name+"."+run_type))
        if os.path.isdir(case_root) and overwrite:
            logger.info("---- removing the existing case -------")
            shutil.rmtree(case_root)
        if not os.path.isdir(case_root):
            # read_only = False should not be required here
            with Case(base_case_root, read_only=False) as basecase:
                logger.info("---- cloning the base case in {}".format(case_root))
                basecase.create_clone(case_root, keepexe=True, user_mods_dirs=user_mods_dirs)


        with Case(case_root, read_only=False) as case:
            case.set_value("STOP_OPTION", "ndays")
            case.set_value("STOP_N", run_length)
            case.set_value("REST_N", 100)
            case.set_value("CONTINUE_RUN", False)
            
            if run_type == "ad":
                case.set_value("CLM_FORCE_COLDSTART","on")
                case.set_value("CLM_ACCELERATED_SPINUP","on")
                case.set_value("RUN_REFDATE", "0018-01-01")
                case.set_value("RUN_STARTDATE", "0018-01-01")                
                case.set_value("JOB_WALLCLOCK_TIME", "12:00:00", subgroup="case.run")
            else:
                case.set_value("CLM_FORCE_COLDSTART","off")
                case.set_value("CLM_ACCELERATED_SPINUP","off")
                case.set_value("RUN_TYPE", "hybrid")
                
            if run_type == "postad":
                self.set_ref_case(case)
                
            if run_type == "transient":
                self.set_ref_case(case)
                case.set_value("STOP_OPTION","nmonths")
                case.set_value("STOP_N", self.diff_month())
                case.set_value("REST_N", "12")
                case.set_value("DATM_YR_ALIGN",self.start_year)
                case.set_value("DATM_YR_START",self.start_year)
                case.set_value("DATM_YR_END",self.end_year)
            else:
                # for the spinup we want the start and end on year boundaries
                if self.start_month == 1:
                    case.set_value("DATM_YR_ALIGN",self.start_year)
                    case.set_value("DATM_YR_START",self.start_year)
                elif self.start_year + 1 <= self.end_year:
                    case.set_value("DATM_YR_ALIGN",self.start_year+1)
                    case.set_value("DATM_YR_START",self.start_year+1)
                if self.end_month == 12:
                    case.set_value("DATM_YR_END",self.end_year)
                else:
                    case.set_value("DATM_YR_END",self.end_year-1)

                    
                self.modify_user_nl(case_root, run_type)
                
            case.create_namelists()
            case.submit()

    def set_ref_case(self, case):
        rundir = case.get_value("RUNDIR")
        case_root = case.get_value("CASEROOT")
        if case_root.endswith(".postad"):
            ref_case_root = case_root.replace(".postad",".ad")
            root = ".ad"
        else:
            ref_case_root = case_root.replace(".transient",".postad")
            root = ".postad"
        expect(os.path.isdir(ref_case_root), "ERROR: spinup must be completed first, could not find directory {}".format(ref_case_root))
        with Case(ref_case_root) as refcase:
            refrundir = refcase.get_value("RUNDIR")
        case.set_value("RUN_REFDIR", refrundir)
        case.set_value("RUN_REFCASE", os.path.basename(ref_case_root))
        for reffile in glob.iglob(refrundir + "/{}{}.*.nc".format(self.name, root)):
            m = re.search("(\d\d\d\d-\d\d-\d\d)-\d\d\d\d\d.nc", reffile)
            if m:
                refdate = m.group(1)
            symlink_force(reffile, os.path.join(rundir,os.path.basename(reffile)))
        for rpfile in glob.iglob(refrundir + "/rpointer*"):
            safe_copy(rpfile, rundir)
        if not os.path.isdir(os.path.join(rundir, "inputdata")) and os.path.isdir(os.path.join(refrundir,"inputdata")):
            symlink_force(os.path.join(refrundir,"inputdata"),os.path.join(rundir,"inputdata"))
        case.set_value("RUN_REFDATE", refdate)
        if case_root.endswith(".postad"):
            case.set_value("RUN_STARTDATE", refdate)
        else:
            case.set_value("RUN_STARTDATE", "{yr:04d}-{mo:02d}-01".format(yr=self.start_year, mo=self.start_month))
                
        
    def modify_user_nl(self, case_root, run_type):
        user_nl_fname = os.path.join(case_root, "user_nl_clm")

        if run_type != "transient":
            user_nl_lines = [
                "hist_fincl2 = ''",
                "hist_mfilt = 20",
                "hist_nhtfrq = -8760",
                "hist_empty_htapes = .true.",
                "hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'"]
                        
        with open(user_nl_fname, "a") as fd:
            for line in user_nl_lines:
                fd.write("{}\n".format(line))

        

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
 
    #-- filter lines with atm/cdep
    df = df[df['object'].str.contains("atm/cdeps/")]

    #-- split the object str to extract site name
    df=df['object'].str.split("/", expand=True)
 
    #-- groupby site name
    grouped_df = df.groupby(7)
    for key, item in grouped_df:
        #-- check if it is a valid neon site
        if any(key in x for x in valid_neon_sites):
            site_name = key

            tmp_df = grouped_df.get_group(key)
 
            #-- filter files only ending with YYYY-MM.nc
            tmp_df = tmp_df[tmp_df[8].str.contains('\d\d\d\d-\d\d.nc')]
            latest_version = tmp_df[6].iloc[-1]
            tmp_df = tmp_df[tmp_df[6].str.contains(latest_version)]
            #-- remove .nc from the file names
            tmp_df[8] = tmp_df[8].str.replace('.nc','')
 
            tmp_df2 = tmp_df[8].str.split("-",  expand=True)
            # ignore any prefix in file name and just get year
            tmp_df2[0] = tmp_df2[0].str.slice(-4)
            #-- figure out start_year and end_year
            start_year = int(tmp_df2[0].iloc[0])
            end_year = int(tmp_df2[0].iloc[-1])
 
            #-- figure out start_month and end_month
            start_month = int(tmp_df2[1].iloc[0])
            end_month = int(tmp_df2[1].iloc[-1])

            logger.debug ("Valid neon site " + site_name+" found!")
            logger.debug ("File version {}".format(latest_version))
            logger.debug ('start_year={}'.format(start_year))
            logger.debug ('end_year={}'.format(end_year))
            logger.debug ('start_month={}'.format(start_month))
            logger.debug ('end_month={}'.format(end_month))
 
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

def main(description):
    cesmroot = path_to_ctsm_root()
    # Get the list of supported neon sites from usermods
    valid_neon_sites = glob.glob(os.path.join(cesmroot,"cime_config","usermods_dirs","NEON","[!d]*"))
    valid_neon_sites = [v.split('/')[-1] for v in valid_neon_sites]

    site_list, case_root, run_type, overwrite, run_length, base_case_root = get_parser(sys.argv, description, valid_neon_sites)
    if case_root:
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

    for neon_site in available_list: 
        if neon_site.name in site_list:
            if not base_case_root:
                base_case_root = neon_site.build_base_case(cesmroot, case_root, res,
                                                      compset, overwrite)
            logger.info ("-----------------------------------")
            logger.info ("Running CTSM for neon site : {}".format(neon_site.name))
            neon_site.run_case(base_case_root, run_type, run_length, overwrite)

if __name__ == "__main__":                                                                                                                                  
        main(__doc__) 


