#! /usr/bin/env python
 
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This is a wrapper script for running CTSM simulation for one or more
neon sites. 

This script is only for neon site and we will develop a more general
code later.

This script will do the following:
    1) Make sure CIME and other dependencies is checked out.
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
#- [ ] Add debug or verbose option
 
#- [ ] Check if it would be better to use cime case obj
#- [ ] Make sure both AD and SASU are not on at the same time
#- [ ]  wget the fields available and run for those available:
#- [ ] Switch to check manage_externals status instead of running it always
#- [ ] query the end year from neon?
 
#- [ ] Matrix spin-up if Eric merged it
 
#Import libraries
from __future__ import print_function
 
import os
import sys
import glob
import shutil
import logging
import requests
import argparse
import subprocess
 
import pandas as pd
 
from datetime import date
from getpass import getuser
 

#cesmroot = '/glade/work/jedwards/sandboxes/ctsm/'
#cesmroot = '/glade/scratch/negins/ctsm_wrapper_oop/'
#cesmroot = '/Users/negins/Desktop/cheyenne/ctsm_wrapper_oop/'

cesmroot = os.environ.get('CESM_ROOT')

cesmroot = '/home/negins/ctsm_cime/'

_LIBDIR = os.path.join(cesmroot,"cime","scripts","Tools")
sys.path.append(_LIBDIR)
_LIBDIR = os.path.join(cesmroot,"cime","scripts","lib")
sys.path.append(_LIBDIR)

import CIME.build as build
from standard_script_setup import *
from CIME.case             import Case
from CIME.utils            import safe_copy
from argparse              import RawTextHelpFormatter
from CIME.locked_files          import lock_file, unlock_file

 
myname = getuser()
 
#-- valid neon site options
valid_neon_sites = ['ABBY','BARR','BART','BLAN',
                    'BONA','CLBJ','CPER','DCFS',
                    'DEJU','DELA','DSNY','GRSM',
                    'GUAN','HARV','HEAL','JERC',
                    'JORN','KONA','KONZ','LAJA',
                    'LENO','MLBS','MOAB','NIWO',
                    'NOGP','OAES','ONAQ','ORNL',
                    'OSBS','PUUM','RMNP','SCBI',
                    'SERC','SJER','SOAP','SRER',
                    'STEI','STER','TALL','TEAK',
                    'TOOL','TREE','UKFS','UNDE',
                    'WOOD','WREF','YELL'
                   ]
 
 
def get_parser():                                                                                                                                           
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
 
    parser.print_usage = parser.print_help
 
    parser.add_argument('--neon_site',
                help='4-letter neon site code.', 
                action="store",
                dest="site_name",
                required=False,
                choices=valid_neon_sites,
                default=["OSBS"],
                #action='append', 
                nargs='+')
    parser.add_argument('--surf_dir',
                help='''
                Directory of single point surface dataset.
                [default: %(default)s]
                ''', 
                action="store", 
                dest="surf_dir",
                type =str,
                required=False,
                default="/glade/scratch/"+myname+"/single_point/")
    parser.add_argument('--out_dir',
                help='''
                Directory to write updated single point surface dataset.
                [default: %(default)s] 
                ''',
                action="store", 
                dest="out_dir",
                type =str,
                required=False,
                default="/glade/scratch/"+myname+"/neon_sims/")
    parser.add_argument('-a','--all',
                help='''
                Flag for running all neon sites, instead of selected few sites.
                [default: %(default)s] 
                ''',
                action="store_true", 
                dest="all_flag",
                default = False)
    parser.add_argument('--start_date',
                help='''
                Start date for running CTSM simulation.
                [default: %(default)s]
                ''', 
                action="store",
                dest="start_date",
                required = False,
                type = int,
                default = 2018)
    parser.add_argument('--end_year',
                help='''
                End year for running CTSM simulation.
                [default: %(default)s]
                ''', 
                action="store",
                dest="end_year",
                required = False,
                type = int,
                default = 2020)
    parser.add_argument('--spinup',
                help='''
                AD spin-up
                [default: %(default)s]
                ''', 
                action="store_true",
                dest="ad_flag",
                required = False,
                default = True)
    parser.add_argument('--sasu','--matrix',
                help='''
                Matrix (SASU) spin-up
                [default: %(default)s]
                ''', 
                action="store_true",
                dest="sasu_flag",
                required = False,
                default = False)
    parser.add_argument('-d','--debug', 
                help='Debug mode will print more information. ', 
                action="store_true", 
                dest="debug", 
                default=False)
                                                                                                                                                            
 
    return parser

class NeonCase (Case):
    """
    A case to encapsulate Neon single point cases.
    ...
    Attributes
    ----------
    plat : float
        latitude
    plon : float
        longitude
    site_name: str -- default = None
        Site name
    Methods
    -------
    create_tag
        create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.
    create_domain_at_point
        Create domain file at a single point.
    create_landuse_at_point:
        Create landuse file at a single point.
    create_surfdata_at_point:
        Create surface dataset at a single point.
    create_datmdomain_at_point:
        Create DATM domain file at a single point.
    """

    def __init__(self, case_root, read_only, compset="I1PtClm51Bgc", res="CLM_USRDAT"):
        super().__init__(case_root, read_only)

        self.compset =compset
        self._compsetname =compset
        self.res = res
        #self._res = res




    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item])
                    for item in self.__dict__))

    def set_compset(self):
        self._compsetname = self.compset

    def build_base_case (self,base_root, base_case_name, res, compset, overwrite,
                        user_mods_dir):

        print (">>>>> BUILDING BASE CASE...<<<<<")
        print ('base_root:',base_root)

        case_root = os.path.join(base_root,base_case_name)
        print ('user_mods_dir:',user_mods_dir)
        print ('case_root:',case_root)
        if overwrite and os.path.isdir(case_root):
            shutil.rmtree(case_root)
            print ("Removing the case...")

        print ("Creating the case....")

        #a_neon_case = NeonCase(case_root, read_only=False)
        #print (a_neon_case)
        #if not os.path.isdir(case_root):
        #print ("not a case")
        #a_neon_case.create(os.path.basename(case_root), cesmroot, compset, res,
        #            machine_name="cheyenne", driver="nuopc",
        #            run_unsupported=True, answer="r",walltime="01:00:00")
        #print (a_neon_case)

        #self.create(os.path.basename(case_root), srcroot=cesmroot, compset_name=self.compset, grid_name=res,
        print (self.__str__())
        if not os.path.isdir(case_root):
            self.create(os.path.basename(case_root), cesmroot, self.compset, res,
                        run_unsupported=True, answer="r",walltime="04:00:00",
                        user_mods_dirs = ["NEON/HARV"], driver="nuopc")

        self.base_case = os.path.basename(case_root)
        self.case_setup()

        build.case_build(case_root, case=self, save_build_provenance=False)

        #exit()
        #with NeonCase(case_root, read_only=False) as neon_case:
        #    print (neon_case.__str__())
        #    if not os.path.isdir(case_root):
        #        neon_case.create(os.path.basename(case_root), cesmroot, compset, res,
        #                    run_unsupported=True, answer="r",walltime="04:00:00",
        #                    user_mods_dirs = ["NEON/HARV"], driver="nuopc")
    #def clone_base_case (self, case_root, case_name, 


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


def check_neon_listing():
    """
    A function to download and parse neon listing file.
    """
    listing_file = 'listing.csv'
    url = 'https://neon-ncar.s3.data.neonscience.org/listing.csv'
 
    download_file(url, listing_file)
    available_list= parse_neon_listing(listing_file)
    return available_list

def parse_neon_listing(listing_file):
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
            print ("Valid neon site",site_name," found!")
            #print(grouped_df.get_group(key), "\n\n")
            tmp_df = grouped_df.get_group(key)
 
            #-- filter files only ending with .nc
            tmp_df = tmp_df[tmp_df[8].str.contains('.nc')]
 
            #-- remove .nc from the file names
            tmp_df[8] = tmp_df[8].str.replace('.nc','')
 
            tmp_df2 = tmp_df[8].str.split("-",  expand=True)
            print (tmp_df2)
            print ("----------------")
        
            start_year = tmp_df2[0].iloc[0]
            end_year = tmp_df2[0].iloc[-1]
 
            start_month = tmp_df2[1].iloc[0]  
            end_month = tmp_df2[1].iloc[-1]  
            print ('start_year=',start_year)
            print ('start_year=',end_year)
            print ('start_month=',start_month)
            print ('end_month=',end_month)
 
            a_neon_site = NeonSite(site_name, start_year, end_year, start_month, end_month)
            print (a_neon_site)
            available_list.append(a_neon_site)
 
    return available_list



#-----------------------------------
def build_base_case2(base_root, basecasename, res, compset, overwrite,
                    user_mods_dir):

    print(">>>>> BUILDING BASE CASE...")
    print ('base_root:',base_root)

    case_root = os.path.join(base_root,basecasename)
    print ('user_mods_dir:',user_mods_dir)
    print ('case_root:',case_root)
    if overwrite and os.path.isdir(case_root):
        shutil.rmtree(case_root)
    print ("Creating the case....")


    print ("Creating the case....")

    print ("user_mods_dir:",user_mods_dir)
    #a_neon_case = NeonCase(case_root, read_only=False)
    #print (a_neon_case)
    #if not os.path.isdir(case_root):
    #print ("not a case")
    #a_neon_case.create(os.path.basename(case_root), cesmroot, compset, res,
    #            machine_name="cheyenne", driver="nuopc",
    #            run_unsupported=True, answer="r",walltime="01:00:00")
    #print (a_neon_case)
    #exit()
    with NeonCase(case_root, read_only=False) as neon_case:
        if not os.path.isdir(case_root):
            neon_case.create(os.path.basename(case_root), cesmroot, compset, res,
                        run_unsupported=True, answer="r",walltime="04:00:00",
                        user_mods_dirs = ["NEON/HARV"], driver="nuopc")

    #with Case(case_root, read_only=False) as case:
    #    if not os.path.isdir(case_root):
    #        #case.create(os.path.basename(case_root), cesmroot, compset, res,
    #        #            driver="nuopc",
    #        #            run_unsupported=True, answer="r",walltime="01:00:00",user_mods_dir=user_mods_dir)
    #        case.create(os.path.basename(case_root), cesmroot, compset, res,
    #                    run_unsupported=True, answer="r",walltime="04:00:00",
    #                    user_mods_dirs = ["NEON/HARV"], driver="nuopc")





    #        print (case.__str__())
    #        print ("case is created!")
    #        # make sure that changing the casename will not affect these variables                                           
    #        case.set_value("EXEROOT",case.get_value("EXEROOT", resolved=True))
    #        case.set_value("RUNDIR",case.get_value("RUNDIR",resolved=True)+".00")
    #     
    #        case.set_value("RUN_TYPE","startup")
    #        case.set_value("GET_REFCASE",False)
         
         
    #    rundir = case.get_value("RUNDIR")
    #    case_root = case.get_value("case_root")
    #    print(">> base case_setup...")
    #    case.case_setup()
    #    print(">> base case_build...")
    #    build.case_build(case_root, case=case)
         
        return case_root

def main():


    args = get_parser().parse_args()

    if args.debug:   
        logging.basicConfig(level=logging.DEBUG)

    #--  specify site/sites to run ctsm
    site_list=args.site_name

    #--  if -a or --all run all neon sites:
    all_flag = args.all_flag
    if all_flag:     
        site_list = valid_neon_sites

    out_dir = args.out_dir
    print ("out_dir:", out_dir)

    #if (not os.path.isdir(out_dir)):
    #    os.mkdir(out_dir)

    neon_site= "ABBY"

    # What to do?
    #1) Create a base case
    #2) Clone base case
    #3) Run base case for your neon

    # But what to do now?
    #) just create a case for your neon

    base_case_name = "base_case"

    base_root = os.path.join("/home/negins/","neon_0725",base_case_name)

    print ("base_root:",base_root)

    res = "CLM_USRDAT"
    compset = "I1PtClm51Bgc"

    #res = "f19_g17"
    #compset = "I1850Clm50Sp"


    #user_mods_dir = os.path.join(cesmroot,"cime_config","usermods_dirs","scam_"+iop)
    user_mods_dir = [os.path.join(cesmroot,"cime_config","usermods_dirs","NEON",neon_site)]

    #ensemble_startval = "001" # The startval strings should be the same length, or else
    #basecase_startval = "000"
    #project = "I1PtClm51Bgc"

    overwrite=True

    # Create and build the base case that all PPE cases are cloned from                                                                                     
    #case_root = build_base_case2(base_root, base_case_name, res,
    #                           compset, overwrite, user_mods_dir)

    print(">>>>> BUILDING BASE CASE...")
    print ('base_root:',base_root)

    case_root = os.path.join(base_root,base_case_name)
    print ('user_mods_dir:',user_mods_dir)
    print ('case_root:',case_root)
    if overwrite and os.path.isdir(case_root):
        shutil.rmtree(case_root)
    print ("Creating the case....")

    print ("user_mods_dir:",user_mods_dir)

    #neon_case = NeonCase(case_root, read_only=False)
    #print (neon_case)

    with NeonCase(case_root, read_only=False) as neon_case:
        if not os.path.isdir(case_root):
            print ("not a case")
            print ("~~~~~~~~~~~~~")
            print (case_root)
            print (os.path.basename(case_root))
            print (cesmroot)
            print ("~~~~~~~~~~~~~")
            print ("~~~~~~~~~~~~~")
            print (neon_case)
            print ("~~~~~~~~~~~~~")

        
            #neon_case.create(os.path.basename(case_root), cesmroot, compset, res,
            #    run_unsupported=True, answer="r",walltime="04:00:00",
            #    user_mods_dirs = ["NEON/HARV"], driver="nuopc")

            #neon_case.build_base_case (base_root, base_case_name, res, compset, overwrite, user_mods_dir)
            neon_case.build_base_case (base_root, base_case_name, res, compset, overwrite, user_mods_dir)
            #neon_case.clone_base_case (base_root, )


    #print (neon_case)
    exit()

    with NeonCase(case_root, read_only=False) as neon_case:
        if not os.path.isdir(case_root):
            neon_case.create(os.path.basename(case_root), cesmroot, compset, res,
                        run_unsupported=True, answer="r",walltime="04:00:00",
                        user_mods_dirs = ["NEON/HARV"], driver="nuopc")
            print (neon_case.__str__())
            print ("case is created!")


    print ("Done!")



if __name__ == "__main__":                                                                                                                                  
        main() 


