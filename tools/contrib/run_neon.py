#! /usr/bin/env python

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This is a wrapper script for running CTSM simulation for one or more
neon sites. 

This script will do the following:
    1) Make sure CIME and other dependencies is checked out.
    2) Make the case for the specific neon site(s).
    3) Make changes to the case, for:
        a. AD spinup
        b. SASU spinup
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

#QUESTIONS (NS):
#- [ ] Matrix spin-up if Eric merged it


#Import libraries
from __future__ import print_function

import os
import sys
import shutil
import logging
import argparse
import requests
import subprocess

import pandas as pd

from datetime import date
from getpass import getuser


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

"""
def execute(cmd):

    print ('\n',' >>  ',*cmd,'\n')
    process = subprocess.Popen(cmd, shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)

    for stdout_line in iter(process.stdout.readline, ""):
        yield stdout_line 
    process.stdout.close()
    return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
"""

def execute(command):
    """
    Function for running a command on shell.

    Args:
        command (str):
            command that we want to run.

    Raises:
        Error with the return code from shell.
    """
    print ('\n',' >>  ',*command,'\n')

    try:
        subprocess.check_call(command, stdout=sys.stdout, stderr=subprocess.STDOUT)

    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))


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
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item])
                    for item in (self.__dict__)))


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

    if (not os.path.isdir(out_dir)):
        os.mkdir(out_dir)


    #-- check neon listing file for available data:
    available_list = check_neon_listing()


    #1 - run manage externals to make sure CIME exist:

    #-- First go to top level directory:
    clone_dir = os.path.abspath(os.path.join(__file__, "../../.."))
    print (clone_dir)

    os.chdir(clone_dir)
    logging.debug("Current working directory: {0}".format(os.getcwd()))

    execute(['./manage_externals/checkout_externals','-vv'])

    #--  Looping over neon sites

    for neon_site in available_list:
        if neon_site.name in site_list:

            print ("-----------------------------------")
            print ("Running CTSM for neon site : ", neon_site.name)
            print ("-----------------------------------")

            case_dir = os.path.join(out_dir, 'NEON_'+neon_site.name)

            #-- remove the case if it exists
            if (os.path.isdir(case_dir)):
                shutil.rmtree(case_dir)

            user_mod = 'NEON/'+neon_site.name


            #3 - Create case for neon site

            os.chdir('cime/scripts/')
            print("Current working directory: {0}".format(os.getcwd()))

            execute(['./create_newcase','--case', case_dir, 
                       '--compset','I1PtClm51Bgc',
                       '--res','CLM_USRDAT',
                       '--driver','nuopc',
                       '--user-mods',user_mod,
                       '--handle-preexisting-dirs','u'])

            #4 - Case setup
            os.chdir (case_dir)
            print("Current working directory: {0}".format(os.getcwd()))

            execute(['./case.setup'])

            #5 - Make the changes
            execute(['./xmlchange','DATM_YR_END=2020'])
            execute(['./xmlchange','STOP_OPTION=nyears'])


            # Spin-up AD mode:
            if args.ad_flag:

            #-- cycle over available input data
                #execute(['./xmlchange','DATM_CLMNCEP_YR_START='+neon_site.start_year])
                #execute(['./xmlchange','DATM_CLMNCEP_YR_END='+neon_site.end_year])

                execute(['./xmlchange','DATM_YR_START='+neon_site.start_year])
                execute(['./xmlchange','DATM_YR_END=2020'])
                #execute(['./xmlchange','DATM_YR_END='+neon_site.end_year])
                execute(['./xmlchange','STOP_OPTION=nyears'])

                execute(['./xmlchange','CLM_FORCE_COLDSTART=on'])
                execute(['./xmlchange','CLM_ACCELERATED_SPINUP=on'])
                execute(['./xmlchange','STOP_N=400'])
                execute(['./xmlchange','REST_N=100'])

                execute(['./xmlchange','RUN_REFDATE=0018-01-01'])
                execute(['./xmlchange','RUN_STARTDATE=0018-01-01'])

                execute(['./xmlchange','CONTINUE_RUN=FALSE'])
                execute(['./xmlchange','RESUBMIT=0'])


                #execute(['./xmlchange','CLM_FORCE_COLDSTART=on'])
                #execute(['./xmlchange','CLM_ACCELERATED_SPINUP=on'])
                #execute(['./xmlchange','STOP_N=400'])
                #execute(['./xmlchange','REST_N=100'])
                #execute(['./xmlchange','RUN_REFDATE=0018-01-01'])
                #execute(['./xmlchange','RUN_STARTDATE=0018-01-01'])
                #execute(['./xmlchange','RUN_STARTDATE=0018-01-01'])
                #execute(['./xmlchange','CONTINUE_RUN=FALSE'])
                #execute(['./xmlchange','RESUBMIT=0'])

                print("Current working directory: {0}".format(os.getcwd()))
                execute(['./preview_namelists'])

            #6- TODO: Changes for Matrix ...


            #on cheyenne:
            # execute(['qcmd --','./case.build'])
            execute(['./case.build'])
            execute(['./case.submit'])
            os.chdir(clone_dir)

            #postAD spin-up

            #transient
            #exit()


if __name__ == "__main__": 
    main()


