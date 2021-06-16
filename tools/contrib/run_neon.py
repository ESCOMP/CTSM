#! /usr/bin/env python


"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This is a wrapper script for running CTSM simulation for one for more
neon sites. 

This script will do the following:
    1) 

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

#Import libraries
from __future__ import print_function

import os
import sys
#import glob
import shutil
import argparse
#import requests
import subprocess

#import numpy as np
#import pandas as pd
#import xarray as xr

from datetime import date
from getpass import getuser


myname = getuser()

neon_sites = ['ABBY','BARR','BART','BLAN',
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
              'WOOD','WREF','YELL'] 


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
                choices=neon_sites,
                default=["OSBS"],
                #action='append', 
                nargs='+')
    parser.add_argument('--surf_dir',
                help='Directory of single point surface dataset. [default: %(default)s]', 
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
                default="/glade/scratch/"+myname+"/single_point_neon_updated/")
    parser.add_argument('-a','--all',
                help='''
                Flag for running all neon sites, instead of selected few sites.
                [default: %(default)s] 
                ''',
                action="store_true", 
                dest="all_flag",
                default = False)
    return parser 

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


def execute_2(command):
    print ('\n',' >>  ',*command,'\n')

    try:
        subprocess.check_call(command, stdout=sys.stdout, stderr=subprocess.STDOUT)

    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))



def main():

    args = get_parser().parse_args()

    #--  specify site/sites from which to extract data

    site_list=args.site_name

    #--  if -a or --all Run all neon sites:
    all_flag = args.all_flag
    if all_flag:
        site_list = neon_sites


    out_dir = args.out_dir

    if (not os.path.isdir(out_dir)):
        os.mkdir(out_dir)


    #--  Looping over neon_sites

    for site_name in site_list:

        case_dir = os.path.join(out_dir, 'NEON_'+site_name)

        if (os.path.isdir(case_dir)):
            shutil.rmtree(case_dir)

        user_mod = 'NEON/'+site_name

        print ("-----------------------------------")
        print ("Running CTSM for neon site : ", site_name)
        print ("-----------------------------------")

        #1 - clone CTSM repository/tag?

        #2 - run manage externals to make sure CIME exist:

        # First go to top level directory:
        os.chdir('../..')
        print("Current working directory: {0}".format(os.getcwd()))

        execute_2(['./manage_externals/checkout_externals','-vv'])

        #3 - Create case for neon site

        os.chdir('cime/scripts/')

        print("Current working directory: {0}".format(os.getcwd()))

        execute_2(['./create_newcase','--case', case_dir, 
                   '--compset','I1PtClm51Bgc',
                   '--res','CLM_USRDAT',
                   '--driver','nuopc',
                   '--user-mods',user_mod,
                   '--handle-preexisting-dirs','u'])

        #4 - Case setup
        os.chdir (case_dir)
        print("Current working directory: {0}".format(os.getcwd()))

        execute_2(['./case.setup'])

        #5 - Make the changes 
        #execute_2(['./xmlchange','DATM_YR_END=2020'])
        #execute_2(['./xmlchange','STOP_OPTION=nyears'])

        execute_2(['./case.build'])
        execute_2(['./case.submit'])


if __name__ == "__main__": 
    main()


