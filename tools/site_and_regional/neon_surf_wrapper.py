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

Instructions for running on Cheyenne/Casper:
load the following into your local environment
    module load python
    ncar_pylib

"""
# TODO
# Automatic downloading of missing files if they are missing
#-[ ] Download neon sites and dom pft file
#-[ ] Make sure verbose works for printing out commands running

#  Import libraries
from __future__ import print_function

import os
import sys
import tqdm
import logging
import argparse
import subprocess

import pandas as pd
#import tqdm as tqdm



def get_parser():                                                                                                                                                                   
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                           formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.print_usage = parser.print_help

    parser.add_argument('-v','--verbose', 
                help='Verbose mode will print more information. ', 
                action="store_true", 
                dest="verbose", 
                default=False)


    return parser


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
        subprocess.check_call(command, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)

    except subprocess.CalledProcessError as e:
        #raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        #print (e.ouput)
        print (e)






def main():

    args = get_parser().parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)


    neon_sites = pd.read_csv('neon_sites_dompft.csv')


    for i, row in tqdm.tqdm(neon_sites.iterrows()):
        lat = row['Lat']
        lon = row['Lon']
        site = row['Site']
        pft = row['pft']
        print ("Now processing site :", site)
        command = ['./subset_data.py','point','--lat',str(lat),'--lon',str(lon),'--site',site,'--dompft',str(pft),'--crop']
        execute(command)

        command = ['./modify_singlept_site_neon.py','--neon_site',site]
        execute(command)

if __name__ == "__main__": 
    main()

