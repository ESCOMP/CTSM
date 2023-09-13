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

    parser.add_argument('--16pft',
                        help='Create and/or modify 16-PFT surface datasets (e.g. for a FATES run) ',
                        action="store_true",
                        dest="pft_16",
                        default=True)

    parser.add_argument('-m', '--mixed',
                        help='Do not overwrite surface dataset to be just one dominant PFT at 100%',
                        action="store_true",
                        dest="mixed",
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


    plumber2_sites = pd.read_csv('PLUMBER2_sites.csv', skiprows=3)  


    for i, row in tqdm.tqdm(plumber2_sites.iterrows()):
        lat = row['Lat']
        lon = row['Lon']
        site = row['Site']
        pft1 = row['pft1']
        pctpft1 = row['pft1-%']
        pft2 = row['pft2']
        pctpft2 = row['pft2-%']
        # overwrite missing values from .csv file
        if pft2 == -999:
           pft2 = 0
           pctpft2 = 0

        clmsite = "1x1_PLUMBER2_"+site
        print ("Now processing site :", site)

        if args.mixed and args.pft_16:
            # use surface dataset with 16 pfts, and don't overwrite with 100% 1 dominant PFT
            # don't set crop flag
            # don't set a dominant pft
            subset_command = ['./subset_data','point','--lat',str(lat),'--lon',str(lon),
                              '--site',clmsite, '--create-surface','--uniform-snowpack',
                              '--cap-saturation','--verbose','--overwrite']
            #modify_command = ['./modify_singlept_site_neon.py', '--neon_site', site, '--surf_dir',
            #           'subset_data_single_point', '--16pft']
        elif args.pft_16:
            # use surface dataset with 16 pfts, but overwrite to 100% 1 dominant PFT
            # don't set crop flag
            # set dominant pft
            subset_command = ['./subset_data','point','--lat',str(lat),'--lon',str(lon),
                              '--site',clmsite,'--dompft',str(pft1),str(pft2), 
                              '--pctpft', str(pctpft1),str(pctpft2), '--create-surface',
                              '--uniform-snowpack','--cap-saturation','--verbose','--overwrite']
        elif args.mixed:
            # use surface dataset with 78 pfts, and don't overwrite with 100% 1 dominant PFT
            # NOTE: FATES will currently not run with a 78-PFT surface dataset
            # set crop flag
            # don't set dominant pft
            subset_command = ['./subset_data','point','--lat',str(lat),'--lon',str(lon),
                              '--site',clmsite,'--crop','--create-surface',
                              '--uniform-snowpack','--cap-saturation','--verbose','--overwrite']
        else:
            # use surface dataset with 78 pfts, and overwrite to 100% 1 dominant PFT
            # NOTE: FATES will currently not run with a 78-PFT surface dataset
            # set crop flag
            # set dominant pft
            subset_command = ['./subset_data', 'point', '--lat', str(lat), '--lon', str(lon),
                              '--site', clmsite,'--crop', '--dompft', str(pft1),str(pft2), 
                              '--pctpft', str(pctpft1),str(pctpft2), '--create-surface',
                              '--uniform-snowpack', '--cap-saturation', '--verbose', '--overwrite']
        execute(subset_command)

if __name__ == "__main__": 
    main()

