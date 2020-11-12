#!/usr/bin/env python

# 2020-11-08                Negin Sobhani

import subprocess
import os
import argparse
import glob
import re
import datetime 

def mk_dir(dir):
        cmd_line  = "mkdir "+ dir 
        #print cmd_line
        subprocess.call(["mkdir", dir]) 

def get_parser():
        """Get parser object for this script."""
        from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
        parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
        parser = argparse.ArgumentParser(description='Setting the flags and input files.')

        parser.add_argument('-r','--res'
                    , help='Resolution is the supported resolution(s) to use for files.'
                    , action="store"
                    , dest="res"
                    , required=False
                    , default="4x5")
        #parser.add_argument('-y','--years', help='Simulation year(s) to run over.' ,action="store", dest="years",
        #        default="2000", type = valid_date)
        parser.add_argument('-sy','--start_year'
                    , help='Start year for simulation year(s) to run over.'
                    , action="store"
                    , dest="start_year" 
                    , default="2000")
        
        parser.add_argument('-ey','--end_year'  , 
                    help='End year for simulation year(s) to run over.'   ,
                    action="store",
                    dest="years",
                    default="2000")

        parser.add_argument('-ge','--glc_nec',
                    help='Number of glacier elevation classes to use' ,
                    action="store",
                    dest="glc_nec",
                    default = "10")
        
        parser.add_argument('--rundir', 
                    help='Directory to run in.' ,
                    action="store",
                    dest="run_dir", 
                    required = False, 
                    default =os.getcwd())
        parser.add_argument('--ssp_rcp',
                    help='Shared Socioeconomic Pathway and Representative Concentration Pathway Scenario name(s).' ,
                    action="store",
                    dest="ssp_rcp", 
                    required = False,
                    choices = ["hist","SSP1-2.6","SSP3-7.0","SSP5-3.4","SSP2-4.5","SSP1-1.9","SSP4-3.4","SSP4-6.0","SSP5-8.5"],
                    default = "hist")
        parser.add_argument('-l','--dinlc', 
                    help='/path/of/root/of/input/data',  
                    action="store",
                    dest="input_path",
                    default="/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/")
        parser.add_argument('-d','--debug', help='Just print out what would happen if ran'
                ,action="store_true", dest="debug", default=False)
        return parser

def valid_date(s):
        try:
            return datetime.strptime(s, "%Y")
        except ValueError:
            msg = "Not a valid date: '{0}'.".format(s)
            raise argparse.ArgumentTypeError(msg)

def name_nl  (year, res, glc_nec, num_pft):
    """
    Build namelist file name.
    """

    namelist_fname = "surfdata_"+res+"hist_"+num_pft+"pfts_CMIP6_simyr2000_c201028"+".namelist"
    print ("namelist file created is : ", namelist_fname)

    return namelist_fname 

def build_nl (year, res, glc_nec, num_pft, input_path):
    """
    Build the namelist/control file for ****. 
    """

    namelist_fname = name_nl (year, res, glc_nec, num_pft)
    namelist_file = open (namelist_fname,'w')

    nl_template = ( \
            "&clmexp\n"                                                                                                                                                                  
            "nglcec         =  "+glc_nec + "\n"
            "mksrf_fsoitex  =  "+input_path+"mksrf_soitex.10level.c010119.nc"+"\n"
            "mksrf_forganic =  "+input_path+"mksrf_organic_10level_5x5min_ISRIC-WISE-NCSCD_nlev7_c120830.nc"+"\n"
            "outnc_double   =  .true. \n"
            "all_urban      = .false.\n"
            "no_inlandwet   = .true.\n"
            "mksrf_fdynuse  = ''\n"
            "fdyndat        = ''\n"
            "numpft         =  "+num_pft+"\n"
            "/\n"
            )

    namelist_file.write (nl_template)
    namelist_file.close ()




def main ():
    print ('Testing gen_mksurf_namelist')

    args         = get_parser().parse_args()
    year         = args.years
    res          = args.res
    glc_nec      = args.glc_nec
    num_pft      = "78"
    input_path   = args.input_path

    ssp_rcp      = args.ssp_rcp
    
    build_nl (year, res, glc_nec, num_pft, input_path)

    #label = subprocess.check_output(["git", "describe"]).strip()


if __name__ == "__main__":
    main()

