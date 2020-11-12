#!/usr/bin/env python

# 2020-11-08                Negin Sobhani

import subprocess
import os
import argparse
import re
from datetime import datetime

valid_opts = {
        'res' :
        ['512x1024','360x720cru','128x256','64x128','48x96','94x192','0.23x0.31','0.47x0.63','0.9x1.25','1.9x2.5','2.5x3.33',
        '4x5,10x15','0.125nldas2','5x5_amazon','1x1_camdenNJ','1x1_vancouverCAN','1x1_mexicocityMEX',
        '1x1_asphaltjungleNJ','1x1_brazil,1x1_urbanc_alpha','1x1_numaIA,1x1_smallvilleIA','0.1x0.1','0.25x0.25','0.5x0.5',
        '3x3min','5x5min','10x10min','0.33x0.33','0.125x0.125','ne4np4,ne16np4','ne30np4.pg2','ne30np4.pg3','ne30np4','ne60np4','ne120np4']
        }

class ctsm_case:
    def __init__ (self, res, glc_nec, ssp_rcp):
        self.res = res
        self.glc_nec = glc_nec
        self_ssp_rcp = ssp_rcp
    def name_nl  (self):
        print ("testing")

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
        parser.add_argument('-y','--year',
                    help='Simulation year to run over.', 
                    action="store",
                    dest="sim_year",
                    required=False,
                    choices = ['1885','1895','1980','1982','2000'],
                    default='2000')
        parser.add_argument('--range',
                    help='Simulation years to run over.', 
                    action="store",
                    dest="sim_range",
                    required=False,
                    choices = ['1850-2000','1850-2005','1850-2100'])
        parser.add_argument('-ge','--glc_nec',
                    help='Number of glacier elevation classes to use' ,
                    action="store",
                    dest="glc_nec",
# it can be any number greater than 0 
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
        parser.add_argument('-l','--dinlc',  #--raw_dir or --rawdata_dir
                    help='/path/of/root/of/input/data',  
                    action="store",
                    dest="input_path",
                    default="/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/")
        parser.add_argument('-d','--debug', 
                    help='Just print out what would happen if ran', 
                    action="store_true", 
                    dest="debug", 
                    default=False)
        return parser

def name_nl  (start_year,end_year, res, ssp_rcp, num_pft):
    """
    Build namelist file name.
    """

    time_stamp = datetime.today().strftime('%y%m%d')
    namelist_fname = "surfdata_"+res+"_"+ssp_rcp+"_"+num_pft+"pfts_CMIP6_"+start_year+'-'+end_year+"_c"+time_stamp+".namelist"
    print ("namelist file created is : ", namelist_fname)

    return namelist_fname 

def build_nl (start_year, end_year, res, ssp_rcp, glc_nec, num_pft, input_path):
    """
    Build the namelist/control file for ****. 
    """

    namelist_fname = name_nl (start_year,end_year, res, ssp_rcp, num_pft)
    namelist_file = open (namelist_fname,'w')

    print (num_pft)
    label = tag_describe()
    print (label)
    nl_template = ( \
            "&clmexp\n"                                                                                                                                                                  
            "nglcec         =  "+glc_nec + "\n"
            "mksrf_fsoitex  =  "+input_path+"mksrf_soitex.10level.c010119.nc"+"\n"
            "mksrf_forganic =  "+input_path+"mksrf_organic_10level_5x5min_ISRIC-WISE-NCSCD_nlev7_c120830.nc"+"\n"
            "mksrf_flakwat  =  "+input_path+"mksrf_LakePnDepth_3x3min_simyr2004_csplk_c151015.nc"+"\n"
            "mksrf_fwetlnd  =  "+input_path+"mksrf_lanwat.050425.nc"+"\n"
            "mksrf_fmax     =  "+input_path+"mksrf_fmax_3x3min_USGS_c120911.nc"+"\n"
            "mksrf_fglacier =  "+input_path+"mksrf_glacier_3x3min_simyr2000.c120926.nc"+"\n"
            "mksrf_fvocef   =  "+input_path+"mksrf_vocef_0.5x0.5_simyr2000.c110531.nc" +"\n"
            "mksrf_furbtopo =  "+input_path+"mksrf_topo.10min.c080912.nc"+"\n"
            "mksrf_fgdp     =  "+input_path+"mksrf_gdp_0.5x0.5_AVHRR_simyr2000.c130228.nc"+"\n"
            "mksrf_fpeat    =  "+input_path+"mksrf_peatf_0.5x0.5_AVHRR_simyr2000.c130228.nc"+"\n"
            "mksrf_fsoildepth ="+input_path+"mksf_soilthk_5x5min_ORNL-Soil_simyr1900-2015_c170630.nc"+"\n"
            "mksrf_fabm     =  "+input_path+"mksrf_abm_0.5x0.5_AVHRR_simyr2000.c130201.nc"+"\n"
            "outnc_double   =  .true. \n"
            "all_urban      = .false.\n"
            "no_inlandwet   = .true.\n"
            "mksrf_furban   =  "+input_path+"mksrf_urban_0.05x0.05_simyr2000.c170724.nc"+"\n"
            "gitdescribe    =  "+label+"\n"
            "mksrf_fdynuse  = ''\n"
            "fdyndat        = ''\n"
            #"numpft         =  "+num_pft+"\n"
            "/\n"
            )

    namelist_file.write (nl_template)
    namelist_file.close ()


def tag_describe ():
    label = subprocess.check_output(["git", "describe"]).strip()
    return label.decode()


def main ():
    print ('Testing gen_mksurf_namelist')

    args         = get_parser().parse_args()

    res          = args.res
    glc_nec      = args.glc_nec
    num_pft      = "78"
    input_path   = args.input_path
    ssp_rcp      = args.ssp_rcp

    sim_year         = args.sim_year
    if args.sim_range:
        run_type = "transient"
        sim_range = args.sim_range
    else:
        run_type = "timeslice"

    if (run_type =="timeslice"):
        start_year = sim_year
        end_year   = sim_year
    elif (run_type =="transient"):
        start_year, end_year = sim_range.split("-")


    build_nl (start_year, end_year, res, ssp_rcp, glc_nec, num_pft, input_path)

if __name__ == "__main__":
    main()

