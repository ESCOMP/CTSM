#!/usr/bin/env python

# 2020-11-08                Negin Sobhani

import subprocess
import os
import argparse
import re
from datetime import datetime
import logging
import sys

valid_opts = {
        'res' :
        ['512x1024','360x720cru','128x256','64x128','48x96','94x192','0.23x0.31','0.47x0.63','0.9x1.25','1.9x2.5','2.5x3.33',
        '4x5,10x15','0.125nldas2','5x5_amazon','1x1_camdenNJ','1x1_vancouverCAN','1x1_mexicocityMEX',
        '1x1_asphaltjungleNJ','1x1_brazil,1x1_urbanc_alpha','1x1_numaIA,1x1_smallvilleIA','0.1x0.1','0.25x0.25','0.5x0.5',
        '3x3min','5x5min','10x10min','0.33x0.33','0.125x0.125','ne4np4,ne16np4','ne30np4.pg2','ne30np4.pg3','ne30np4','ne60np4','ne120np4']
        ,'ssp_rcp': ["hist","SSP1-2.6","SSP3-7.0","SSP5-3.4","SSP2-4.5","SSP1-1.9","SSP4-3.4","SSP4-6.0","SSP5-8.5"]
        }

class ctsm_case:
    def __init__ (self, res, glc_nec, ssp_rcp):
        self.res = res
        self.glc_nec = glc_nec
        self_ssp_rcp = ssp_rcp
    def name_nl  (self):
        print ("testing")

def get_parser():
## Add default values in the help page.
        """Get parser object for this script."""
        from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
        parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
        parser = argparse.ArgumentParser(description='Setting the flags and input files.')

        parser.add_argument('-r','--res',
                    help='Resolution is the supported resolution(s) to use for files.',
                    action="store", 
                    dest="res", 
                    choices=valid_opts['res'], 
                    required=False, 
                    default="4x5")
        parser.add_argument('-sy','--start_year',
                    help='Simulation start year.', 
                    action="store",
                    dest="start_year",
                    required=False,
                    type = start_year_type,
                    default=2000)
        parser.add_argument('-ey','--end_year',
## Add information about if this is optional
## Add info for help page note -- by default is start_year
                    help='Simulation end year.', 
                    action="store",
                    dest="end_year",
                    required=False,
                    type = int)
                    #type = end_year_type,
                    #default="start_year")
# Possibly remove year --years and range options
# comment them out
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
###############################################
        parser.add_argument('-ge','--glc_nec',
                    help='Number of glacier elevation classes to use' ,
                    action="store",
                    dest="glc_nec",
                    type = glc_nec_type,
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
                    choices=valid_opts['ssp_rcp'], 
                    default = "hist")
        parser.add_argument('-l','--dinlc',  #--raw_dir or --rawdata_dir
                    help='/path/of/root/of/input/data',  
                    action="store",
                    dest="input_path",
                    default="/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/")
#maybe a verbose option and removing debug
        parser.add_argument('-d','--debug', 
                    help='Just print out what would happen if ran', 
                    action="store_true", 
                    dest="debug", 
                    default=False)
        parser.add_argument('-vic','--vic', 
                    help='Add the fields required for the VIC model', 
                    action="store_true", 
                    dest="vic_flag", 
                    default=False)
        parser.add_argument('-glc','--glc', 
                    help='Add the optional 3D glacier fields for verification of the glacier model', 
                    action="store_true", 
                    dest="glc_flag", 
                    default=False)
        parser.add_argument('-hirespft','--hirespft', 
                    help='If you want to use the high-resolution pft dataset rather \n' 
                        'than the default lower resolution dataset \n'
#add error check for hi-res and years if they are 1850 and 2005.
                        '(low resolution is at quarter-degree, high resolution at 3minute)\n'
                        '(hires only available for 1850 and 2005))', 
                    action="store_true", 
                    dest="hres_flag", 
                    default=False)
        parser.add_argument('-nocrop','--nocrop', 
                    help='Create datasets without the extensive list of prognostic crop types', 
                    action="store_false", 
                    dest="crop_flag", 
                    default=True)
        parser.add_argument('-f','--fast', 
                    help='Toggle fast mode which does not user the large mapping file', 
                    action="store_true", 
                    dest="fast_flag", 
                    default=False)
        return parser

def name_nl  (start_year,end_year, res, ssp_rcp, num_pft):
    """
    Build namelist file name.
    """

    time_stamp = datetime.today().strftime('%y%m%d')
    namelist_fname = "surfdata_"+res+"_"+ssp_rcp+"_"+num_pft+"pfts_CMIP6_"+start_year.__str__()+'-'+end_year.__str__()+"_c"+time_stamp+".namelist"
    print ("namelist file created is : ", namelist_fname)

    return namelist_fname 

def build_nl (start_year, end_year, res, ssp_rcp, glc_nec, num_pft, input_path, run_type, vic_flag, glc_flag ):
    """
    Build the namelist/control file for ****. 
    """

    namelist_fname = name_nl (start_year,end_year, res, ssp_rcp, num_pft)
    namelist_file = open (namelist_fname,'w')

    label = tag_describe()
    
    if (run_type == "transient"):
        use_transient = ".true."
    else: 
        use_transient = ".false"
    
    nl_template = ( \
            "&clmexp\n"                                                                                                                                                                  
            "nglcec           = "+glc_nec + "\n"
            "mksrf_fsoitex    = "+input_path+"mksrf_soitex.10level.c010119.nc"+"\n"
            "mksrf_forganic   = "+input_path+"mksrf_organic_10level_5x5min_ISRIC-WISE-NCSCD_nlev7_c120830.nc"+"\n"
            "mksrf_flakwat    = "+input_path+"mksrf_LakePnDepth_3x3min_simyr2004_csplk_c151015.nc"+"\n"
            "mksrf_fwetlnd    = "+input_path+"mksrf_lanwat.050425.nc"+"\n"
            "mksrf_fmax       = "+input_path+"mksrf_fmax_3x3min_USGS_c120911.nc"+"\n"
            "mksrf_fglacier   = "+input_path+"mksrf_glacier_3x3min_simyr2000.c120926.nc"+"\n"
            "mksrf_fvocef     = "+input_path+"mksrf_vocef_0.5x0.5_simyr2000.c110531.nc" +"\n"
            "mksrf_furbtopo   = "+input_path+"mksrf_topo.10min.c080912.nc"+"\n"
            "mksrf_fgdp       = "+input_path+"mksrf_gdp_0.5x0.5_AVHRR_simyr2000.c130228.nc"+"\n"
            "mksrf_fpeat      = "+input_path+"mksrf_peatf_0.5x0.5_AVHRR_simyr2000.c130228.nc"+"\n"
            "mksrf_fsoildepth = "+input_path+"mksf_soilthk_5x5min_ORNL-Soil_simyr1900-2015_c170630.nc"+"\n"
            "mksrf_fabm       = "+input_path+"mksrf_abm_0.5x0.5_AVHRR_simyr2000.c130201.nc"+"\n"
            "outnc_double     = .true. \n"
            "all_urban        = .false.\n"
            "no_inlandwet     = .true. \n"
            "mksrf_furban     = "+input_path+"mksrf_urban_0.05x0.05_simyr2000.c170724.nc"+"\n"
            "gitdescribe      = "+label+"\n"
            "mksrf_ftopostats = "+input_path+"mksrf_topostats_1km-merge-10min_HYDRO1K-merge-nomask_simyr2000.c130402.nc"+"\n"
            "mksrf_fvegtyp    = "+input_path + "pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629/mksrf_landuse_histclm50_LUH2_1850.c170629.nc"+"\n"
            "mksrf_fsoicol    = "+input_path + "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_soilcolor_CMIP6_simyr2005.c170623.nc"+"\n"
            "mksrf_flai       = "+input_path + "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_lai_78pfts_simyr2005.c170413.nc" +"\n"
            "mksrf_fdynuse    = ''\n"
            "fdyndat          = ''\n"
            "numpft           = "+num_pft+"\n"
            "dst_mesh_file    = \n"
            "\n&transient\n"
            "use_transient    = "+use_transient + "\n"
            "mksrf_fdynuse    = \n"
            "start_year       = "+start_year.__str__() + "\n"
            "end_year         = "+end_year.__str__() + "\n"
            "mksrf_dyn_lu     = "+input_path+"\n"
            "\n&vic\n"
            "use_vic          = "+vic_flag.__str__()+"\n"
            "mksrf_fvic       = "+input_path+"mksrf_vic_0.9x1.25_GRDC_simyr2000.c130307.nc\n"
            "outnc_vic        = \n"
            "\n&glc\n"
            "use_glc          = "+glc_flag.__str__()+"\n"
            "outnc_3dglc      = \n"
            "/\n"
            )

    namelist_file.write (nl_template)
    namelist_file.close ()


def tag_describe ():
    label = subprocess.check_output(["git", "describe"]).strip()
    return label.decode()

def glc_nec_type(x):
    x = int(x)
    if (x <= 0) or (x >= 100):
        raise argparse.ArgumentTypeError("ERROR: glc_nec must be between 1 and 99.")
    return x

def start_year_type(x):
    x = int(x)
    if (x <= 850) or (x >= 2105):
        raise argparse.ArgumentTypeError("ERROR: Simulation start year should be between 850 and 2105.")
    return x

def check_endyear(start_year, end_year):
    if (end_year < start_year):
        print ( "ERROR: end_year should be bigger than the start_year : ", start_year, ".")
        sys.exit()

def main ():
    print ('Testing gen_mksurf_namelist')

    args         = get_parser().parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)


    res          = args.res
    glc_nec      = args.glc_nec.__str__()
    input_path   = args.input_path
    ssp_rcp      = args.ssp_rcp
    crop_flag    = args.crop_flag
    vic_flag     = args.vic_flag
    glc_flag     = args.glc_flag

    start_year   = args.start_year
    end_year     = args.end_year

    # determine end_year
    check_endyear (start_year, end_year)
    if not end_year:
        end_year = start_year

    if (end_year > start_year):
        run_type = "transient"
    else:
        run_type = "timeslice"


    if crop_flag:
        num_pft      = "78"
    else:
        num_pft      = "16"

    logging.debug(' crop_flag = '+ crop_flag.__str__()+ ' num_pft ='+ num_pft)


    # different path for each range of years for transient cases. 
    # defualt should be picked based on the year. 1850 - 2015 -->
    # /glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629/

    #850-1850 
    #pftcftdynharv.0.25x0.25.LUH2.histsimyr0850-1849.c171012

    #based on sscp it should choose also 
    ################################################
    #sim_year         = args.sim_year
    #if args.sim_range:
    #    run_type = "transient"
    #    sim_range = args.sim_range
    #else:
    #    run_type = "timeslice"

    #if (run_type =="timeslice"):
    #    start_year = sim_year
    #    end_year   = sim_year
    #elif (run_type =="transient"):
    #    start_year, end_year = sim_range.split("-")


    build_nl (start_year, end_year, res, ssp_rcp, glc_nec, num_pft, input_path, run_type, vic_flag, glc_flag)

if __name__ == "__main__":
    main()

