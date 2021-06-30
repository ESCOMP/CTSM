#!/usr/bin/env python

# 2020-11-08                Negin Sobhani

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This Python script is part of the simplified toolchain for creating
the surface dataset for ctsm cases.
This script should be used as the first step of the new toolchain. 
It will automatically creating namelist (control  file) that is 
needed for creating surface dataset and other relevant files for
running CTSM cases. 
For transient cases, it will also create a txt file that includes the
landuse files for every year. 

-------------------------------------------------------------------
Instructions for running on Cheyenne/Casper:

load the following into your local environment:

    module load python
    ncar_pylib
-------------------------------------------------------------------
To see the available options:
    ./gen_mksurf_namelist.py --help

To run the script:
    ./gen_mksurf_namelist.py
 
To remove NPL from your environment on Cheyenne/Casper:
    deactivate
-------------------------------------------------------------------
"""

#TODO (NS)

# -[x] Add default values in the help page.
# -[x] Add info for help page note for end_year -- by default is start_year
# -[ ] Possibly remove year --years and range options
#      Currently comment them out.

# -[ ] maybe a verbose option and removing debug
# -[x] --debug mode is not working...

# -[ ] add error check for hi-res and years if they are 1850 and 2005.

# -[x] different path for each range of years for transient cases. 
#      default should be picked based on the year. 1850 - 2015 -->
#       /glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629/
#      850-1850 --> 
#       pftcftdynharv.0.25x0.25.LUH2.histsimyr0850-1849.c171012

#QUESTIONS:
# -[ ] Do we need --crop to accept y/n or just --crop is enough?
# -[ ] Add information about if this is optional.
# Everything is optional because they have defaults...

#  Import libraries
from __future__ import print_function

import os
import re
import sys
import logging
import argparse
import subprocess
import tqdm

from datetime import datetime

## valid options for resolution and SSP scenarios:
valid_opts = {
        'res' :
        ['512x1024','360x720cru','128x256','64x128','48x96','94x192','0.23x0.31','0.47x0.63','0.9x1.25','1.9x2.5','2.5x3.33',
        '4x5','10x15','0.125nldas2','5x5_amazon','1x1_camdenNJ','1x1_vancouverCAN','1x1_mexicocityMEX',
        '1x1_asphaltjungleNJ','1x1_brazil,1x1_urbanc_alpha','1x1_numaIA,1x1_smallvilleIA','0.1x0.1','0.25x0.25','0.5x0.5',
        '3x3min','5x5min','10x10min','0.33x0.33','0.125x0.125','ne4np4,ne16np4','ne30np4.pg2','ne30np4.pg3','ne30np4','ne60np4','ne120np4']
        ,'ssp_rcp': ["hist","SSP1-2.6","SSP3-7.0","SSP5-3.4","SSP2-4.5","SSP1-1.9","SSP4-3.4","SSP4-6.0","SSP5-8.5"]
        }


def get_parser():
        """
            Get parser object for this script.
        """
        parser = argparse.ArgumentParser(description=__doc__,
                           formatter_class=argparse.RawDescriptionHelpFormatter)

        parser.print_usage = parser.print_help


        parser.add_argument('--sy','--start_year',
                    help='Simulation start year. [default: %(default)s] ', 
                    action="store",
                    dest="start_year",
                    required=False,
                    type = start_year_type,
                    default=2000)
        parser.add_argument('--ey','--end_year',
                    help='Simulation end year.  [default: start_year] ', 
                    action="store",
                    dest="end_year",
                    required=False,
                    type = int)
                    #type = end_year_type,
                    #default="start_year")

        # We decided to use start year and end year instead of years or range:
        # These lines should be ultimately removed if everyone agrees. 
        #parser.add_argument('-y','--year',
        #            help='Simulation year to run over.', 
        #            action="store",
        #            dest="sim_year",
        #            required=False,
        #            choices = ['1885','1895','1980','1982','2000'],
        #            default='2000')
        #parser.add_argument('--range',
        #            help='Simulation years to run over.', 
        #            action="store",
        #            dest="sim_range",
        #            required=False,
        #            choices = ['1850-2000','1850-2005','1850-2100'])

        parser.add_argument('--glc_nec',
                    help='''
                    Number of glacier elevation classes to use. 
                    [default: %(default)s] 
                    ''' ,
                    action="store",
                    dest="glc_nec",
                    type = glc_nec_type,
                    default = "10")
        parser.add_argument('--rundir', 
                    help='''
                    Directory to run in. 
                    [default: %(default)s] 
                    ''' ,
                    action="store",
                    dest="run_dir", 
                    required = False, 
                    default =os.getcwd())
        parser.add_argument('--ssp_rcp',
                    help='''
                    Shared Socioeconomic Pathway and Representative
                    Concentration Pathway Scenario name(s).
                    [default: %(default)s]
                    ''' ,
                    action="store",
                    dest="ssp_rcp",
                    required = False,
                    choices=valid_opts['ssp_rcp'],
                    default = "hist")

        ##############################################
        # In mksurfdata.pl these options are -l --dinlc
        # But the group decided --raw_dir is more descriptive.
        # If everyone agrees, the commented out line should be removed. 
        #parser.add_argument('-l','--dinlc',  #--raw_dir or --rawdata_dir

        parser.add_argument('--raw_dir','--rawdata_dir',
                    help='''
                    /path/of/root/of/input/data', 
                    [default: %(default)s]
                    ''',
                    action="store",
                    dest="input_path",
                    default="/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/")
        parser.add_argument('-d','--debug', 
                    help='Debug mode will print more information. ', 
                    action="store_true", 
                    dest="debug", 
                    default=False)
        parser.add_argument('-v','--verbose', 
                    help='Verbose mode will print more information. ', 
                    action="store_true", 
                    dest="verbose", 
                    default=False)
        parser.add_argument('--vic', 
                    help='''
                    Add the fields required for the VIC model.
                    [default: %(default)s]
                    ''', 
                    action="store_true", 
                    dest="vic_flag", 
                    default=False)
        parser.add_argument('--glc', 
                    help='''
                    Add the optional 3D glacier fields for verification of the glacier model.
                    [default: %(default)s]
                    ''', 
                    action="store_true", 
                    dest="glc_flag", 
                    default=False)
        parser.add_argument('--hirespft', 
                    help='''
                    If you want to use the high-resolution pft dataset rather
                    than the default lower resolution dataset.
                    (Low resolution is at quarter-degree, high resolution at 3-minute)
                    [Note: hires only available for 1850 and 2005.]
                    ''', 
                    action="store_true", 
                    dest="hres_flag", 
                    default=False)
        parser.add_argument('--crop', 
                    help='''
                    Create datasets with the extensive list of prognostic crop types.
                    [default: %(default)s]
                    ''', 
                    action="store",
                    type = str2bool,  
                    nargs = '?',
                    const = True,
                    required = False,
                    dest="crop_flag", 
                    default=False)
        parser.add_argument('-f','--fast', 
                    help='Toggle fast mode which does not user the large mapping file', 
                    action="store_true", 
                    dest="fast_flag", 
                    default=False)
        parser.add_argument('-r','--res',
                    help='''
                    Resolution is the supported resolution(s) to use for files.
                    [default: %(default)s]
                    ''',
                    action="store", 
                    dest="res", 
                    choices=valid_opts['res'], 
                    required=False, 
                    default="4x5")
        return parser


def str2bool(v):
    """
    Function for converting different forms of
    command line boolean strings to boolean value.

    Args:
        v (str): String bool input

    Raises:
        if the argument is not an acceptable boolean string
        (such as yes or no ; true or false ; y or n ; t or f ; 0 or 1).
        argparse.ArgumentTypeError: The string should be one of the mentioned values.

    Returns:
        bool: Boolean value corresponding to the input.
    """
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected. [true or false] or [y or n]')

def tag_describe ():
    """
    Function for converting different forms of
    command line boolean strings to boolean value.

    Args:

    Raises:

    Returns:
        label.decode (str) : ouput of running 'git describe' in shell
    """
    label = subprocess.check_output(["git", "describe"]).strip()
    return label.decode()

def glc_nec_type(x):
    """
    Function for defining acceptable glc_nec input.

    Args:
        x (str) : glc_nec value from command line args.

    Raises:
        Error if value of glc_nec is not in the range 
        of 1-99. 

    Returns:
        x (int) : Acceptable glc_nec value.
    """
    x = int(x)
    if (x <= 0) or (x >= 100):
        raise argparse.ArgumentTypeError("ERROR: glc_nec must be between 1 and 99.")
    return x

def start_year_type(x):
    """
    Function for defining acceptable start_year input.

    Args:
        x (str) : start_year string from command line args.

    Raises:
        Error if value of glc_start_year is not in the range 
        of 850-2015. 

    Returns:
        x (int) : Acceptable start_year value.
    """
    x = int(x)
    if (x < 850) or (x > 2105):
        raise argparse.ArgumentTypeError(
                "ERROR: Simulation start year should be between 850 and 2105.")
    return x

class CtsmCase:
    """
    A class for encapsulate different ctsm cases.

    ...

    Attributes
    ---------
    res : str
        resolution from a list of acceptable options.
    glc_nec : str
        number of glacier elevation classes.  
    ssp_rcp : str
        Shared Socioeconomic Pathway and Representative
        Concentration Pathway Scenario name.
    crop_flag : bool 
        Crop flag for determining number of pfts
    input_path : str
        Raw data input path
    vic_flag : bool
        Flag for VIC model
    glc_flag : bool
        Flag for 3D glacier model
    start_year : str
        Simulation start year
    end_year : str
        Simulation end year 

    Methods
    -------
    check_endyear:
        Check if end_year is bigger than start year
        in a ctsm case.
    check_run_type:
        Determine if a ctsm case is transient or
        time-slice. 
    def check_run_type:
        Determine num_pft based on crop_flag for a
        ctsm case.
    landuse_filename:
        Build the land-use filename for a transient 
        case.
    create_landuse:
        Create land-use txt file a transient case. 
    name_nl
        Build the name of the namelist/control file
        for a ctsm case.
    build_nl:
        Build the namelist/control file for a ctsm
        case.
    """

    def __init__ (self, res, glc_nec, ssp_rcp, crop_flag, input_path, vic_flag, glc_flag, start_year, end_year=None):
        self.res = res
        self.glc_nec = glc_nec
        self.ssp_rcp = ssp_rcp
        self.crop_flag = crop_flag
        self.input_path = input_path
        self.vic_flag = vic_flag
        self.glc_flag = glc_flag
        self.start_year = start_year
        self.end_year = end_year if end_year is not None else start_year

        #-- check if end year value is a valid value
        self.check_endyear()

        #-- Determine if the case is transient
        self.check_run_type()

        #-- determine the num_pft
        self.check_num_pft()
    

    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item])
                    for item in self.__dict__))

    def check_endyear (self):
        if (int(self.end_year) < int(self.start_year)):
            sys.exit('ERROR: end_year should be bigger than the start_year : '
                     + self.start_year.__str__() + '.')

    def check_run_type (self):
        if (self.end_year > self.start_year):
            self.run_type = "transient"
        else:               
            self.run_type = "timeslice"
        logging.debug(' run_type  = '+ self.run_type)

    def check_num_pft (self):
        #-- determine the num_pft
        if self.crop_flag:
            self.num_pft      = "78"
        else:
            self.num_pft      = "16"
        logging.debug(' crop_flag = '+ self.crop_flag.__str__()+ ' => num_pft ='+ self.num_pft)


    def name_nl  (self):
        """
        Build namelist file name.
        """
        time_stamp = datetime.today().strftime('%y%m%d')
        namelist_fname = "surfdata_"+ \
            self.res+"_"+ \
            self.ssp_rcp+"_"+ \
            self.num_pft+ \
            "pfts_CMIP6_"+ \
            self.start_year.__str__()+'-' + \
            self.end_year.__str__()+ \
            "_c"+time_stamp+".namelist"

        self.namelist_fname  = namelist_fname

    def landuse_filename(self):
        if (self.run_type == 'transient'):
            lu_fname = "landuse_timeseries_hist_78pfts_simyr"+str(self.start_year)+"-"+str(self.end_year)+".txt"
        else:
            lu_fname = ""
        self.lu_fname = lu_fname

    def create_landuse(self):
        self.landuse_filename()
        lu_file = open (self.lu_fname,'w')

        for yr in tqdm.tqdm(range(self.start_year, self.end_year+1)):

            #-- choose different files for years of 850-1850:
            if (849 < yr < 1850):
                lu_input_fname = os.path.join(self.input_path, 
                            "pftcftdynharv.0.25x0.25.LUH2.histsimyr0850-1849.c171012",
                            "mksrf_landuse_histclm50_LUH2_"+str(yr)+".c171012.nc"
                            )
            else : 
                lu_input_fname = os.path.join(self.input_path, 
                            "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412",
                            "mksrf_landuse_histclm50_LUH2_"+str(yr)+".c170412.nc"
                            )
            #-- Check if the land-use input file exist:
            if not os.path.isfile(lu_input_fname):
                print ("lu_input_fname:", lu_input_fname)
                print ("Error: land-use input file does not exist for year:", yr,".")

            # TODO: make the space/tab exactly the same as pl code:
            lu_line = lu_input_fname + "\t\t\t"+ str(yr) + "\n"
            lu_file.write (lu_line)
        print ("Successfully created land use file : ", self.lu_fname,".")
        print ("-------------------------------------------------------")


    def build_nl (self):
        """
        Build the namelist/control file for a ctsm class.  
        """

        self.landuse_filename()
        if (self.run_type == "transient"):
            self.create_landuse()

        self.name_nl()
        namelist_file = open (self.namelist_fname,'w')

        label = tag_describe()
        
        if (self.run_type == "transient"):
            use_transient = ".true."
        else: 
            use_transient = ".false"
       
        dst_mesh = which_mesh (self.res)

        print ('dst mesh is :', dst_mesh)

        nl_template = ( \
                "&clmexp\n"                                                                                                                                                                  
                "nglcec           = "+self.glc_nec + "\n"
                "mksrf_fsoitex    = "+self.input_path+"mksrf_soitex.10level.c201018.nc"+"\n"
                "mksrf_forganic   = "+self.input_path+"mksrf_organic_10level_5x5min_ISRIC-WISE-NCSCD_nlev7_c120830.nc"+"\n"
                "mksrf_flakwat    = "+self.input_path+"mksrf_LakePnDepth_3x3min_simyr2004_csplk_c151015.nc"+"\n"
                "mksrf_fwetlnd    = "+self.input_path+"mksrf_lanwat.050425.nc"+"\n"
                "mksrf_fmax       = "+self.input_path+"mksrf_fmax_3x3min_USGS_c120911.nc"+"\n"
                "mksrf_fglacier   = "+self.input_path+"mksrf_glacier_3x3min_simyr2000.c120926.nc"+"\n"
                "mksrf_fvocef     = "+self.input_path+"mksrf_vocef_0.5x0.5_simyr2000.c110531.nc" +"\n"
                "mksrf_furbtopo   = "+self.input_path+"mksrf_topo.10min.c080912.nc"+"\n"
                "mksrf_fgdp       = "+self.input_path+"mksrf_gdp_0.5x0.5_AVHRR_simyr2000.c130228.nc"+"\n"
                "mksrf_fpeat      = "+self.input_path+"mksrf_peatf_0.5x0.5_AVHRR_simyr2000.c130228.nc"+"\n"
                "mksrf_fsoildepth = "+self.input_path+"mksf_soilthk_5x5min_ORNL-Soil_simyr1900-2015_c170630.nc"+"\n"
                "mksrf_fabm       = "+self.input_path+"mksrf_abm_0.5x0.5_AVHRR_simyr2000.c130201.nc"+"\n"
                "outnc_double     = .true. \n"
                "all_urban        = .false.\n"
                "no_inlandwet     = .true. \n"
                "mksrf_furban     = "+self.input_path+"mksrf_urban_0.05x0.05_simyr2000.c170724.nc"+"\n"
                "gitdescribe      = "+label+"\n"
                "mksrf_ftopostats = "+self.input_path+"mksrf_topostats_1km-merge-10min_HYDRO1K-merge-nomask_simyr2000.c130402.nc"+"\n"
                "mksrf_fvegtyp    = "+self.input_path + "pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629/mksrf_landuse_histclm50_LUH2_1850.c170629.nc"+"\n"
                "mksrf_fsoicol    = "+self.input_path + "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_soilcolor_CMIP6_simyr2005.c170623.nc"+"\n"
                "mksrf_flai       = "+self.input_path + "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_lai_78pfts_simyr2005.c170413.nc" +"\n"
                "fdyndat          = ''\n"
                "numpft           = "+self.num_pft+"\n"
                "dst_mesh_file    = "+self.input_path+dst_mesh+"\n"
                "\n&transient\n"
                "use_transient    = "+use_transient + "\n"
                "start_year       = "+self.start_year.__str__() + "\n"
                "end_year         = "+self.end_year.__str__() + "\n"
                "mksrf_dyn_lu     = "+self.input_path+ "pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629" + "\n"
                "mksrf_fdynuse    = "+self.lu_fname + "\n"
                "\n&vic\n"
                "use_vic          = "+self.vic_flag.__str__()+"\n"
                "mksrf_fvic       = "+self.input_path+"mksrf_vic_0.9x1.25_GRDC_simyr2000.c130307.nc\n"
                "outnc_vic        = \n"
                "\n&glc\n"
                #"use_glc          = "+self.glc_flag.__str__()+"\n"
                "outnc_3dglc      = \n"
                "/\n"
                )

        print ("Successfully created namelist file : ", self.namelist_fname,".")
        print ("--------------------------------------------------------")
        namelist_file.write (nl_template)
        namelist_file.close ()

def which_mesh(res):
   '''
   Figure out the dst mesh file for each res
   '''
   switcher = { 
       "512x1024" : "lnd/clm2/mappingdata/grids/SCRIPgrid_512x1024_nomask_c110308.nc" , 
       "128x256" : "lnd/clm2/mappingdata/grids/SCRIPgrid_128x256_nomask_c110308.nc" ,
       "94x192" : "lnd/clm2/mappingdata/grids/SCRIPgrid_94x192_nomask_c110308.nc" ,
       "64x128" : "lnd/clm2/mappingdata/grids/SCRIPgrid_64x128_nomask_c110308.nc" ,
       "48x96"  : "lnd/clm2/mappingdata/grids/SCRIPgrid_48x96_nomask_c110308.nc" ,
       "32x64"  : "lnd/clm2/mappingdata/grids/SCRIPgrid_32x64_nomask_c110308.nc",
       "8x16"   : "lnd/clm2/mappingdata/grids/SCRIPgrid_8x16_nomask_c110308.nc" , 
       "0.23x0.31" : "lnd/clm2/mappingdata/grids/SCRIPgrid_0.23x0.31_nomask_c110308.nc",
       "0.47x0.63" : "lnd/clm2/mappingdata/grids/SCRIPgrid_0.47x0.63_nomask_c170914.nc",
       "0.9x1.25" : "lnd/clm2/mappingdata/grids/0.9x1.25_c110307.nc",
       "1.9x2.5" : "lnd/clm2/mappingdata/grids/1.9x2.5_c110308.nc",
       "2.5x3.33" : "lnd/clm2/mappingdata/grids/SCRIPgrid_2.5x3.33_nomask_c110308.nc",
       "4x5" : "lnd/clm2/mappingdata/grids/SCRIPgrid_4x5_nomask_c110308.nc",
       "10x15" : "lnd/clm2/mappingdata/grids/SCRIPgrid_10x15_nomask_c110308.nc",
       "C384" : "atm/cam/coords/C384_SCRIP_desc.181018.nc",
       "C192" : "atm/cam/coords/C192_SCRIP_desc.181018.nc",
       "C96" : "atm/cam/coords/C96_SCRIP_desc.181018.nc" ,
       "C48" : "atm/cam/coords/C48_SCRIP_desc.181018.nc" ,
       "C24" : "atm/cam/coords/C24_SCRIP_desc.181018.nc" ,
       "ne240np4" : "lnd/clm2/mappingdata/grids/SCRIPgrid_ne240np4_nomask_c091227.nc",
       "ne120np4" : "lnd/clm2/mappingdata/grids/SCRIPgrid_ne120np4_nomask_c101123.nc",
       "ne60np4" : "lnd/clm2/mappingdata/grids/SCRIPgrid_ne60np4_nomask_c100408.nc",
       "ne30np4" : "lnd/clm2/mappingdata/grids/SCRIPgrid_ne30np4_nomask_c101123.nc",
       "ne16np4" : "lnd/clm2/mappingdata/grids/SCRIPgrid_ne16np4_nomask_c110512.nc",
       "360x720cru" : "lnd/clm2/mappingdata/grids/SCRIPgrid_360x720_nomask_c120830.nc", 
   }

   return switcher.get(res, "nothing") 

def main ():

    args         = get_parser().parse_args()

    if args.debug or args.verbose:
        logging.basicConfig(level=logging.DEBUG)


    res          = args.res
    glc_nec      = args.glc_nec.__str__()
    input_path   = args.input_path
    ssp_rcp      = args.ssp_rcp
    crop_flag    = args.crop_flag
    vic_flag     = args.vic_flag
    glc_flag     = args.glc_flag
    hres_flag    = args.hres_flag

    start_year   = args.start_year
    end_year     = args.end_year

    #-- determine end_year if not given as an argument:
    if not end_year:
        end_year = start_year

    #-- check if the input path exist
    if not os.path.exists(input_path):
        sys.exit('ERROR: \n'+
                 '\t raw_dir does not exist on this machine. \n'+
                 '\t Please point to the correct raw_dir using --raw_dir'+
                 'or --rawdata_dir flags.')

    ctsm_case = CtsmCase(res, glc_nec, ssp_rcp, crop_flag, input_path,
                         vic_flag, glc_flag, start_year, end_year)

    ctsm_case.name_nl()

    logging.debug('--------------------------')
    logging.debug(' ctsm case : %s', ctsm_case)
    logging.debug('--------------------------')

    ctsm_case.build_nl()

if __name__ == "__main__":
    main()

