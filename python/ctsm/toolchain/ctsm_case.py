# 2020-11-08                Negin Sobhani
"""
This module includes the definition for CtsmCase class for the purpose of gen_mksurf_namelist.
"""

# -- Import libraries
# -- Import Python Standard Libraries

import os
import re
import sys
import logging

from datetime import datetime

from ctsm.git_utils import tag_describe

# -- import local classes for this script
logger = logging.getLogger(__name__)

class CtsmCase:
    """
    A class to encapsulate different ctsm cases.

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
    check_num_pft:
        Determine num_pft based on crop_flag for a
        ctsm case.
    build_landuse_filename:
        Build the land-use filename for a transient
        case.
    create_landuse_file:
        Create land-use txt file for a transient case.
    build_namelist_filename
        Build the name of the namelist/control file
        for a ctsm case.
    create_namelist_file:
        Build the namelist/control file for a ctsm
        case.
    """
    # pylint: disable=too-many-instance-attributes


    def __init__(
        self,
        res,
        glc_nec,
        ssp_rcp,
        crop_flag,
        input_path,
        vic_flag,
        glc_flag,
        start_year,
        end_year,
        hres_flag
    ):
        self.res = res
        self.glc_nec = glc_nec
        self.ssp_rcp = ssp_rcp
        self.crop_flag = crop_flag
        self.input_path = input_path
        self.vic_flag = vic_flag
        self.glc_flag = glc_flag
        self.start_year = start_year
        self.end_year = end_year
        self.hres_flag = hres_flag
        self.lu_fname = None
        self.namelist_fname =None
        self.ssp_val=None
        self.rcp_val=None


        # -- check if end year value is a valid value
        self.check_endyear()

        # -- Determine if the case is transient
        self.check_run_type()

        # -- determine the num_pft
        self.check_num_pft()

    def __str__(self):
        return (
            str(self.__class__)
            + "\n"
            + "\n".join(
                (str(item) + " = " + str(self.__dict__[item]) for item in self.__dict__)
            )
        )

    def check_endyear(self):
        """
        check if end_year is valid.

        Raises:
            Error is end_year is smaller than start_year
        """
        if self.end_year < self.start_year:
            sys.exit(
                "ERROR: end_year should be bigger than the start_year : "
                + self.start_year.__str__()
                + "."
            )

    def check_run_type(self):
        """
        Determine if a ctsm case is transient or
        time-slice.
        """
        if self.end_year > self.start_year:
            self.run_type = "transient"
        else:
            self.run_type = "timeslice"
        logger.debug(" run_type  = %s", self.run_type)

    def check_num_pft(self):
        """
        determine the num_pft
        """
        if self.crop_flag:
            self.num_pft = "78"
        else:
            self.num_pft = "16"
        logger.debug(
            " crop_flag = " + self.crop_flag.__str__() + " => num_pft = " + self.num_pft
        )

    def build_landuse_filename(self):
        """
        Build the land-use filename for a transient
        case.
        """
        if self.run_type == "transient":
            lu_fname = (
                "landuse_timeseries_hist_"
                + self.num_pft.__str__()
                + "pfts_simyr"
                + self.start_year.__str__()
                + "-"
                + self.end_year.__str__()
                + ".txt"
            )
        else:
            lu_fname = ""
        self.lu_fname = lu_fname

    def create_landuse_file(self):
        """
        Create land-use txt file for a transient case.
        """
        self.build_landuse_filename()
        lu_file = open(self.lu_fname, "w")

        for yr in range(self.start_year, self.end_year + 1):

            # -- choose different files for years of 850-1850:
            if 849 < yr < 1850:
                lu_input_fname = os.path.join(
                    self.input_path,
                    "pftcftdynharv.0.25x0.25.LUH2.histsimyr0850-1849.c171012",
                    "mksrf_landuse_histclm50_LUH2_" + str(yr) + ".c171012.nc",
                )
            elif 1849 < yr < 2016:
                lu_input_fname = os.path.join(
                    self.input_path,
                    "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412",
                    "mksrf_landuse_histclm50_LUH2_" + str(yr) + ".c170412.nc",
                )
            elif 2015 < yr < 2106:
                self.decode_ssp_rcp()
                lu_input_fname = os.path.join(
                    self.input_path,
                    "pftcftdynharv.0.25x0.25."
                    + self.ssp_rcp
                    + ".simyr2016-2100.c181217",
                    "mksrf_landuse_SSP"
                    + self.ssp_val
                    + "RCP"
                    + self.rcp_val
                    + "_clm5_"
                    + str(yr)
                    + ".c181217.nc",
                )
            else:
                logger.warning("year:", yr, " not valid.")

            # -- Check if the land-use input file exist:
            if not os.path.isfile(lu_input_fname):
                logger.debug("lu_input_fname:", lu_input_fname)
                logger.warning("land-use input file does not exist for year:", yr, ".")

            # TODO: make the space/tab exactly the same as pl code:
            lu_line = lu_input_fname + "\t\t\t" + str(yr) + "\n"

            # -- Each line is written twice in the original pl code:
            lu_file.write(lu_line)
            lu_file.write(lu_line)

            logger.debug("year : %s", yr)
            logger.debug(lu_line)

        print ("Successfully created land use file : ", self.lu_fname, ".")
        print("-------------------------------------------------------")

    def build_namelist_filename(self):
        """
        Build namelist file name.
        """
        time_stamp = datetime.today().strftime("%y%m%d")
        namelist_fname = (
            "surfdata_"
            + self.res
            + "_"
            + self.ssp_rcp
            + "_"
            + self.num_pft
            + "pfts_CMIP6_"
            + self.start_year.__str__()
            + "-"
            + self.end_year.__str__()
            + "_c"
            + time_stamp
            + ".namelist"
        )

        self.namelist_fname = namelist_fname

    def create_namelist_file(self):
        """
        Build the namelist/control file for a ctsm case.
        """

        self.build_landuse_filename()
        if self.run_type == "transient":
            self.create_landuse_file()

        self.build_namelist_filename()
        namelist_file = open(self.namelist_fname, "w")

        label = tag_describe()

        dst_mesh = which_mesh(self.res)

        logger.debug("dst mesh is : %s", dst_mesh)

        if self.run_type == "transient":
            use_transient = ".true."
        else:
            use_transient = ".false"

        # pylint: disable=line-too-long

        nl_template = (
            "&clmexp\n"
            "nglcec           = " + self.glc_nec + "\n"
            "mksrf_fsoitex    = "
            + self.input_path
            + "mksrf_soitex.10level.c201018.nc"
            + "\n"
            "mksrf_forganic   = "
            + self.input_path
            + "mksrf_organic_10level_5x5min_ISRIC-WISE-NCSCD_nlev7_c120830.nc"
            + "\n"
            "mksrf_flakwat    = "
            + self.input_path
            + "mksrf_LakePnDepth_3x3min_simyr2004_csplk_c151015.nc"
            + "\n"
            "mksrf_fwetlnd    = " + self.input_path + "mksrf_lanwat.050425.nc" + "\n"
            "mksrf_fmax       = "
            + self.input_path
            + "mksrf_fmax_3x3min_USGS_c120911.nc"
            + "\n"
            "mksrf_fglacier   = "
            + self.input_path
            + "mksrf_glacier_3x3min_simyr2000.c120926.nc"
            + "\n"
            "mksrf_fvocef     = "
            + self.input_path
            + "mksrf_vocef_0.5x0.5_simyr2000.c110531.nc"
            + "\n"
            "mksrf_furbtopo   = "
            + self.input_path
            + "mksrf_topo.10min.c080912.nc"
            + "\n"
            "mksrf_fgdp       = "
            + self.input_path
            + "mksrf_gdp_0.5x0.5_AVHRR_simyr2000.c130228.nc"
            + "\n"
            "mksrf_fpeat      = "
            + self.input_path
            + "mksrf_peatf_0.5x0.5_AVHRR_simyr2000.c130228.nc"
            + "\n"
            "mksrf_fsoildepth = "
            + self.input_path
            + "mksf_soilthk_5x5min_ORNL-Soil_simyr1900-2015_c170630.nc"
            + "\n"
            "mksrf_fabm       = "
            + self.input_path
            + "mksrf_abm_0.5x0.5_AVHRR_simyr2000.c130201.nc"
            + "\n"
            "outnc_double     = .true. \n"
            "all_urban        = .false.\n"
            "no_inlandwet     = .true. \n"
            "mksrf_furban     = "
            + self.input_path
            + "mksrf_urban_0.05x0.05_simyr2000.c170724.nc"
            + "\n"
            "gitdescribe      = " + label + "\n"
            "mksrf_ftopostats = "
            + self.input_path
            + "mksrf_topostats_1km-merge-10min_HYDRO1K-merge-nomask_simyr2000.c130402.nc"
            + "\n"
            "mksrf_fvegtyp    = "
            + self.input_path
            + "pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629/mksrf_landuse_histclm50_LUH2_1850.c170629.nc"
            + "\n"
            "mksrf_fsoicol    = "
            + self.input_path
            + "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_soilcolor_CMIP6_simyr2005.c170623.nc"
            + "\n"
            "mksrf_flai       = "
            + self.input_path
            + "pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_lai_78pfts_simyr2005.c170413.nc"
            + "\n"
            "fdyndat          = ''\n"
            "numpft           = " + self.num_pft + "\n"
            "dst_mesh_file    = " + self.input_path + dst_mesh + "\n"
            "\n&transient\n"
            "use_transient    = " + use_transient + "\n"
            "start_year       = " + self.start_year.__str__() + "\n"
            "end_year         = " + self.end_year.__str__() + "\n"
            "mksrf_dyn_lu     = "
            + self.input_path
            + "pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629"
            + "\n"
            "mksrf_fdynuse    = " + self.lu_fname + "\n"
            "\n&vic\n"
            "use_vic          = " + self.vic_flag.__str__() + "\n"
            "mksrf_fvic       = "
            + self.input_path
            + "mksrf_vic_0.9x1.25_GRDC_simyr2000.c130307.nc\n"
            "outnc_vic        = \n"
            "\n&glc\n"
            "use_glc          = " + self.glc_flag.__str__() + "\n"
            "outnc_3dglc      = \n"
            "/\n"
        )
        # pylint: enable=line-too-long

        print("Successfully created namelist file : ", self.namelist_fname, ".")
        print("--------------------------------------------------------")
        namelist_file.write(nl_template)
        namelist_file.close()

    def decode_ssp_rcp(self):
        """
        Decode ssp_rcp strings because
        the raw filenames and folder names
        are not consistent.

        For example:
        folder names have ssp_rcp as SSP1-2.6

        but the files in these folders have
        ssp_rcp as SSP1RCP26
        """

        if self.ssp_rcp is not "hist":
            temp = re.sub("[^0-9]", "", self.ssp_rcp)
            self.ssp_val = temp[0]
            self.rcp_val = temp[1:]
        else:
            sys.exit(
                "ERROR: \n"
                + "\t Please choose a ssp_rcp scenario for years beyond 2015 using --ssp_rcp flag."
            )


def which_mesh(res):
    """
    Figure out the dst mesh file for each res
    """
    switcher = {
        "512x1024": "lnd/clm2/mappingdata/grids/SCRIPgrid_512x1024_nomask_c110308.nc",
        "128x256": "lnd/clm2/mappingdata/grids/SCRIPgrid_128x256_nomask_c110308.nc",
        "94x192": "lnd/clm2/mappingdata/grids/SCRIPgrid_94x192_nomask_c110308.nc",
        "64x128": "lnd/clm2/mappingdata/grids/SCRIPgrid_64x128_nomask_c110308.nc",
        "48x96": "lnd/clm2/mappingdata/grids/SCRIPgrid_48x96_nomask_c110308.nc",
        "32x64": "lnd/clm2/mappingdata/grids/SCRIPgrid_32x64_nomask_c110308.nc",
        "8x16": "lnd/clm2/mappingdata/grids/SCRIPgrid_8x16_nomask_c110308.nc",
        "0.23x0.31": "lnd/clm2/mappingdata/grids/SCRIPgrid_0.23x0.31_nomask_c110308.nc",
        "0.47x0.63": "lnd/clm2/mappingdata/grids/SCRIPgrid_0.47x0.63_nomask_c170914.nc",
        "0.9x1.25": "lnd/clm2/mappingdata/grids/0.9x1.25_c110307.nc",
        "1.9x2.5": "lnd/clm2/mappingdata/grids/1.9x2.5_c110308.nc",
        "2.5x3.33": "lnd/clm2/mappingdata/grids/SCRIPgrid_2.5x3.33_nomask_c110308.nc",
        "4x5": "lnd/clm2/mappingdata/grids/SCRIPgrid_4x5_nomask_c110308.nc",
        "10x15": "lnd/clm2/mappingdata/grids/SCRIPgrid_10x15_nomask_c110308.nc",
        "C384": "atm/cam/coords/C384_SCRIP_desc.181018.nc",
        "C192": "atm/cam/coords/C192_SCRIP_desc.181018.nc",
        "C96": "atm/cam/coords/C96_SCRIP_desc.181018.nc",
        "C48": "atm/cam/coords/C48_SCRIP_desc.181018.nc",
        "C24": "atm/cam/coords/C24_SCRIP_desc.181018.nc",
        "ne240np4": "lnd/clm2/mappingdata/grids/SCRIPgrid_ne240np4_nomask_c091227.nc",
        "ne120np4": "lnd/clm2/mappingdata/grids/SCRIPgrid_ne120np4_nomask_c101123.nc",
        "ne60np4": "lnd/clm2/mappingdata/grids/SCRIPgrid_ne60np4_nomask_c100408.nc",
        "ne30np4": "lnd/clm2/mappingdata/grids/SCRIPgrid_ne30np4_nomask_c101123.nc",
        "ne16np4": "lnd/clm2/mappingdata/grids/SCRIPgrid_ne16np4_nomask_c110512.nc",
        "360x720cru": "lnd/clm2/mappingdata/grids/SCRIPgrid_360x720_nomask_c120830.nc",
    }

    return switcher.get(res, "nothing")
