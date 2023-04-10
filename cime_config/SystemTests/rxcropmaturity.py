"""
CTSM-specific test that first performs a GDD-generating run, then calls
Python code to generate the maturity requirement file. This is then used
in a sowing+maturity forced run, which finally is tested to ensure
correct behavior.
"""

import os
import re
import subprocess
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files

logger = logging.getLogger(__name__)

class RXCROPMATURITY(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)
        
        
        """
        Set run start and length
        """
        
        self._case.set_value("RUN_STARTDATE", "1995-01-01")
        self._case.set_value("STOP_OPTION", "nyears")
        self._case.set_value("STOP_N", 10)
        
        
        """
        Get and set sowing and harvest dates
        """
        
        # Get sowing and harvest dates for this resolution
        lnd_grid = self._case.get_value("LND_GRID")
        blessed_crop_dates_dir="/glade/work/samrabin/crop_dates_blessed"
        if lnd_grid == "10x15":
            self._sdatefile = os.path.join(
                blessed_crop_dates_dir,
                "sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f10_f10_mg37.2000-2000.20230330_165301.fill1.nc")
            self._hdatefile = os.path.join(
                blessed_crop_dates_dir,
                "hdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f10_f10_mg37.2000-2000.20230330_165301.fill1.nc")
        elif lnd_grid == "1.9x2.5":
            self._sdatefile = os.path.join(
                blessed_crop_dates_dir,
                "sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f19_g17.2000-2000.20230102_175625.fill1.nc")
            self._hdatefile = os.path.join(
                blessed_crop_dates_dir,
                "hdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f19_g17.2000-2000.20230102_175625.fill1.nc")
        else:
            print("ERROR: RXCROPMATURITY currently only supports 10x15 and 1.9x2.5 resolutions")
            raise
        if not os.path.exists(self._sdatefile):
            print(f"ERROR: Sowing date file not found: {self._sdatefile}")
            raise
        if not os.path.exists(self._hdatefile):
            print(f"ERROR: Harvest date file not found: {self._sdatefile}")
            raise
        
        # Set sowing dates file for run
        logger.info("  modify user_nl files: generate GDDs")
        self._modify_user_nl_gengdds()
        
        
        """
        If needed, generate a surface dataset file with no crops missing years
        """
        
        # Create out-of-the-box lnd_in to obtain fsurdat and flanduse_timeseries
        case.create_namelists(component='lnd')
        
        # Is flanduse_timeseries defined? If so, where is it? 
        lnd_in_path = os.path.join(self._get_caseroot(), 'CaseDocs', 'lnd_in')
        self._flanduse_timeseries_in = None
        with open (lnd_in_path,'r') as lnd_in:
            for line in lnd_in:
                flanduse_timeseries_in = re.match(r" *flanduse_timeseries *= *'(.*)'", line)
                if fsurdat_in:
                    self._flanduse_timeseries_in = flanduse_timeseries_in.group(1)
                    break
        
        # If flanduse_timeseries is defined, we need to make our own version for
        # this test (if we haven't already).
        if self._flanduse_timeseries_in is not None:
            
            # Download files from the server
            case.check_all_input_data()
            
            # fsurdat should be defined. Where is it?
            self._fsurdat_in = None
            with open (lnd_in_path,'r') as lnd_in:
                for line in lnd_in:
                    fsurdat_in = re.match(r" *fsurdat *= *'(.*)'", line)
                    if fsurdat_in:
                        self._fsurdat_in = fsurdat_in.group(1)
                        break
            if self._fsurdat_in is None:
                print("fsurdat not defined")
                raise
            
            # Where we will save the fsurdat version for this test
            self._fsurdat_out = os.path.join(self._get_caseroot(), 'fsurdat.nc')
            
            # Make fsurdat for this test, if not already done
            if not os.path.exists(self._fsurdat_out):
                logger.info("  run make_surface_for_gdden")
                self._run_make_surface_for_gdden()
            
            # Modify namelist
            logger.info("  modify user_nl files: new fsurdat")
            self._modify_user_nl_newfsurdat()
            
            
    def _run_make_surface_for_gdden(self):
        tool_path = os.path.join(self._ctsm_root,
                                 'python', 'ctsm', 'crop_calendars',
                                 'make_surface_for_gddgen.py')

        self._case.load_env(reset=True)
        conda_env = ". "+self._get_caseroot()+"/.env_mach_specific.sh; "
        # Preprend the commands to get the conda environment for python first
        conda_env += self._get_conda_env()
        # Source the env
        try:
            command = f"{conda_env}python3 {tool_path} "\
                + "--flanduse-timeseries {self._flanduse_timeseries_in} "\
                + "--fsurdat {self._fsurdat_in} "\
                + "--outfile {self._fsurdat_out}"
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            print("ERROR while getting the conda environment and/or ")
            print("running the fsurdat_modifier tool: ")
            print("(1) If your ctsm_pylib environment is out of date or you ")
            print("have not created the ctsm_pylib environment, yet, you may ")
            print("get past this error by running ./py_env_create ")
            print("in your ctsm directory and trying this test again. ")
            print("(2) If conda is not available, install and load conda, ")
            print("run ./py_env_create, and then try this test again. ")
            print("(3) If (1) and (2) are not the issue, then you may be ")
            print("getting an error within the fsurdat_modifier tool itself. ")
            print("Default error message: ")
            print(error.output)
            raise
        except:
            print("ERROR trying to run fsurdat_modifier tool.")
            raise

    def _modify_user_nl_gengdds(self):
        nl_additions = [
            "stream_meshfile_cropcal = '{}'".format(self._case.get_value("LND_DOMAIN_MESH")),
            "stream_fldFileName_sdate = '{}'".format(self._sdatefile),
            "generate_crop_gdds = .true.",
            "use_mxmat = .false.",
            "stream_year_first_cropcal = 2000",
            "stream_year_last_cropcal = 2000",
            "model_year_align_cropcal = 2000",
            " ",
            "! (h1) Daily outputs for GDD generation and figure-making",
            "hist_fincl2 = 'HUI', 'GDDACCUM', 'GDDHARV'",
            "hist_nhtfrq(2) = -24",
            "hist_mfilt(2) = 365",
            "hist_type1d_pertape(2) = 'PFTS'",
            "hist_dov2xy(2) = .false.",
            " ",
            "! (h2) Annual outputs for GDD generation (checks)",
            "hist_fincl3 = 'GRAINC_TO_FOOD_PERHARV', 'GRAINC_TO_FOOD_ANN', 'SDATES', 'SDATES_PERHARV', 'SYEARS_PERHARV', 'HDATES', 'GDDHARV_PERHARV', 'GDDACCUM_PERHARV', 'HUI_PERHARV', 'SOWING_REASON_PERHARV', 'HARVEST_REASON_PERHARV'",
            "hist_nhtfrq(3) = 17520",
            "hist_mfilt(3) = 999",
            "hist_type1d_pertape(3) = 'PFTS'",
            "hist_dov2xy(3) = .false.",
        ]
        for addition in nl_additions:
            append_to_user_nl_files(caseroot = self._get_caseroot(),
                                    component = "clm",
                                    contents = addition)

    def _modify_user_nl_newfsurdat(self):
        nl_additions = [
            "fsurdat = '{}'".format(self._fsurdat_out),
            "do_transient_crops = .false.",
            "flanduse_timeseries = ''",
            "use_init_interp = .true.",
        ]
        for addition in nl_additions:
            append_to_user_nl_files(caseroot = self._get_caseroot(),
                                    component = "clm",
                                    contents = addition)

    def _get_conda_env(self):
        #
        # Add specific commands needed on different machines to get conda available
        # Use semicolon here since it's OK to fail
        #
        # Execute the module unload/load when "which conda" fails
        # eg on cheyenne
        try:
            subprocess.run( "which conda", shell=True, check=True)
            conda_env = " "
        except subprocess.CalledProcessError:
            # Remove python and add conda to environment for cheyennne
            conda_env = "module unload python; module load conda;"

        # Activate the python environment
        conda_env += " conda activate ctsm_pylib"
        # End above to get to actual command
        conda_env += " && "

        return( conda_env )
