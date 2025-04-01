"""
CTSM-specific test that first performs a GDD-generating run, then calls
Python code to generate the maturity requirement file. This is then used
in a sowing+maturity forced run, which finally is tested to ensure
correct behavior.

Currently only supports 0.9x1.25, 1.9x2.5, and 10x15 resolutions. Eventually,
this test should be able to generate its own files at whatever resolution it's
called at. Well, really, the ultimate goal would be to give CLM the files
at the original resolution (for GGCMI phase 3, 0.5°) and have the stream
code do the interpolation. However, that wouldn't act on harvest dates
(which are needed for generate_gdds.py). I could have Python interpolate
those, but this would cause a potential inconsistency.
"""

import os
import re
import systemtest_utils as stu
import subprocess
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.test_utils.user_nl_utils import append_to_user_nl_files
from CIME.case import Case
import shutil, glob

logger = logging.getLogger(__name__)


class RXCROPMATURITYSHARED(SystemTestsCommon):
    def __init__(self, case):
        # initialize an object interface to the SMS system test
        SystemTestsCommon.__init__(self, case)

        # Is this a real RXCROPMATURITY test or not?
        casebaseid = self._case.get_value("CASEBASEID")
        full_test = "RXCROPMATURITY_" in casebaseid
        skipgen_test = "RXCROPMATURITYSKIPGEN_" in casebaseid

        # Ensure run length is at least 5 years. Minimum to produce one complete growing season
        # (i.e., two complete calendar years) actually 4 years, but that only gets you 1 season
        # usable for GDD generation, so you can't check for season-to-season consistency.
        stop_n = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        stop_n_orig = stop_n
        stop_option_orig = stop_option
        if "nsecond" in stop_option:
            stop_n /= 60
            stop_option = "nminutes"
        if "nminute" in stop_option:
            stop_n /= 60
            stop_option = "nhours"
        if "nhour" in stop_option:
            stop_n /= 24
            stop_option = "ndays"
        if "nday" in stop_option:
            stop_n /= 365
            stop_option = "nyears"
        if "nmonth" in stop_option:
            stop_n /= 12
            stop_option = "nyears"
        error_message = None
        if "nyear" not in stop_option:
            error_message = (
                f"STOP_OPTION ({stop_option_orig}) must be nsecond(s), nminute(s), "
                + "nhour(s), nday(s), nmonth(s), or nyear(s)"
            )
        elif full_test and stop_n < 5:
            error_message = (
                "RXCROPMATURITY must be run for at least 5 years; you requested "
                + f"{stop_n_orig} {stop_option_orig[1:]}"
            )
        elif skipgen_test and stop_n < 3:
            # First year is discarded because crops are already in the ground at restart, and those
            # aren't affected by the new crop calendar inputs. The second year is useable, but we
            # need a third year so that all crops planted in the second year have a chance to
            # finish.
            error_message = (
                "RXCROPMATURITYSKIPGEN (both-forced part) must be run for at least 3 years; you requested "
                + f"{stop_n_orig} {stop_option_orig[1:]}"
            )
        if error_message is not None:
            logger.error(error_message)
            raise RuntimeError(error_message)

        # Get the number of complete years that will be run
        self._run_Nyears = int(stop_n)

        # Only allow RXCROPMATURITY to be called with test cropMonthOutput
        if casebaseid.split("-")[-1] != "cropMonthOutput":
            error_message = (
                "Only call RXCROPMATURITY with test cropMonthOutput "
                + "to avoid potentially huge sets of daily outputs."
            )
            logger.error(error_message)
            raise RuntimeError(error_message)

        # Get files with prescribed sowing and harvest dates
        self._get_rx_dates()

        # Get cultivar maturity requirement file to fall back on if not generating it here
        self._gdds_file = None
        self._fallback_gdds_file = os.path.join(
            os.path.dirname(self._sdatefile), "gdds_20230829_161011.nc"
        )

        # Which conda environment should we use?
        self._get_conda_env()

    def _run_phase(self, skip_gen=False, h1_inst=False):
        # Modeling this after the SSP test, we create a clone to be the case whose outputs we don't
        # want to be saved as baseline.

        # -------------------------------------------------------------------
        # (1) Set up GDD-generating run
        # -------------------------------------------------------------------
        # Create clone to be GDD-Generating case
        logger.info("RXCROPMATURITY log:  cloning setup")
        case_rxboth = self._case
        caseroot = self._case.get_value("CASEROOT")
        clone_path = f"{caseroot}.gddgen"
        self._path_gddgen = clone_path
        if os.path.exists(self._path_gddgen):
            shutil.rmtree(self._path_gddgen)
        logger.info("RXCROPMATURITY log:  cloning")
        case_gddgen = self._case.create_clone(clone_path, keepexe=True)
        logger.info("RXCROPMATURITY log:  done cloning")

        os.chdir(self._path_gddgen)
        self._set_active_case(case_gddgen)

        # Set up stuff that applies to both tests
        self._setup_all(h1_inst)

        # Add stuff specific to GDD-Generating run
        logger.info("RXCROPMATURITY log:  modify user_nl files: generate GDDs")
        self._append_to_user_nl_clm(
            [
                "stream_fldFileName_cultivar_gdds = ''",
                "generate_crop_gdds = .true.",
                "use_mxmat = .false.",
                " ",
                "! (h2) Daily outputs for GDD generation and figure-making",
                "hist_fincl3 = 'GDDACCUM', 'GDDHARV'",
                "hist_nhtfrq(3) = -24",
                "hist_mfilt(3) = 365",
                "hist_type1d_pertape(3) = 'PFTS'",
                "hist_dov2xy(3) = .false.",
            ]
        )

        # If flanduse_timeseries is defined, we need to make a static version for this test. This
        # should have every crop in most of the world.
        self._get_flanduse_timeseries_in(case_gddgen)
        if self._flanduse_timeseries_in is not None:

            # Download files from the server, if needed
            case_gddgen.check_all_input_data()

            # Copy needed file from original to gddgen directory
            shutil.copyfile(
                os.path.join(caseroot, ".env_mach_specific.sh"),
                os.path.join(self._path_gddgen, ".env_mach_specific.sh"),
            )

            # Make custom version of surface file
            logger.info("RXCROPMATURITY log:  run fsurdat_modifier")
            self._run_fsurdat_modifier()

        # -------------------------------------------------------------------
        # (2) Perform GDD-generating run and generate prescribed GDDs file
        # -------------------------------------------------------------------
        logger.info("RXCROPMATURITY log:  Start GDD-Generating run")

        # As per SSP test:
        # "No history files expected, set suffix=None to avoid compare error"
        # We *do* expect history files here, but anyway. This works.
        self._skip_pnl = False

        # If not generating GDDs, only run a few days of this.
        if skip_gen:
            with Case(self._path_gddgen, read_only=False) as case:
                case.set_value("STOP_N", 5)
                case.set_value("STOP_OPTION", "ndays")

        self.run_indv(suffix=None, st_archive=True)
        if skip_gen:
            # Interpolate an existing GDD file. Needed to check obedience to GDD inputs.
            self._run_interpolate_gdds()
        else:
            self._run_generate_gdds(case_gddgen)

        # -------------------------------------------------------------------
        # (3) Set up and perform Prescribed Calendars run
        # -------------------------------------------------------------------
        os.chdir(caseroot)
        self._set_active_case(case_rxboth)

        # Set up stuff that applies to both tests
        self._setup_all(h1_inst)

        # Add stuff specific to Prescribed Calendars run
        logger.info("RXCROPMATURITY log:  modify user_nl files: Prescribed Calendars")
        self._append_to_user_nl_clm(
            [
                "generate_crop_gdds = .false.",
                f"stream_fldFileName_cultivar_gdds = '{self._gdds_file}'",
            ]
        )

        self.run_indv()

        # -------------------------------------------------------------------
        # (4) Check Prescribed Calendars run
        # -------------------------------------------------------------------
        logger.info("RXCROPMATURITY log:  output check: Prescribed Calendars")
        self._run_check_rxboth_run(skip_gen)

    # Get sowing and harvest dates for this resolution.
    def _get_rx_dates(self):
        # Eventually, I want to remove these hard-coded resolutions so that this test can generate
        # its own sowing and harvest date files at whatever resolution is requested.
        lnd_grid = self._case.get_value("LND_GRID")
        input_data_root = self._case.get_value("DIN_LOC_ROOT")
        processed_crop_dates_dir = f"{input_data_root}/lnd/clm2/cropdata/calendars/processed"
        if lnd_grid == "10x15":
            self._sdatefile = os.path.join(
                processed_crop_dates_dir,
                "sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f10_f10_mg37.2000-2000.20230330_165301.nc",
            )
            self._hdatefile = os.path.join(
                processed_crop_dates_dir,
                "hdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f10_f10_mg37.2000-2000.20230330_165301.nc",
            )
        elif lnd_grid == "1.9x2.5":
            self._sdatefile = os.path.join(
                processed_crop_dates_dir,
                "sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f19_g17.2000-2000.20230102_175625.nc",
            )
            self._hdatefile = os.path.join(
                processed_crop_dates_dir,
                "hdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f19_g17.2000-2000.20230102_175625.nc",
            )
        elif lnd_grid == "0.9x1.25":
            self._sdatefile = os.path.join(
                processed_crop_dates_dir,
                "sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f09_g17.2000-2000.20230520_134417.nc",
            )
            self._hdatefile = os.path.join(
                processed_crop_dates_dir,
                "hdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f09_g17.2000-2000.20230520_134418.nc",
            )
        else:
            error_message = "ERROR: RXCROPMATURITY currently only supports 0.9x1.25, 1.9x2.5, and 10x15 resolutions"
            logger.error(error_message)
            raise RuntimeError(error_message)

        # Ensure files exist
        error_message = None
        if not os.path.exists(self._sdatefile):
            error_message = f"ERROR: Sowing date file not found: {self._sdatefile}"
        elif not os.path.exists(self._hdatefile):
            error_message = f"ERROR: Harvest date file not found: {self._sdatefile}"
        if error_message is not None:
            logger.error(error_message)
            raise RuntimeError(error_message)

    def _setup_all(self, h1_inst):
        logger.info("RXCROPMATURITY log:  _setup_all start")

        # Get some info
        self._ctsm_root = self._case.get_value("COMP_ROOT_DIR_LND")
        run_startdate = self._case.get_value("RUN_STARTDATE")
        self._run_startyear = int(run_startdate.split("-")[0])

        # Set sowing dates file (and other crop calendar settings) for all runs
        logger.info("RXCROPMATURITY log:  modify user_nl files: all tests")
        self._modify_user_nl_allruns(h1_inst)
        logger.info("RXCROPMATURITY log:  _setup_all done")

    # Make a surface dataset that has every crop in every gridcell
    def _run_fsurdat_modifier(self):

        # fsurdat should be defined. Where is it?
        self._fsurdat_in = None
        with open(self._lnd_in_path, "r") as lnd_in:
            for line in lnd_in:
                fsurdat_in = re.match(r" *fsurdat *= *'(.*)'", line)
                if fsurdat_in:
                    self._fsurdat_in = fsurdat_in.group(1)
                    break
        if self._fsurdat_in is None:
            error_message = "fsurdat not defined"
            logger.error(error_message)
            raise RuntimeError(error_message)

        # Where we will save the fsurdat version for this test
        path, ext = os.path.splitext(self._fsurdat_in)
        dir_in, filename_in_noext = os.path.split(path)
        self._fsurdat_out = os.path.join(
            self._path_gddgen, f"{filename_in_noext}.all_crops_everywhere{ext}"
        )

        # Make fsurdat for this test, if not already done
        if not os.path.exists(self._fsurdat_out):
            tool_path = os.path.join(
                self._ctsm_root,
                "tools",
                "modify_input_files",
                "fsurdat_modifier",
            )

            # Create configuration file for fsurdat_modifier
            self._cfg_path = os.path.join(
                self._path_gddgen,
                "modify_fsurdat_allcropseverywhere.cfg",
            )
            self._create_config_file_evenlysplitcrop()

            command = f"python3 {tool_path} {self._cfg_path} "
            stu.run_python_script(
                self._get_caseroot(),
                self._this_conda_env,
                command,
                tool_path,
            )

        # Modify namelist
        logger.info("RXCROPMATURITY log:  modify user_nl files: new fsurdat")
        self._append_to_user_nl_clm(
            [
                "fsurdat = '{}'".format(self._fsurdat_out),
                "do_transient_crops = .false.",
                "flanduse_timeseries = ''",
                "use_init_interp = .true.",
            ]
        )

    def _create_config_file_evenlysplitcrop(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        cfg_template_path = os.path.join(
            self._ctsm_root, "tools/modify_input_files/modify_fsurdat_template.cfg"
        )

        with open(self._cfg_path, "w", encoding="utf-8") as cfg_out:
            # Copy template, replacing some lines
            with open(cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *evenly_split_cropland *=", line):
                        line = f"evenly_split_cropland = True"
                    elif re.match(r" *fsurdat_in *=", line):
                        line = f"fsurdat_in = {self._fsurdat_in}"
                    elif re.match(r" *fsurdat_out *=", line):
                        line = f"fsurdat_out = {self._fsurdat_out}"
                    elif re.match(r" *process_subgrid_section *=", line):
                        line = f"process_subgrid_section = True"
                    cfg_out.write(line)

            # Add new lines
            cfg_out.write("\n")
            cfg_out.write("[modify_fsurdat_subgrid_fractions]\n")
            cfg_out.write("PCT_CROP    = 100.0\n")
            cfg_out.write("PCT_NATVEG  = 0.0\n")
            cfg_out.write("PCT_GLACIER = 0.0\n")
            cfg_out.write("PCT_WETLAND = 0.0\n")
            cfg_out.write("PCT_LAKE    = 0.0\n")
            cfg_out.write("PCT_OCEAN   = 0.0\n")
            cfg_out.write("PCT_URBAN   = 0.0 0.0 0.0\n")

    def _run_check_rxboth_run(self, skip_gen):

        output_dir = os.path.join(self._get_caseroot(), "run")

        if skip_gen:
            first_usable_year = self._run_startyear + 1
            last_usable_year = first_usable_year
        else:
            first_usable_year = self._run_startyear + 2
            last_usable_year = self._run_startyear + self._run_Nyears - 2

        tool_path = os.path.join(
            self._ctsm_root, "python", "ctsm", "crop_calendars", "check_rxboth_run.py"
        )
        command = (
            f"python3 {tool_path} "
            + f"--directory {output_dir} "
            + f"-y1 {first_usable_year} "
            + f"-yN {last_usable_year} "
            + f"--rx-sdates-file {self._sdatefile} "
            + f"--rx-gdds-file {self._gdds_file} "
        )
        stu.run_python_script(
            self._get_caseroot(),
            self._this_conda_env,
            command,
            tool_path,
        )

    def _modify_user_nl_allruns(self, h1_inst):
        nl_additions = [
            "cropcals_rx = .true.",
            "cropcals_rx_adapt = .false.",
            "stream_meshfile_cropcal = '{}'".format(self._case.get_value("LND_DOMAIN_MESH")),
            "stream_fldFileName_swindow_start = '{}'".format(self._sdatefile),
            "stream_fldFileName_swindow_end   = '{}'".format(self._sdatefile),
            "stream_year_first_cropcal_swindows = 2000",
            "stream_year_last_cropcal_swindows = 2000",
            "model_year_align_cropcal_swindows = 2000",
            " ",
            "! (h1) Annual outputs on sowing or harvest axis",
            "hist_fincl2 = 'GRAINC_TO_FOOD_PERHARV', 'GRAINC_TO_FOOD_ANN', 'SDATES', 'SDATES_PERHARV', 'SYEARS_PERHARV', 'HDATES', 'GDDHARV_PERHARV', 'GDDACCUM_PERHARV', 'HUI_PERHARV', 'SOWING_REASON_PERHARV', 'HARVEST_REASON_PERHARV'",
            "hist_nhtfrq(2) = 17520",
            "hist_mfilt(2) = 999",
            "hist_type1d_pertape(2) = 'PFTS'",
            "hist_dov2xy(2) = .false.",
        ]
        if h1_inst:
            nl_additions.append("hist_avgflag_pertape(2) = 'I'")
        self._append_to_user_nl_clm(nl_additions)

    def _run_generate_gdds(self, case_gddgen):
        self._generate_gdds_dir = os.path.join(self._path_gddgen, "generate_gdds_out")
        os.makedirs(self._generate_gdds_dir)

        # Get arguments to generate_gdds.py
        dout_sr = case_gddgen.get_value("DOUT_S_ROOT")
        input_dir = os.path.join(dout_sr, "lnd", "hist")
        first_season = self._run_startyear + 2
        last_season = self._run_startyear + self._run_Nyears - 2
        sdates_file = self._sdatefile
        hdates_file = self._hdatefile

        # It'd be much nicer to call generate_gdds.main(), but I can't import generate_gdds.
        # See https://github.com/ESCOMP/CTSM/issues/2603
        tool_path = os.path.join(
            self._ctsm_root, "python", "ctsm", "crop_calendars", "generate_gdds.py"
        )
        command = " ".join(
            [
                f"python3 {tool_path}",
                f"--input-dir {input_dir}",
                f"--first-season {first_season}",
                f"--last-season {last_season}",
                f"--sdates-file {sdates_file}",
                f"--hdates-file {hdates_file}",
                f"--output-dir generate_gdds_out",
                f"--skip-crops miscanthus,irrigated_miscanthus,switchgrass,irrigated_switchgrass",
            ]
        )
        stu.run_python_script(
            self._get_caseroot(),
            self._this_conda_env,
            command,
            tool_path,
        )

        # Where were the prescribed maturity requirements saved?
        generated_gdd_files = glob.glob(os.path.join(self._generate_gdds_dir, "gdds_*.nc"))
        if len(generated_gdd_files) != 1:
            error_message = f"ERROR: Expected one matching prescribed maturity requirements file; found {len(generated_gdd_files)}: {generated_gdd_files}"
            logger.error(error_message)
            raise RuntimeError(error_message)
        self._gdds_file = generated_gdd_files[0]

    def _run_interpolate_gdds(self):
        # File where interpolated GDDs should be saved
        self._gdds_file = os.path.join(self._get_caseroot(), "interpolated_gdds.nc")

        # It'd be much nicer to call interpolate_gdds.main(), but I can't import interpolate_gdds.
        # See https://github.com/ESCOMP/CTSM/issues/2603
        tool_path = os.path.join(
            self._ctsm_root, "python", "ctsm", "crop_calendars", "interpolate_gdds.py"
        )
        command = " ".join(
            [
                f"python3 {tool_path}",
                f"--input-file {self._fallback_gdds_file}",
                f"--target-file {self._sdatefile}",
                f"--output-file {self._gdds_file}",
                "--overwrite",
            ]
        )
        stu.run_python_script(
            self._get_caseroot(),
            self._this_conda_env,
            command,
            tool_path,
        )

    def _get_conda_env(self):
        conda_setup_commands = stu.cmds_to_setup_conda(self._get_caseroot())

        # If npl conda environment is available, use that (It has dask, which
        # enables chunking, which makes reading daily 1-degree netCDF files
        # much more efficient.
        if "npl " in os.popen(conda_setup_commands + "conda env list").read():
            self._this_conda_env = "npl"
        else:
            self._this_conda_env = "ctsm_pylib"

    def _append_to_user_nl_clm(self, additions):
        caseroot = self._get_caseroot()
        append_to_user_nl_files(caseroot=caseroot, component="clm", contents=additions)

    # Is flanduse_timeseries defined? If so, where is it?
    def _get_flanduse_timeseries_in(self, case):
        case.create_namelists(component="lnd")
        self._lnd_in_path = os.path.join(self._path_gddgen, "CaseDocs", "lnd_in")
        self._flanduse_timeseries_in = None
        with open(self._lnd_in_path, "r") as lnd_in:
            for line in lnd_in:
                flanduse_timeseries_in = re.match(r" *flanduse_timeseries *= *'(.*)'", line)
                if flanduse_timeseries_in:
                    self._flanduse_timeseries_in = flanduse_timeseries_in.group(1)
                    break


class RXCROPMATURITY(RXCROPMATURITYSHARED):
    def run_phase(self):
        self._run_phase()
