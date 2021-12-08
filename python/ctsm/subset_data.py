#!/usr/bin/env python3
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
Instructions for running on Cheyenne/Casper:
load the following into your local environment
    module load python
    ncar_pylib
-------------------------------------------------------------------
To see the available options for single point cases:
    ./subset_data.py point --help
To see the available options for regional cases:
    ./subset_data.py reg --help
-------------------------------------------------------------------
This script extracts domain files, surface dataset, and DATM files
at either a single point or a region using a global dataset. Currently this
script subsets default surface, landuse, and DATM files, which can be seen in
the defaults.cfg file.

To run a single-point or regional case using this data, you must update the
variable(s) `fsurdat` and/or `landuse` in the user_nl_clm namelist file to be
the full path to the subset files. This script will automatically create this
file using the flag --create-user-mods.
To use subset climate data, the namelist file user_nl_datm_streams must also
be updated - this script will automatically create this file with
--create-user-mods. This flag will also create necessary single-point xml
commands in the file shell_commands.

To use the created user mods with a case use --user-mods-dir PATH/TO/USER/MODS
in the ./create.newcase call.

By default, this script only extracts surface dataset. For extracting other
files, the appropriate flags should be used.
-------------------------------------------------------------------
"""

# TODO
# Automatic downloading of missing files if they are missing
# default 78 pft vs 16 pft

#  Import libraries
from __future__ import print_function

import sys
import os
import string
import logging
import subprocess
import argparse
import configparser

import numpy as np
import xarray as xr

from datetime import date
from getpass import getuser
from logging.handlers import RotatingFileHandler
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..","..",'python'))
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root


from ctsm.site_and_regional.base_case import BaseCase
from ctsm.site_and_regional.single_point_case import SinglePointCase
from ctsm.site_and_regional.regional_case import RegionalCase

# -- Globals and Default Values ---
DEFAULTS_FILE = "default_data.cfg"
myname = getuser()


def get_parser():
    """
    Get parser object for this script.
    """
    parser = ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help
    subparsers = parser.add_subparsers(
        help="Two possible ways to run this script, either:", dest="run_type")
    pt_parser = subparsers.add_parser(
        "point", help="Run script for a single point.")
    rg_parser = subparsers.add_parser("reg", help="Run script for a region.")

    # First add arguments specific to a point or regional parser
    pt_parser.add_argument(
        "--lat",
        help="Single point latitude. [default: %(default)s]",
        action="store",
        dest="plat",
        required=False,
        type=plat_type,
        default=42.5,
    )
    pt_parser.add_argument(
        "--lon",
        help="Single point longitude. [default: %(default)s]",
        action="store",
        dest="plon",
        required=False,
        type=plon_type,
        default=287.8,
    )
    pt_parser.add_argument(
        "--site",
        help="Site name or tag. [default: %(default)s]",
        action="store",
        dest="site_name",
        required=False,
        type=str,
        default="single_point",
    )
    rg_parser.add_argument(
        "--lat1",
        help="Region start latitude. [default: %(default)s]",
        action="store",
        dest="lat1",
        required=False,
        type=plat_type,
        default=-40,
    )
    rg_parser.add_argument(
        "--lat2",
        help="Region end latitude. [default: %(default)s]",
        action="store",
        dest="lat2",
        required=False,
        type=plat_type,
        default=15,
    )
    rg_parser.add_argument(
        "--lon1",
        help="Region start longitude. [default: %(default)s]",
        action="store",
        dest="lon1",
        required=False,
        type=plon_type,
        default=275.0,
    )
    rg_parser.add_argument(
        "--lon2",
        help="Region end longitude. [default: %(default)s]",
        action="store",
        dest="lon2",
        required=False,
        type=plon_type,
        default=330.0,
    )
    rg_parser.add_argument(
        "--reg",
        help="Region name or tag. [default: %(default)s]",
        action="store",
        dest="reg_name",
        required=False,
        type=str,
        default="region",
    )

    # Now add arguments shared between pt_parser and rg_parser
    for subparser in [pt_parser, rg_parser]:
        subparser.add_argument(
            "--create-domain",
            help="Create CLM domain file at single point or region. [default: %(default)s]",
            action="store",
            dest="create_domain",
            type=str2bool,
            nargs="?",
            const=True,
            required=False,
            default=False,
        )
        subparser.add_argument(
            "--create-surface",
            help="Create surface data file at single point or region. [default: %(default)s]",
            action="store",
            dest="create_surfdata",
            type=str2bool,
            nargs="?",
            const=True,
            required=False,
            default=True,
        )
        subparser.add_argument(
            "--create-landuse",
            help="Create landuse data file at single point or region. [default: %(default)s]",
            action="store",
            dest="create_landuse",
            type=str2bool,
            nargs="?",
            const=True,
            required=False,
            default=False,
        )
        subparser.add_argument(
            "--create-datm",
            help="Create DATM forcing data at single point or region. [default: %(default)s]",
            action="store",
            dest="create_datm",
            type=str2bool,
            nargs="?",
            const=True,
            required=False,
            default=False,
        )
        subparser.add_argument(
            "--datm-syr",
            help="Start year for creating DATM forcing. [default: %(default)s]",
            action="store",
            dest="datm_syr",
            required=False,
            type=int,
            default=1901,
        )
        subparser.add_argument(
            "--datm-eyr",
            help="End year for creating DATM forcing. [default: %(default)s]",
            action="store",
            dest="datm_eyr",
            required=False,
            type=int,
            default=2014,
        )
        subparser.add_argument(
            "--crop",
            help="Create datasets using the extensive list of prognostic crop types. [default: %(default)s]",
            action="store_true",
            dest="crop_flag",
            default=False,
        )
        subparser.add_argument(
            "--dompft",
            help="Dominant PFT type. [default: %(default)s] ",
            action="store",
            dest="dom_pft",
            type=int,
            default=7,
        )
        subparser.add_argument(
            "--no-unisnow",
            help="Turn off the flag for create uniform snowpack. [default: %(default)s]",
            action="store_false",
            dest="uni_snow",
            default=True,
        )
        subparser.add_argument(
            "--no-overwrite-single-pft",
            help="Turn off the flag for making the whole grid 100%% single PFT. [default: %(default)s]",
            action="store_false",
            dest="overwrite_single_pft",
            default=True,
        )
        subparser.add_argument(
            "--zero-nonveg",
            help="Set all non-vegetation landunits to zero. [default: %(default)s]",
            action="store",
            dest="zero_nonveg",
            type=bool,
            default=True,
        )
        subparser.add_argument(
            "--no-saturation-excess",
            help="Turn off the flag for saturation excess. [default: %(default)s]",
            action="store",
            dest="no_saturation_excess",
            type=bool,
            default=True,
        )
        subparser.add_argument(
            "--create-user-mods",
            help="Create a user mods directory for running CTSM. [default: %(default)s]",
            action="store",
            dest="create_user_mods",
            type=str2bool,
            default=False,

        )
        subparser.add_argument(
            "--user-mods-dir",
            help="Path to user mods directory. [default: %(default)s]",
            action="store",
            dest="user_mods_dir",
            type=str,
            default="/glade/scratch/" + myname + "/subset_data/user_mods",
        )
        subparser.add_argument(
            "--outdir",
            help="Output directory. [default: %(default)s]",
            action="store",
            dest="out_dir",
            type=str,
            default="/glade/scratch/" + myname + "/subset_data/",
        )

    newline = "\n"
    parser.epilog = f"""==================================={newline}{newline}{pt_parser.format_help()}{newline}{newline}{rg_parser.format_help()}"""

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
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError(
            "Boolean value expected. [true or false] or [y or n]"
        )


def plat_type(x):
    """
    Function to define lat type for the parser
    and
    raise error if latitude is not between -90 and 90.
    Args:
        x(str): latitude
    Raises:
        Error when x (latitude) is not between -90 and 90.
    Returns:
        x (float): latitude in float
    """
    x = float(x)
    if (x < -90) or (x > 90):
        raise argparse.ArgumentTypeError(
            "ERROR: Latitude should be between -90 and 90."
        )
    return x


def plon_type(x):
    """
    Function to define lon type for the parser and
    convert negative longitudes and
    raise error if lon is not between -180 and 360.
    Args:
        x (str): longitude
    Raises:
        Error: when latitude is <-180 and >360.
    Returns:
        x(float): converted longitude between 0 and 360
    """
    x = float(x)
    if (-180 < x) and (x < 0):
        print("lon is :", x)
        x = x % 360
        print("after modulo lon is :", x)
    if (x < 0) or (x > 360):
        raise argparse.ArgumentTypeError(
            "ERROR: Latitude of single point should be between 0 and 360 or -180 and 180."
        )
    return x


def get_git_sha():
    """
    Returns Git short SHA for the currect directory.
    """
    try:

        #os.abspath(__file__)
        sha = (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
            .strip()
            .decode()
        )
    except subprocess.CalledProcessError:
        sha = "NOT-A-GIT-REPOSITORY"
    return sha


def setup_logging(log_file, log_level):
    """
    Setup logging to log to console and log file.
    """

    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # setup log file
    one_mb = 1000000
    handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=one_mb, backupCount=10
    )

    fmt = logging.Formatter(
        "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%y-%m-%d %H:%M:%S",
    )

    handler.setFormatter(fmt)
    root_logger.addHandler(handler)

    # setup logging to console
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    root_logger.addHandler(stream_handler)

    # redirect stdout/err to log file
    StreamToLogger.setup_stdout()
    StreamToLogger.setup_stderr()


class StreamToLogger(object):
    """
    Custom class to log all stdout and stderr streams.
    modified from:
    https://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/
    """

    def __init__(
        self, stream, logger, log_level=logging.INFO, also_log_to_stream=False
    ):
        self.logger = logger
        self.stream = stream
        self.log_level = log_level
        self.linebuf = ""
        self.also_log_to_stream = also_log_to_stream

    @classmethod
    def setup_stdout(cls, also_log_to_stream=True):
        """
        Setup logger for stdout
        """
        stdout_logger = logging.getLogger("STDOUT")
        sl = StreamToLogger(sys.stdout, stdout_logger,
                            logging.INFO, also_log_to_stream)
        sys.stdout = sl

    @classmethod
    def setup_stderr(cls, also_log_to_stream=True):
        """
        Setup logger for stdout
        """
        stderr_logger = logging.getLogger("STDERR")
        sl = StreamToLogger(
            sys.stderr, stderr_logger, logging.ERROR, also_log_to_stream
        )
        sys.stderr = sl

    def write(self, buf):
        temp_linebuf = self.linebuf + buf
        self.linebuf = ""
        for line in temp_linebuf.splitlines(True):
            if line[-1] == "\n":
                self.logger.log(self.log_level, line.rstrip())
            else:
                self.linebuf += line

    def flush(self):
        if self.linebuf != "":
            self.logger.log(self.log_level, self.linebuf.rstrip())
        self.linebuf = ""


def main():

    defaults = configparser.ConfigParser()
    defaults.read(DEFAULTS_FILE)

    # Parse defaults
    dir_inputdata = defaults.get('main', 'clmforcingindir')
    dir_input_datm = defaults.get('datm_gswp3', 'dir')
    domain_file = defaults.get('domain', 'file')
    fdomain_in = os.path.join(dir_inputdata, domain_file)
    fdatmdomain_in = os.path.join(defaults.get('datm_gswp3', 'dir'), defaults.get('datm_gswp3', 'domain'))
    fsurfdat_78pft = os.path.join(dir_inputdata, defaults.get('surfdat', 'dir'),
                                  defaults.get('surfdat', 'surfdat_78pft'))
    fsurfdat_16pft = os.path.join(dir_inputdata, defaults.get('surfdat', 'dir'),
                                  defaults.get('surfdat', 'surfdat_16pft'))
    landuse_78pft = os.path.join(dir_inputdata, defaults.get('landuse', 'dir'),
                                  defaults.get('landuse', 'landuse_78pft'))
    landuse_16pft = os.path.join(dir_inputdata, defaults.get('landuse', 'dir'),
                                  defaults.get('landuse', 'landuse_16pft'))

    datm_solardir = defaults.get('datm_gswp3', 'solardir')
    datm_precdir = defaults.get('datm_gswp3', 'precdir')
    datm_tpqwdir = defaults.get('datm_gswp3', 'tpqwdir')
    datm_solartag = defaults.get('datm_gswp3', 'solartag')
    datm_prectag = defaults.get('datm_gswp3', 'prectag')
    datm_tpqwtag = defaults.get('datm_gswp3', 'tpqwtag')
    datm_solarname = defaults.get('datm_gswp3', 'solarname')
    datm_precname = defaults.get('datm_gswp3', 'precname')
    datm_tpqwname = defaults.get('datm_gswp3', 'tpqwname')

    args = get_parser().parse_args()

    # --------------------------------- #

    today = date.today()
    today_string = today.strftime("%Y%m%d")

    pwd = os.getcwd()

    log_file = os.path.join(pwd, today_string + ".log")

    log_level = logging.DEBUG
    setup_logging(log_file, log_level)
    log = logging.getLogger(__name__)

    print("User = " + myname)
    print("Current directory = " + pwd)

    # --------------------------------- #

    # print help and exit when no option is chosen
    if (args.run_type != "point" and args.run_type != "reg"):
        get_parser().print_help()
        quit()

    # -- Specify which types of data to subset
    create_domain = args.create_domain
    create_surfdata = args.create_surfdata
    create_landuse = args.create_landuse
    create_datm = args.create_datm

    crop_flag = args.crop_flag
    if crop_flag:
        num_pft = "78"
        fsurf_in = fsurfdat_78pft
        fluse_in = landuse_78pft
    else:
        num_pft = "16"
        fsurf_in = fsurfdat_16pft
        fluse_in = landuse_16pft

    print("crop_flag = " + crop_flag.__str__() + " => num_pft =" + num_pft)

    # -- Start and ending years for DATM data
    datm_syr = args.datm_syr
    datm_eyr = args.datm_eyr

    # --  Modify landunit structure
    overwrite_single_pft = args.overwrite_single_pft
    dominant_pft = args.dom_pft
    zero_nonveg_landunits = args.zero_nonveg
    uniform_snowpack = args.uni_snow
    no_saturation_excess = args.no_saturation_excess

    # --  Set input and output filenames
    # --  Specify input and output directories
    dir_output = args.out_dir
    if not os.path.isdir(dir_output):
        os.mkdir(dir_output)

    dir_output_datm = os.path.join(dir_output, "datmdata/")
    if create_datm:
        if not os.path.isdir(dir_output_datm):
            os.mkdir(dir_output_datm)
        print("dir_input_datm  : ", dir_input_datm)
        print("dir_output_datm : ", dir_output_datm)

    if args.create_user_mods:
        if not os.path.isdir(args.user_mods_dir):
            os.mkdir(args.user_mods_dir)

        cesmroot = path_to_ctsm_root()

        # -- Create user_nl_clm file
        nl_clm_base = os.path.join(cesmroot, "cime_config/user_nl_clm")
        nl_clm = os.path.join(args.user_mods_dir, "user_nl_clm")
        with open(nl_clm_base, 'r') as basefile, open(nl_clm, 'w') as userfile:
            for line in basefile:
                userfile.write(line)

        # -- Create user_nl_datm_streams file
        if create_datm:
            nl_datm_base = os.path.join(cesmroot, "components/cdeps/datm/cime_config/user_nl_datm_streams")
            nl_datm = os.path.join(args.user_mods_dir, "user_nl_datm_streams")
            with open(nl_datm_base, 'r') as basefile, open(nl_datm, 'w') as userfile:
                for line in basefile:
                    userfile.write(line)



    if args.run_type == "point":
        print(
            "----------------------------------------------------------------------------"
        )
        print(
            "This script extracts a single point from the global CTSM datasets."
        )

        # --  Specify point to extract
        plon = args.plon
        plat = args.plat
        site_name = args.site_name

        # --  Create SinglePoint Object
        single_point = SinglePointCase(
            plat,
            plon,
            site_name,
            create_domain,
            create_surfdata,
            create_landuse,
            create_datm,
            overwrite_single_pft,
            dominant_pft,
            zero_nonveg_landunits,
            uniform_snowpack,
            no_saturation_excess,
        )
        single_point.create_tag()

        # --  Create CTSM domain file
        if create_domain:
            # --  Specify land domain file  ---------------------------------
            fdomain_out = dir_output + single_point.add_tag_to_filename(
                fdomain_in, single_point.tag
            )
            print(fdomain_out)
            single_point.fdomain_in = fdomain_in
            single_point.fdomain_out = fdomain_out
            print("fdomain_in  :", fdomain_in)  #
            print("fdomain_out :", fdomain_out)  #
            single_point.create_domain_at_point()

        # --  Create CTSM surface data file
        if create_surfdata:
            # --  Specify surface file  ---------------------------------
            fsurf_out = dir_output + single_point.create_fileout_name(
                fsurf_in, single_point.tag
            )
            single_point.fsurf_in = fsurf_in
            single_point.fsurf_out = fsurf_out
            print("fsurf_in   :", fsurf_in)
            print("fsurf_out  :", fsurf_out)
            single_point.create_surfdata_at_point()
            if args.create_user_mods:
                nl_clm = open(os.path.join(args.user_mods_dir, "user_nl_clm"), "a")
                line = 'fsurdat = ' + "'" + fsurf_out + "'"
                nl_clm.write('\n' + line + '\n')
                nl_clm.close()

        # --  Create CTSM transient landuse data file
        if create_landuse:
            # --  Specify surface file  ---------------------------------
            fluse_out = dir_output + single_point.create_fileout_name(
                fluse_in, single_point.tag
            )
            single_point.fluse_in = fluse_in
            single_point.fluse_out = fluse_out
            print("fluse_in   :", fluse_in)
            print("fluse_out  :", fluse_out)
            single_point.create_landuse_at_point()
            if args.create_user_mods:
                nl_clm = open(os.path.join(args.user_mods_dir, "user_nl_clm"), "a")
                line = 'landuse = ' + fluse_out
                nl_clm.write('\n' + line + '\n')
                nl_clm.close()

        # --  Create single point atmospheric forcing data
        if create_datm:
            # --  Specify datm domain file  ---------------------------------
            fdatmdomain_out = dir_output_datm + single_point.add_tag_to_filename(
                fdatmdomain_in, single_point.tag
            )
            single_point.fdatmdomain_in = fdatmdomain_in
            single_point.fdatmdomain_out = fdatmdomain_out
            print("fdatmdomain_in   : ", fdatmdomain_in)
            print("fdatmdomain out  : ", fdatmdomain_out)
            single_point.create_datmdomain_at_point()
            single_point.datm_syr = datm_syr
            single_point.datm_eyr = datm_eyr
            single_point.dir_input_datm = dir_input_datm
            single_point.dir_output_datm = dir_output_datm
            single_point.dir_solar = datm_solardir
            single_point.dir_prec = datm_precdir
            single_point.dir_tpqw = datm_tpqwdir
            single_point.tag_solar = datm_solartag
            single_point.tag_prec = datm_prectag
            single_point.tag_tpqw = datm_tpqwtag
            single_point.name_solar = datm_solarname
            single_point.name_prec = datm_precname
            single_point.name_tpqw = datm_tpqwname
            single_point.create_datm_at_point(args.create_user_mods, nl_datm)

        if args.create_user_mods:
            shell_commands_file = open(os.path.join(args.user_mods_dir,
            "shell_commands"), 'w')
            single_point.write_shell_commands(shell_commands_file)


        print("Successfully ran script for single point.")
        exit()

    elif args.run_type == "reg":
        print("Running the script for the region")

        # --  Specify region to extract
        lat1 = args.lat1
        lat2 = args.lat2

        lon1 = args.lon1
        lon2 = args.lon2

        reg_name = args.reg_name

        region = RegionalCase(
            lat1,
            lat2,
            lon1,
            lon2,
            reg_name,
            create_domain,
            create_surfdata,
            create_landuse,
            create_datm,
        )

        print(region)

        region.create_tag()

        # --  Set time stamp
        command = 'date "+%y%m%d"'
        x2 = subprocess.Popen(command, stdout=subprocess.PIPE, shell="True")
        x = x2.communicate()
        timetag = x[0].strip()
        print(timetag)

        # --  Specify land domain file  ---------------------------------
        fdomain_out = (
            dir_output + "domain.lnd.fv1.9x2.5_gx1v7." + region.tag + "_170518.nc"
        )
        # SinglePointCase.set_fdomain (fdomain)
        region.fdomain_in = fdomain_in
        region.fdomain_out = fdomain_out

        # --  Specify surface data file  --------------------------------
        fsurf_out = (
            dir_output
            + "surfdata_1.9x2.5_78pfts_CMIP6_simyr1850_"
            + region.tag
            + "_c170824.nc"
        )
        region.fsurf_in = fsurf_in
        region.fsurf_out = fsurf_out

        # --  Specify landuse file  -------------------------------------
        fluse_out = (
            dir_output
            + "landuse.timeseries_1.9x2.5_hist_78pfts_CMIP6_simyr1850-2015_"
            + region.tag
            + ".c170824.nc"
        )
        region.fluse_in = fluse_in
        region.fluse_out = fluse_out

        # --  Create CTSM domain file
        if create_domain:
            region.create_domain_at_reg()

        # --  Create CTSM surface data file
        if create_surfdata:
            region.create_surfdata_at_reg()

        # --  Create CTSM transient landuse data file
        if create_landuse:
            region.create_landuse_at_reg()
        print("Successfully ran script for a regional case.")
        exit()
