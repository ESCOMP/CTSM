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
from ctsm.site_and_regional.regional_case import RegionalCase
from ctsm.site_and_regional.single_point_case import SinglePointCase
from ctsm.site_and_regional.base_case import BaseCase
from ctsm.path_utils import path_to_ctsm_root
from ctsm import add_cime_to_path

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
_CTSM_PYTHON = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", 'python'))
sys.path.insert(1, _CTSM_PYTHON)

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
        default="",
    )
    pt_parser.add_argument(
        "--outdir",
        help="Output directory. [default: %(default)s]",
        action="store",
        dest="out_dir",
        type=str,
        default="/glade/scratch/" + myname + "/single_point/",
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
        default="",
    )
    rg_parser.add_argument(
        "--outdir",
        help="Output directory. [default: %(default)s]",
        action="store",
        dest="out_dir",
        type=str,
        default="/glade/scratch/" + myname + "/region/",
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

        # os.abspath(__file__)
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

    # parse arguments
    args = get_parser().parse_args()

    # parse defaults file
    defaults = configparser.ConfigParser()
    defaults.read(DEFAULTS_FILE)

    # --------------------------------- #

    today = date.today()
    today_string = today.strftime("%Y%m%d")

    pwd = os.getcwd()

    log_file = os.path.join(pwd, today_string + ".log")

    log_level = logging.DEBUG
    setup_logging(log_file, log_level)

    print("User = " + myname)
    print("Current directory = " + pwd)

    # --------------------------------- #

    # print help and exit when no option is chosen
    if (args.run_type != "point" and args.run_type != "reg"):
        get_parser().print_help()
        quit()

    # if the crop flag is on - we need to use a different landuse and surface
    # data file
    if args.crop_flag:
        num_pft = "78"
        fsurf_in = defaults.get('surfdat', 'surfdat_78pft')
        fluse_in = defaults.get('landuse', 'landuse_78pft')
    else:
        num_pft = "16"
        fsurf_in = defaults.get('surfdat', 'surfdat_16pft')
        fluse_in = defaults.get('landuse', 'landuse_16pft')

    print("crop_flag = " + args.crop_flag.__str__() + " => num_pft =" + num_pft)

    # --  Specify input and output directories and files

    # top-level output directory
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    # datm data
    dir_output_datm = "datmdata"
    dir_input_datm = defaults.get('datm_gswp3', 'dir')
    if args.create_datm:
        if not os.path.isdir(os.path.join(args.out_dir, dir_output_datm)):
            os.mkdir(os.path.join(args.out_dir, dir_output_datm))
        print("dir_input_datm  : ", dir_input_datm)
        print("dir_output_datm : ", os.path.join(args.out_dir, dir_output_datm))

    # -- Set up user mods directories and base files
    if args.create_user_mods:
        if not os.path.isdir(args.user_mods_dir):
            os.mkdir(args.user_mods_dir)

        cesmroot = path_to_ctsm_root()

        # -- Create empty user_nl_clm file
        if args.create_surfdata or args.create_landuse:
            nl_clm_base = os.path.join(cesmroot, "cime_config/user_nl_clm")
            nl_clm = os.path.join(args.user_mods_dir, "user_nl_clm")
            with open(nl_clm_base, 'r') as basefile, open(nl_clm, 'w') as userfile:
                for line in basefile:
                    userfile.write(line)

        # -- Create empty user_nl_datm_streams file
        if args.create_datm:
            nl_datm_base = os.path.join(
                cesmroot, "components/cdeps/datm/cime_config/user_nl_datm_streams")
            nl_datm = os.path.join(args.user_mods_dir, "user_nl_datm_streams")
            with open(nl_datm_base, 'r') as basefile, open(nl_datm, 'w') as userfile:
                for line in basefile:
                    userfile.write(line)

    # Default files
    dir_inputdata = defaults.get('main', 'clmforcingindir')
    dir_inputsurf = defaults.get('surfdat', 'dir')
    dir_inputluse = defaults.get('landuse', 'dir')
    fdomain_in = os.path.join(dir_inputdata, defaults.get('domain', 'file'))
    fdatmdomain_in = os.path.join(dir_input_datm, defaults.get('datm_gswp3', 'domain'))
    datm_solardir = defaults.get('datm_gswp3', 'solardir')
    datm_precdir = defaults.get('datm_gswp3', 'precdir')
    datm_tpqwdir = defaults.get('datm_gswp3', 'tpqwdir')
    datm_solartag = defaults.get('datm_gswp3', 'solartag')
    datm_prectag = defaults.get('datm_gswp3', 'prectag')
    datm_tpqwtag = defaults.get('datm_gswp3', 'tpqwtag')
    datm_solarname = defaults.get('datm_gswp3', 'solarname')
    datm_precname = defaults.get('datm_gswp3', 'precname')
    datm_tpqwname = defaults.get('datm_gswp3', 'tpqwname')

    if args.run_type == "point":
        print("----------------------------------------------------------------------------")
        print( "This script extracts a single point from the global CTSM datasets.")

        # --  Create SinglePoint Object
        single_point = SinglePointCase(
            args.plat,
            args.plon,
            args.site_name,
            args.create_domain,
            args.create_surfdata,
            args.create_landuse,
            args.create_datm,
            args.overwrite_single_pft,
            args.dom_pft,
            args.zero_nonveg,
            args.uni_snow,
            args.no_saturation_excess,
            args.out_dir
        )
        single_point.create_tag()

        # --  Create CTSM domain file
        if single_point.create_domain:
            # --  Specify land domain file  ---------------------------------
            single_point.fdomain_in = os.path.join(dir_inputdata, fdomain_in)
            single_point.fdomain_out = single_point.add_tag_to_filename(fdomain_in, single_point.tag)
            print("fdomain_in  :", single_point.fdomain_in)
            print("fdomain_out :", os.path.join(single_point.output_dir, single_point.fdomain_out))
            single_point.create_domain_at_point()

        # --  Create CTSM surface data file
        if single_point.create_surfdata:
            # --  Specify surface file  ---------------------------------
            single_point.fsurf_in = os.path.join(dir_inputdata, dir_inputsurf, fsurf_in)
            single_point.fsurf_out = single_point.create_fileout_name(fsurf_in, single_point.tag)
            print("fsurf_in   :", single_point.fsurf_in)
            print("fsurf_out  :", single_point.fsurf_out)
            single_point.create_surfdata_at_point()

            # write to user_nl_clm if specified
            if args.create_user_mods:
                nl_clm = open(os.path.join(args.user_mods_dir, "user_nl_clm"), "a")
                line = 'fsurdat = ' + "'$CLM_USRDAT_DIR/" + single_point.fsurf_out + "'"
                nl_clm.write('\n' + line + '\n')
                nl_clm.close()

        # --  Create CTSM transient landuse data file
        if single_point.create_landuse:
            # --  Specify surface file  ---------------------------------
            single_point.fluse_in = os.path.join(dir_inputdata, dir_inputluse, fluse_in)
            single_point.fluse_out = single_point.create_fileout_name(fluse_in, single_point.tag)
            print("fluse_in   :", single_point.fluse_in)
            print("fluse_out  :", single_point.fluse_out)
            single_point.create_landuse_at_point()

            # write to user_nl_clm data if specified
            if args.create_user_mods:
                nl_clm = open(os.path.join(args.user_mods_dir, "user_nl_clm"), "a")
                line = 'fsurdat = ' + "'$CLM_USRDAT_DIR/" + single_point.fluse_out + "'"
                nl_clm.write('\n' + line + '\n')
                nl_clm.close()

        # --  Create single point atmospheric forcing data
        if single_point.create_datm:
            # --  Specify datm and subset domain file  ---------------------------------
            single_point.fdatmdomain_in = os.path.join(dir_input_datm, fdatmdomain_in)
            single_point.fdatmdomain_out = os.path.join(dir_output_datm,
                                                        single_point.add_tag_to_filename(single_point.fdatmdomain_in,
                                                                                         single_point.tag))
            print("fdatmdomain_in   : ", single_point.fdatmdomain_in)
            print("fdatmdomain out  : ", single_point.fdatmdomain_out)
            single_point.create_datmdomain_at_point()

            # -- Specify DATM directories, tags, and stream names
            single_point.datm_syr = args.datm_syr
            single_point.datm_eyr = args.datm_eyr
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

        # -- Write shell commands
        if args.create_user_mods:
            shell_commands_file = open(os.path.join(args.user_mods_dir,"shell_commands"), 'w')
            single_point.write_shell_commands(shell_commands_file)

        print("Successfully ran script for single point.")
        exit()

    elif args.run_type == "reg":
        print("----------------------------------------------------------------------------")
        print( "This script extracts a region from the global CTSM datasets.")

        # --  Create Region Object
        region = RegionalCase(
            args.lat1,
            args.lat2,
            args.lon1,
            args.lon2,
            args.reg_name,
            args.create_domain,
            args.create_surfdata,
            args.create_landuse,
            args.create_datm,
        )
        region.create_tag()

        # --  Set time stamp
        command = 'date "+%y%m%d"'
        x2 = subprocess.Popen(command, stdout=subprocess.PIPE, shell="True")
        x = x2.communicate()
        timetag = x[0].strip()
        print(timetag)

        # --  Specify land domain file  ---------------------------------
        fdomain_out = os.path.join(args.out_dir, "domain.lnd.fv1.9x2.5_gx1v7." + region.tag + "_170518.nc")

        # SinglePointCase.set_fdomain (fdomain)
        region.fdomain_in = fdomain_in
        region.fdomain_out = fdomain_out

        # --  Specify surface data file  --------------------------------
        fsurf_out = os.path.join(
                args.out_dir,
                "surfdata_1.9x2.5_78pfts_CMIP6_simyr1850_"
                + region.tag
                + "_c170824.nc"
        )
        region.fsurf_in = fsurf_in
        region.fsurf_out = fsurf_out

        # --  Specify landuse file  -------------------------------------
        fluse_out = os.path.join(
                args.out_dir,
                "landuse.timeseries_1.9x2.5_hist_78pfts_CMIP6_simyr1850-2015_"
                + region.tag
                + ".c170824.nc"
        )
        region.fluse_in = fluse_in
        region.fluse_out = fluse_out

        # --  Create CTSM domain file
        if region.create_domain:
            region.create_domain_at_reg()

        # --  Create CTSM surface data file
        if region.create_surfdata:
            region.create_surfdata_at_reg()

        # --  Create CTSM transient landuse data file
        if region.create_landuse:
            region.create_landuse_at_reg()
        print("Successfully ran script for a regional case.")
        exit()
