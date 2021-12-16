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
To see the available options for single point or regional cases:
    ./subset_data.py --help
-------------------------------------------------------------------
This script extracts domain files, surface dataset, and DATM files
at either a single point or a region using a global dataset. Currently this
script subsets default surface, landuse, and DATM files, which can be seen in
the defaults.cfg file.

To run a single-point or regional case using this data with the NUOPC driver,
you must update the variable(s) `fsurdat` and/or `landuse` in the user_nl_clm namelist
file to be the full path to the subset files. This script will automatically create this
file using the flag --create-user-mods.
To use subset climate data, the namelist file user_nl_datm_streams must also
be updated - this script will automatically create this file with
--create-user-mods. This flag will also create necessary single-point xml
commands in the file shell_commands.

To use the created user mods with a case use --user-mods-dir PATH/TO/USER/MODS
in the ./create.newcase call.

By default, this script only extracts surface dataset. For extracting other
files, the appropriate flags should be used.

To run this script the following packages are required:
        - numpy
        - xarray

-------------------------------------------------------------------
To run the script for a single point:
    ./subset_data.py point

To run the script for a region:
    ./subset_data.py region

To remove NPL from your environment on Cheyenne/Casper:
    deactivate
-------------------------------------------------------------------
"""

# TODO [NS]:
# -[] Automatic downloading of missing files if they are missing

# -- Import libraries

# -- standard libraries
import os
import sys
import logging
import argparse
import textwrap
import configparser

from getpass import getuser
from argparse import ArgumentParser

# -- import local classes for this script
from ctsm.site_and_regional.base_case import DatmFiles
from ctsm.site_and_regional.single_point_case import SinglePointCase
from ctsm.site_and_regional.regional_case import RegionalCase
from ctsm.path_utils import path_to_ctsm_root

# -- import ctsm logging flags
from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)

_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", 'python'))
sys.path.insert(1, _CTSM_PYTHON)

DEFAULTS_FILE = "default_data.cfg"

logger = logging.getLogger(__name__)


def get_parser():
    """
    Get the parser object for subset_data.py script.

    Returns:
        parser (ArgumentParser): ArgumentParser which includes all the parser information.

    """
    parser = ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help
    subparsers = parser.add_subparsers(
        help="Two possible ways to run this script, either:", dest="run_type"
    )
    pt_parser = subparsers.add_parser("point", help="Run script for a single point.")
    rg_parser = subparsers.add_parser("region", help="Run script for a region.")

    # -- single point parser options
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
        "--variable-snow-fraction",
        help="Creating datasets using variable snow fraction. If unset, snow is set to  uniform "
             "fraction.",
        action="store_false",
        dest="uni_snow",
        required=False,
        default=True,
    )
    pt_parser.add_argument(
        "--allow-multiple-pft",
        help="Creating dataset using multiple pft. If unset, it assumes the whole grid is 100% "
             "single PFT set by --dom-pft.",
        action="store_false",
        dest="overwrite_single_pft",
        required=False,
        default=True,
    )
    pt_parser.add_argument(
        "--zero-nonveg",
        help="Creating dataset by setting all non-vegetation landunits to zero.",
        action="store_true",
        dest="zero_nonveg",
        required=False,
        default=False,
    )
    pt_parser.add_argument(
        "--allow-saturation-excess",
        help="Creating dataset allowing saturatated conditions. If unset saturation_excess is set "
             "to zero.",
        action="store_true",
        dest="saturation_excess",
        required=False,
        default=False,
    )
    pt_parser.add_argument(
        "--dompft",
        help="Dominant PFT type if we set the grid to 100% one PFT [default: %(default)s].",
        action="store",
        dest="dom_pft",
        type=int,
        default=7,
    )
    pt_parser.add_argument(
        "--datm-from-tower",
        help="Create DATM forcing data at single point for a tower data.",
        action="store_true",
        dest="datm_tower",
        required=False,
        default=False,
    )
    # -- region-specific parser options
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
        "--create-mesh",
        help="Flag for subsetting mesh file.",
        action="store_true",
        dest="create_mesh",
        required=False,
        default=False,
    )

    # -- common options between both subparsers
    for subparser in [pt_parser, rg_parser]:
        subparser.add_argument(
            "--create-domain",
            help="Create CLM domain file at single point/region.",
            action="store_true",
            dest="create_domain",
            required=False,
            default=False,
        )
        subparser.add_argument(
            "--create-surface",
            help="Create surface data file at single point/region.",
            action="store_true",
            dest="create_surfdata",
            required=False,
            default=False,
        )
        subparser.add_argument(
            "--create-landuse",
            help="Create landuse data file at single point/region.",
            action="store_true",
            dest="create_landuse",
            required=False,
            default=False,
        )
        subparser.add_argument(
            "--create-datm",
            help="Create DATM forcing data at single point/region.",
            action="store_true",
            dest="create_datm",
            required=False,
            default=False,
        )
        subparser.add_argument(
            "--create-user-mods",
            help="Create user mods directories and files for running CTSM with the subset data.",
            action="store_true",
            dest="create_user_mods",
            required=False,
            default=True,
        )
        subparser.add_argument(
            "--datm-syr",
            help="Start year for creating DATM forcing at single point/region. [default: %("
                 "default)s]",
            action="store",
            dest="datm_syr",
            required=False,
            type=int,
            default=1901,
        )
        subparser.add_argument(
            "--datm-eyr",
            help="End year for creating DATM forcing at single point/region. "
                 "[default: %(default)s]",
            action="store",
            dest="datm_eyr",
            required=False,
            type=int,
            default=2014,
        )
        subparser.add_argument(
            "--crop",
            help="Create datasets using the extensive list of prognostic crop types.",
            action="store_true",
            dest="crop_flag",
            required=False,
            default=False,
        )

        if subparser == pt_parser:
            parser_name = "single_point"
        else:
            parser_name = "regional"

        subparser.add_argument(
            "--outdir",
            help="Output directory. \n [default: %(default)s]",
            action="store",
            dest="out_dir",
            type=str,
            default=os.path.join(os.getcwd(), "subset_data_" + parser_name),
        )

    # -- print help for both subparsers
    parser.epilog = textwrap.dedent(
        f"""\
         {pt_parser.format_help()}
         {rg_parser.format_help()}
         """
    )
    return parser


def plat_type(plat):
    """
    Function to define lat type for the parser
    and
    raise error if latitude is not between -90 and 90.
    Args:
        plat(str): latitude
    Raises:
        Error (ArgumentTypeError): when plat (latitude) is not between -90 and 90.
    Returns:
        plat (float): latitude in float
    """
    plat = float(plat)
    if (plat < -90) or (plat > 90):
        raise argparse.ArgumentTypeError("ERROR: Latitude should be between -90 and 90.")
    return plat


def plon_type(plon):
    """
    Function to define lon type for the parser and
    convert negative longitudes and
    raise error if lon is not between -180 and 360.
    Args:
        plon (str): longitude
    Raises:
        Error (ArgumentTypeError): when longitude is <-180 and >360.
    Returns:
        plon(float): converted longitude between 0 and 360
    """
    plon = float(plon)
    if -180 <= plon < 0:
        logger.info("lon is: %f", plon)
        plon = plon % 360
        logger.info("after modulo lon is: %f", plon)
    if plon < 0 or plon > 360:
        raise argparse.ArgumentTypeError("ERROR: Longitude of single point should be between 0 and "
                                         "360 or -180 and 180.")
    return plon


def setup_user_mods(user_mods_dir, cesmroot):
    """
    Sets up the user mods files and directories
    """
    if not os.path.isdir(user_mods_dir):
        os.mkdir(user_mods_dir)

    nl_clm_base = os.path.join(cesmroot, "cime_config/user_nl_clm")
    nl_clm = os.path.join(user_mods_dir, "user_nl_clm")
    with open(nl_clm_base, "r") as basefile, open(nl_clm, "w") as user_file:
        for line in basefile:
            user_file.write(line)

    nl_datm_base = os.path.join(cesmroot, "components/cdeps/datm/cime_config"
                                          "/user_nl_datm_streams")
    nl_datm = os.path.join(user_mods_dir, "user_nl_datm_streams")
    with open(nl_datm_base, "r") as base_file, open(nl_datm, 'w') as user_file:
        for line in base_file:
            user_file.write(line)

def setup_files(args, defaults, cesmroot):
    """
    Sets up the files and folders needed for this program
    """

    if args.user_mods_dir == "":
        args.user_mods_dir = os.path.join(args.out_dir, "user_mods")
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    if args.create_user_mods:
        setup_user_mods(args.user_mods_dir, cesmroot)

    # DATM data
    datm_type = 'datm_gswp3'
    dir_output_datm = "datmdata"
    dir_input_datm = defaults.get(datm_type, "dir")
    if args.create_datm:
        if not os.path.isdir(os.path.join(args.out_dir, dir_output_datm)):
            os.mkdir(os.path.join(args.out_dir, dir_output_datm))
        logger.info("dir_input_datm : %s", dir_input_datm)
        logger.info("dir_output_datm: %s", os.path.join(args.out_dir, dir_output_datm))

    # if the crop flag is on - we need to use a different land use and surface data file
    if args.crop_flag:
        num_pft = "78"
        fsurf_in = defaults.get("surfdat", "surfdat_78pft")
        fluse_in = defaults.get("landuse", "landuse_78pft")
    else:
        num_pft = "16"
        fsurf_in = defaults.get("surfdat", "surfdat_16pft")
        fluse_in = defaults.get("landuse", "landuse_16pft")
    logger.debug("crop_flag = %s => num_pft = %s", args.crop_flag.__str__(), num_pft)

    file_dict = {'main_dir': defaults.get("main", "clmforcingindir"),
                 'fdomain_in': defaults.get("domain", "file"),
                 'fsurf_dir': os.path.join(defaults.get("main", "clmforcingindir"),
                                           os.path.join(defaults.get("surfdat", "dir"))),
                 'fluse_dir': os.path.join(defaults.get("main", "clmforcingindir"),
                                           os.path.join(defaults.get("landuse", "dir"))),
                 'fsurf_in': fsurf_in,
                 'fluse_in': fluse_in,
                 'datm_tuple': DatmFiles(dir_input_datm,
                                         dir_output_datm,
                                         defaults.get(datm_type, "domain"),
                                         defaults.get(datm_type, 'solardir'),
                                         defaults.get(datm_type, 'precdir'),
                                         defaults.get(datm_type, 'tpqwdir'),
                                         defaults.get(datm_type, 'solartag'),
                                         defaults.get(datm_type, 'prectag'),
                                         defaults.get(datm_type, 'tpqwtag'),
                                         defaults.get(datm_type, 'solarname'),
                                         defaults.get(datm_type, 'precname'),
                                         defaults.get(datm_type, 'tpqwname'))
                 }

    return file_dict


def subset_point(args, file_dict: dict):
    """
    Subsets surface, domain, land use, and/or DATM files at a single point
    """

    logger.info("----------------------------------------------------------------------------")
    logger.info("This script extracts a single point from the global CTSM datasets.")

    # --  Create SinglePoint Object
    single_point = SinglePointCase(
        args.plat,
        args.plon,
        args.site_name,
        args.create_domain,
        args.create_surfdata,
        args.create_landuse,
        args.create_datm,
        args.create_user_mods,
        args.overwrite_single_pft,
        args.dom_pft,
        args.zero_nonveg,
        args.uni_snow,
        args.saturation_excess,
        args.out_dir,
    )

    single_point.create_tag()
    logger.debug(single_point)

    # --  Create CTSM domain file
    if single_point.create_domain:
        single_point.create_domain_at_point(file_dict["main_dir"], file_dict["fdomain_in"])

    # --  Create CTSM surface data file
    if single_point.create_surfdata:
        single_point.create_surfdata_at_point(file_dict["fsurf_dir"], file_dict["fsurf_in"],
                                              args.user_mods_dir)

    # --  Create CTSM transient landuse data file
    if single_point.create_landuse:
        single_point.create_landuse_at_point(file_dict["fluse_dir"], file_dict["fluse_in"],
                                             args.user_mods_dir)

    # --  Create single point atmospheric forcing data
    if single_point.create_datm:
        # subset DATM domain file
        single_point.create_datmdomain_at_point(file_dict["datm_tuple"])

        # subset the DATM data
        nl_datm = os.path.join(args.user_mods_dir, "user_nl_datm_streams")
        single_point.create_datm_at_point(file_dict['datm_tuple'], args.datm_syr, args.datm_eyr,
                                          nl_datm)

    # -- Write shell commands
    if single_point.create_user_mods:
        single_point.write_shell_commands(os.path.join(args.user_mods_dir, "shell_commands"))

    logger.info("Successfully ran script for single point.")


def subset_region(args, file_dict: dict):
    """
    Subsets surface, domain, land use, and/or DATM files for a region
    """

    logger.info("----------------------------------------------------------------------------")
    logger.info("This script extracts a region from the global CTSM datasets.")

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
        args.create_user_mods,
        args.out_dir,
    )

    region.create_tag()
    logger.debug(region)

    # --  Create CTSM domain file
    if region.create_domain:
        region.create_domain_at_reg(file_dict["main_dir"], file_dict["fdomain_in"])

    # --  Create CTSM surface data file
    if region.create_surfdata:
        region.create_surfdata_at_reg(file_dict["fsurf_dir"], file_dict["fsurf_in"],
                                      args.user_mods_dir)

    # --  Create CTSM transient landuse data file
    if region.create_landuse:
        region.create_landuse_at_reg(file_dict["fluse_dir"], file_dict["fluse_in"],
                                     args.user_mods_dir)

    logger.info("Successfully ran script for a regional case.")
    sys.exit()


def main():
    """
    Calls functions that subset surface, landuse, domain, and/or DATM files for a region or a
    single point.
    """

    # add logging flags from ctsm_logging
    setup_logging_pre_config()
    parser = get_parser()
    add_logging_args(parser)
    args = parser.parse_args()
    process_logging_args(args)

    # parse defaults file
    cesmroot = path_to_ctsm_root()
    defaults = configparser.ConfigParser()
    defaults.read(os.path.join(cesmroot, 'tools/site_and_regional', DEFAULTS_FILE))

    # --------------------------------- #

    myname = getuser()
    pwd = os.getcwd()
    logger.info("User = %s", myname)
    logger.info("Current directory = %s", pwd)

    # --------------------------------- #

    # print help and exit when no option is chosen
    if args.run_type != "point" and args.run_type != "region":
        get_parser().print_help()
        sys.exit()

    # create files and folders necessary and return dictionary of file/folder locations
    file_dict = setup_files(args, defaults, cesmroot)

    if args.run_type == "point":
        subset_point(args, file_dict)
    elif args.run_type == "reg":
        subset_region(args, file_dict)
