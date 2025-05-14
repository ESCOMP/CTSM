"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
Instructions for running using conda python environments:

../../py_env_create
conda activate ctsm_py
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
import logging
import textwrap
import configparser

from getpass import getuser
import argparse
from argparse import ArgumentParser

# -- import local classes for this script
from ctsm.site_and_regional.base_case import DatmFiles
from ctsm.site_and_regional.single_point_case import SinglePointCase
from ctsm.site_and_regional.regional_case import RegionalCase
from ctsm.args_utils import plat_type, plon_type
from ctsm.path_utils import path_to_ctsm_root
from ctsm.utils import abort
from ctsm.config_utils import check_lon1_lt_lon2
from ctsm.longitude import Longitude

# -- import ctsm logging flags
from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)

DEFAULTS_CONFIG = "tools/site_and_regional/default_data_2000.cfg"

logger = logging.getLogger(__name__)


def _add_lon_type_arg(this_parser):
    lon_type_help_str = (
        "Whether longitudes are in the [-180, 180] format (centered around the Prime/0th"
        " Meridian) or the [0, 360] format (centered around the 180th Meridian)."
        " Choose by specifying the upper limit."
    )
    this_parser.add_argument(
        "--lon-type",
        help=lon_type_help_str,
        required=False,
        default=None,
        type=int,
        choices=[180, 360],
    )
    return this_parser


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
        default=287.8,  # Must be unambiguous: Either < 0 or > 180
    )
    pt_parser = _add_lon_type_arg(pt_parser)
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
        "--uniform-snowpack",
        help="Modify surface data to have a uniform snow fraction.",
        action="store_true",
        dest="uni_snow",
        required=False,
    )
    pt_parser.add_argument(
        "--include-nonveg",
        help="Do not zero non-vegetation land units in the surface data.",
        action="store_true",
        dest="include_nonveg",
        required=False,
    )
    pt_parser.add_argument(
        "--cap-saturation",
        help="Modify surface data to not allow saturation excess.",
        action="store_true",
        dest="cap_saturation",
        required=False,
    )
    pt_parser.add_argument(
        "--evenly_split_cropland",
        help="Introduce equal areas of all crops",
        action="store_true",
        dest="evenly_split_cropland",
        required=False,
    )
    pt_parser.add_argument(
        "--dompft",
        help="Dominant PFT(s): if we set the grid to 100%% one or multiple PFTs \
        [default: %(default)s].",
        action="store",
        dest="dom_pft",
        type=int,
        default=None,
        nargs="*",
    )
    pt_parser.add_argument(
        "--pctpft",
        help="Percetages of each pft (set by --dompft) on the land unit.",
        action="store",
        dest="pct_pft",
        type=float,
        default=None,
        nargs="*",
    )
    pt_parser.add_argument(
        "--cth",
        help="canopy top height for pft",
        action="store",
        dest="cth",
        type=float,
        default=None,
        nargs="*",
    )
    pt_parser.add_argument(
        "--cbh",
        help="canopy bottom height for pft",
        action="store",
        dest="cbh",
        type=float,
        default=None,
        nargs="*",
    )

    # -- region-specific parser options
    rg_parser.add_argument(
        "--lat1",
        help="Region southernmost latitude. [default: %(default)s]",
        action="store",
        dest="lat1",
        required=False,
        type=plat_type,
        default=-40,
    )
    rg_parser.add_argument(
        "--lat2",
        help="Region northernmost latitude. [default: %(default)s]",
        action="store",
        dest="lat2",
        required=False,
        type=plat_type,
        default=15,
    )
    rg_parser.add_argument(
        "--lon1",
        help=("Region westernmost longitude. [default: %(default)s]"),
        action="store",
        dest="lon1",
        required=False,
        type=plon_type,
        default=275.0,  # Must be unambiguous: Either < 0 or > 180
    )
    rg_parser.add_argument(
        "--lon2",
        help=("Region easternmost longitude. [default: %(default)s]"),
        action="store",
        dest="lon2",
        required=False,
        type=plon_type,
        default=330.0,  # Must be unambiguous: Either < 0 or > 180
    )
    rg_parser = _add_lon_type_arg(rg_parser)
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
        help="Subset a mesh file for a region.",
        action="store_true",
        dest="create_mesh",
        required=False,
    )

    # -- common options between both subparsers
    for subparser in [pt_parser, rg_parser]:
        subparser.add_argument(
            "--create-surface",
            help="Create surface data file at single point/region.",
            action="store_true",
            dest="create_surfdata",
            required=False,
        )
        subparser.add_argument(
            "--surf-year",
            help="Year for surface data file at single point/region \
            (and start year for land-use timeseries).",
            action="store",
            dest="surf_year",
            type=int,
            default=2000,
            required=False,
        )
        subparser.add_argument(
            "--create-landuse",
            help="Create landuse data file at a single point/region.",
            action="store_true",
            dest="create_landuse",
            required=False,
        )
        subparser.add_argument(
            "--create-datm",
            help="Create DATM forcing data at a single point/region.",
            action="store_true",
            dest="create_datm",
            required=False,
        )
        subparser.add_argument(
            "--create-domain",
            help="Create CLM domain file for a single point/region \
            Domain files are not needed for NUOPC cases, \
            but are needed to create mesh files that are needed for NUOPC cases.",
            action="store_true",
            dest="create_domain",
            required=False,
        )
        subparser.add_argument(
            "--create-user-mods",
            help="Create user mods directories and files for running CTSM with the subset data.",
            action="store_true",
            dest="create_user_mods",
            required=False,
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
        subparser.add_argument(
            "--user-mods-dir",
            help="User mods directory.",
            action="store",
            dest="user_mods_dir",
            type=str,
            default="",
        )

        subparser.add_argument(
            "--out-surface",
            help="Output surface dataset name \
            (if you want to override the default based on the current date). \n \
            (only valid if outputing a surface dataset)",
            action="store",
            dest="out_surface",
            type=str,
        )
        cesmroot = path_to_ctsm_root()
        defaults_file = os.path.join(cesmroot, DEFAULTS_CONFIG)
        subparser.add_argument(
            "--cfg-file",
            help="Default configure file to use for default filenames.",
            action="store",
            dest="config_file",
            type=str,
            default=defaults_file,
        )
        subparser.add_argument(
            "--overwrite",
            help="Flag to overwrite if the files already exists.",
            action="store_true",
            dest="overwrite",
        )
        subparser.add_argument(
            "--inputdata-dir",
            help="Top level path to the CESM inputdata directory.",
            action="store",
            dest="inputdatadir",
            type=str,
            default="defaults.cfg",
        )
        add_logging_args(subparser)

    # -- print help for both subparsers
    parser.epilog = textwrap.dedent(
        f"""\
         {pt_parser.format_help()}
         {rg_parser.format_help()}
         """
    )
    return parser


def check_surf_year(args):
    """
    Check command-line arguments w/r/t --surf-year
    """
    if args.surf_year != 2000 and not args.create_surfdata:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n --surf-year option is set to something besides the default of 2000
                \n without the --create-surface option"
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    if args.surf_year != 1850 and args.create_landuse:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n --surf-year option is NOT set to 1850 and the --create-landuse option
                \n is selected which requires it to be 1850 (see
                https://github.com/ESCOMP/CTSM/issues/2018)
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    if args.surf_year not in [1850, 2000]:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n --surf-year option can only be set to 1850 or 2000
                """
        )
        raise argparse.ArgumentError(None, err_msg)


def check_args(args):
    """Check the command line arguments"""
    # --------------------------------- #
    # print help and exit when no option is chosen
    if args.run_type not in ("point", "region"):
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n Must supply a positional argument: 'point' or 'region'.
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    args = process_args(args)

    if not any([args.create_surfdata, args.create_landuse, args.create_datm, args.create_domain]):
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n Must supply one of:
                \n --create-surface \n --create-landuse \n --create-datm \n --create-domain \n \n
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    if not os.path.exists(args.config_file):
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n Entered default config file does not exist"
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    if args.out_surface and not args.create_surfdata:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n out-surface option is given without the --create-surface option"
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    if args.create_landuse and not args.create_surfdata:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n --create-landuse option requires the --create-surface option:
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    # Checks related to --surf-year
    check_surf_year(args)

    if args.out_surface and os.path.exists(args.out_surface) and not args.overwrite:
        err_msg = textwrap.dedent(
            """\
                \n ------------------------------------
                \n out-surface filename exists and the overwrite option was not also selected"
                """
        )
        raise argparse.ArgumentError(None, err_msg)

    if args.run_type == "region" and args.create_user_mods:
        if not args.create_mesh:
            err_msg = textwrap.dedent(
                """\
                      \n ------------------------------------
                      \nERROR: For regional cases, you can not create user_mods
                      \nwithout creating the mesh file.

                      \nPlease rerun the script adding --create-mesh to subset the mesh file."
                      """
            )
            raise argparse.ArgumentError(None, err_msg)

    if args.run_type == "region" and args.create_mesh:
        if not args.create_domain:
            err_msg = textwrap.dedent(
                """\
                      \n ------------------------------------
                      \nERROR: For regional cases, you can not create mesh files
                      \nwithout creating the domain file.

                      \nPlease rerun the script adding --create-domain to subset the domain file."
                      """
            )
            raise argparse.ArgumentError(None, err_msg)

    if args.run_type == "region" and args.create_datm:
        err_msg = textwrap.dedent(
            """\
                    \n ------------------------------------
                    \nERROR: For regional cases, you can not subset datm data
                    \n (see https://github.com/ESCOMP/CTSM/issues/2110)
                    \n but you can just use the global data instead
                    """
        )
        raise NotImplementedError(None, err_msg)

    if hasattr(args, "lon1"):
        if (args.lon1 is None) != (args.lon2 is None):
            err_msg = textwrap.dedent(
                """\
                        \n ------------------------------------
                        \nERROR: If providing --lon1, you must also provide --lon2
                        """
            )
            raise argparse.ArgumentError(None, err_msg)
        if args.lon1 is not None:
            check_lon1_lt_lon2(args.lon1, args.lon2, args.lon_type)

    return args


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

    nl_datm_base = os.path.join(cesmroot, "components/cdeps/datm/cime_config/user_nl_datm_streams")
    nl_datm = os.path.join(user_mods_dir, "user_nl_datm_streams")
    with open(nl_datm_base, "r") as base_file, open(nl_datm, "w") as user_file:
        for line in base_file:
            user_file.write(line)


def determine_num_pft(crop):
    """
    A simple function to determine the number of pfts.

    Args:
        crop (bool): crop flag denoting if we are using crop

    Raises:

    Returns:
        num_pft (int) : number of pfts for surface dataset
    """
    if crop:
        num_pft = "78"
    else:
        num_pft = "16"
    logger.debug("crop_flag = %s => num_pft = %s", str(crop), num_pft)
    return num_pft


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

    if args.inputdatadir == "defaults.cfg":
        clmforcingindir = defaults.get("main", "clmforcingindir")
    else:
        clmforcingindir = args.inputdatadir

    if not os.path.isdir(clmforcingindir):
        logger.info("clmforcingindir does not exist: %s", clmforcingindir)
        abort("inputdata directory does not exist")

    # DATM data
    # TODO Issue #2960: Make datm_type a user option at the command
    # line. For reference, this option affects three .cfg files:
    #      tools/site_and_regional/default_data_1850.cfg
    #      tools/site_and_regional/default_data_2000.cfg
    #      python/ctsm/test/testinputs/default_data.cfg
    datm_type = "datm_crujra"  # also available: datm_type = "datm_gswp3"
    dir_output_datm = "datmdata"
    dir_input_datm = os.path.join(clmforcingindir, defaults.get(datm_type, "dir"))
    if args.create_datm:
        if not os.path.isdir(os.path.join(args.out_dir, dir_output_datm)):
            os.mkdir(os.path.join(args.out_dir, dir_output_datm))
        logger.info("dir_input_datm : %s", dir_input_datm)
        logger.info("dir_output_datm: %s", os.path.join(args.out_dir, dir_output_datm))

    # if the crop flag is on - we need to use a different land use and surface data file
    num_pft = determine_num_pft(args.crop_flag)

    fsurf_in = defaults.get("surfdat", "surfdat_" + num_pft + "pft")
    fluse_in = defaults.get("landuse", "landuse_" + num_pft + "pft")
    if args.out_surface:
        fsurf_out = args.out_surface
    else:
        fsurf_out = None

    file_dict = {
        "main_dir": clmforcingindir,
        "fdomain_in": defaults.get("domain", "file"),
        "fsurf_dir": os.path.join(
            clmforcingindir,
            os.path.join(defaults.get("surfdat", "dir")),
        ),
        "fluse_dir": os.path.join(
            clmforcingindir,
            os.path.join(defaults.get("landuse", "dir")),
        ),
        "fsurf_in": fsurf_in,
        "fsurf_out": fsurf_out,
        "fluse_in": fluse_in,
        "datm_tuple": DatmFiles(
            dir_input_datm,
            dir_output_datm,
            defaults.get(datm_type, "domain"),
            defaults.get(datm_type, "solardir"),
            defaults.get(datm_type, "precdir"),
            defaults.get(datm_type, "tpqwdir"),
            defaults.get(datm_type, "solartag"),
            defaults.get(datm_type, "prectag"),
            defaults.get(datm_type, "tpqwtag"),
            defaults.get(datm_type, "solarname"),
            defaults.get(datm_type, "precname"),
            defaults.get(datm_type, "tpqwname"),
        ),
    }

    return file_dict


def subset_point(args, file_dict: dict):
    """
    Subsets surface, domain, land use, and/or DATM files at a single point
    """

    logger.info("----------------------------------------------------------------------------")
    logger.info("This script extracts a single point from the global CTSM datasets.")

    num_pft = int(determine_num_pft(args.crop_flag))

    # --  Create SinglePoint Object
    single_point = SinglePointCase(
        plat=args.plat,
        plon=args.plon,
        site_name=args.site_name,
        create_domain=args.create_domain,
        create_surfdata=args.create_surfdata,
        create_landuse=args.create_landuse,
        create_datm=args.create_datm,
        create_user_mods=args.create_user_mods,
        dom_pft=args.dom_pft,
        evenly_split_cropland=args.evenly_split_cropland,
        pct_pft=args.pct_pft,
        num_pft=num_pft,
        cth=args.cth,
        cbh=args.cbh,
        include_nonveg=args.include_nonveg,
        uni_snow=args.uni_snow,
        cap_saturation=args.cap_saturation,
        out_dir=args.out_dir,
        overwrite=args.overwrite,
    )

    logger.debug(single_point)

    # --  Create CTSM surface data file
    if single_point.create_surfdata:
        single_point.create_surfdata_at_point(
            file_dict["fsurf_dir"],
            file_dict["fsurf_in"],
            args.user_mods_dir,
            specify_fsurf_out=file_dict["fsurf_out"],
        )

    # --  Create CTSM transient landuse data file
    if single_point.create_landuse:
        single_point.create_landuse_at_point(
            file_dict["fluse_dir"], file_dict["fluse_in"], args.user_mods_dir
        )

    # --  Create single point atmospheric forcing data
    if single_point.create_datm:
        # subset DATM domain file
        single_point.create_datmdomain_at_point(file_dict["datm_tuple"])

        # subset the DATM data
        nl_datm = os.path.join(args.user_mods_dir, "user_nl_datm_streams")
        single_point.create_datm_at_point(
            file_dict["datm_tuple"], args.datm_syr, args.datm_eyr, nl_datm
        )

    # -- Write shell commands
    if single_point.create_user_mods:
        shell_commands_file = os.path.join(args.user_mods_dir, "shell_commands")
        single_point.write_shell_commands(shell_commands_file, args.datm_syr, args.datm_eyr)

    logger.info("Successfully ran script for single point.")


def _set_up_regional_case(args):
    """
    Set up regional case
    """
    region = RegionalCase(
        lat1=args.lat1,
        lat2=args.lat2,
        lon1=args.lon1,
        lon2=args.lon2,
        reg_name=args.reg_name,
        create_domain=args.create_domain,
        create_surfdata=args.create_surfdata,
        create_landuse=args.create_landuse,
        create_datm=args.create_datm,
        create_user_mods=args.create_user_mods,
        create_mesh=args.create_mesh,
        out_dir=args.out_dir,
        overwrite=args.overwrite,
    )
    logger.debug(region)
    return region


def subset_region(args, file_dict: dict):
    """
    Subsets surface, domain, land use, and/or DATM files for a region
    """

    logger.info("----------------------------------------------------------------------------")
    logger.info("This script extracts a region from the global CTSM datasets.")

    # --  Create Region Object
    region = _set_up_regional_case(args)

    # --  Create CTSM domain file
    if region.create_domain:
        region.create_domain_at_reg(file_dict["main_dir"], file_dict["fdomain_in"])

    # --  Create CTSM surface data file
    if region.create_surfdata:
        region.create_surfdata_at_reg(
            file_dict["fsurf_dir"],
            file_dict["fsurf_in"],
            args.user_mods_dir,
            specify_fsurf_out=file_dict["fsurf_out"],
        )

    # --  Create CTSM transient landuse data file
    if region.create_landuse:
        region.create_landuse_at_reg(
            file_dict["fluse_dir"], file_dict["fluse_in"], args.user_mods_dir
        )

    # -- Write shell commands
    if region.create_user_mods:
        region.write_shell_commands(os.path.join(args.user_mods_dir, "shell_commands"))

        print("\nFor running this regional case with the created user_mods : ")
        print(
            "./create_newcase --case case --res CLM_USRDAT --compset I2000Clm60BgcCrop",
            "--run-unsupported --user-mods-dirs ",
            args.user_mods_dir,
            "\n\n",
        )

    logger.info("Successfully ran script for a regional case.")


def _detect_lon_type(lon_in):
    if lon_in < 0:
        lon_type = 180
    elif lon_in > 180:
        lon_type = 360
    else:
        msg = "When providing an ambiguous longitude, you must specify --lon-type 180 or 360"
        raise argparse.ArgumentTypeError(msg)
    return lon_type


def process_args(args):
    """
    Process arguments after parsing
    """
    # process logging args (i.e. debug and verbose)
    process_logging_args(args)

    # process longitude args
    lon_args = [var for var in ["plon", "lon1", "lon2"] if hasattr(args, var)]
    lon_arg_values = [getattr(args, var) is not None for var in lon_args]
    if any(lon_arg_values):
        if args.lon_type is None:
            if hasattr(args, "plon"):
                args.lon_type = _detect_lon_type(args.plon)
            else:
                lon1_type = _detect_lon_type(args.lon1)
                lon2_type = _detect_lon_type(args.lon2)
                if lon1_type != lon2_type:
                    raise argparse.ArgumentTypeError(
                        "--lon1 and --lon2 seem to be of different types"
                    )
                args.lon_type = lon1_type
        for var in lon_args:
            val = getattr(args, var)
            if val is None:
                continue
            setattr(args, var, Longitude(val, args.lon_type))
    return args


def main():
    """
    Calls functions that subset surface, landuse, domain, and/or DATM files for a region or a
    single point.
    """

    # --------------------------------- #
    # add logging flags from ctsm_logging
    setup_logging_pre_config()
    parser = get_parser()
    args = parser.parse_args()

    # --------------------------------- #
    args = check_args(args)

    # --------------------------------- #
    # parse defaults file
    cesmroot = path_to_ctsm_root()
    defaults = configparser.ConfigParser()
    defaults.read(args.config_file)

    # --------------------------------- #
    myname = getuser()
    pwd = os.getcwd()
    logger.info("User = %s", myname)
    logger.info("Current directory = %s", pwd)

    # --------------------------------- #
    # create files and folders necessary and return dictionary of file/folder locations
    file_dict = setup_files(args, defaults, cesmroot)

    if args.run_type == "point":
        subset_point(args, file_dict)
    elif args.run_type == "region":
        subset_region(args, file_dict)
