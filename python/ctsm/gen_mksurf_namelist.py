# 2020-11-08                Negin Sobhani

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This Python script is part of the simplified toolchain for creating
the surface dataset for ctsm cases.
This script should be used as the first step of the new toolchain. 
It will automatically create namelist (control  file) that is 
needed for creating surface dataset and requisite intermediate files for
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
 
To remove NPL(ncar_pylib) from your environment on Cheyenne/Casper:
    deactivate
-------------------------------------------------------------------
"""

# TODO (NS)

# -[x] Add default values in the help page.
# -[x] Add info for help page note for end_year -- by default is start_year
# -[x] Possibly remove year --years and range options
#      Currently comment them out.

# -[x] maybe a verbose option and removing debug
# -[x] --debug mode is not working...

# -[ ] add error check for hi-res and years if they are 1850 and 2005.
# -[ ] hirespft data only for 2005? add error-check

# -[x] different path for each range of years for transient cases.
#      default should be picked based on the year. 1850 - 2015 -->
#       /glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/pftcftdynharv.0.25x0.25.LUH2.histsimyr1850-2015.c170629/
#      850-1850 -->
#       pftcftdynharv.0.25x0.25.LUH2.histsimyr0850-1849.c171012


#  Import libraries
from __future__ import print_function

import os
import sys
import logging
import argparse
import subprocess

# -- import local classes for this script
from ctsm.toolchain.ctsm_case import CtsmCase

# -- import ctsm logging flags
from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)

logger = logging.getLogger(__name__)

## valid options for resolution and SSP scenarios:
valid_opts = {
    "res": [
        "512x1024",
        "360x720cru",
        "128x256",
        "64x128",
        "48x96",
        "94x192",
        "0.23x0.31",
        "0.47x0.63",
        "0.9x1.25",
        "1.9x2.5",
        "2.5x3.33",
        "4x5",
        "10x15",
        "0.125nldas2",
        "5x5_amazon",
        "1x1_camdenNJ",
        "1x1_vancouverCAN",
        "1x1_mexicocityMEX",
        "1x1_asphaltjungleNJ",
        "1x1_brazil,1x1_urbanc_alpha",
        "1x1_numaIA,1x1_smallvilleIA",
        "0.1x0.1",
        "0.25x0.25",
        "0.5x0.5",
        "3x3min",
        "5x5min",
        "10x10min",
        "0.33x0.33",
        "0.125x0.125",
        "ne4np4,ne16np4",
        "ne30np4.pg2",
        "ne30np4.pg3",
        "ne30np4",
        "ne60np4",
        "ne120np4",
    ],
    "ssp_rcp": [
        "hist",
        "SSP1-2.6",
        "SSP3-7.0",
        "SSP5-3.4",
        "SSP2-4.5",
        "SSP1-1.9",
        "SSP4-3.4",
        "SSP4-6.0",
        "SSP5-8.5",
    ],
}


def get_parser():
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        "--sy",
        "--start-year",
        help="Simulation start year. [default: %(default)s] ",
        action="store",
        dest="start_year",
        required=False,
        type=start_year_type,
        default=2000,
    )
    parser.add_argument(
        "--ey",
        "--end-year",
        help="Simulation end year.  [default: start_year] ",
        action="store",
        dest="end_year",
        required=False,
        type=int,
    )
    parser.add_argument(
        "--glc-nec",
        help="""
            Number of glacier elevation classes to use.
            [default: %(default)s]
            """,
        action="store",
        dest="glc_nec",
        type=glc_nec_type,
        default="10",
    )
    parser.add_argument(
        "--rundir",
        help="""
            Directory to run in.
            [default: %(default)s]
            """,
        action="store",
        dest="run_dir",
        required=False,
        default=os.getcwd(),
    )
    parser.add_argument(
        "--ssp-rcp",
        help="""
            Shared Socioeconomic Pathway and Representative
            Concentration Pathway Scenario name(s).
            [default: %(default)s]
            """,
        action="store",
        dest="ssp_rcp",
        required=False,
        choices=valid_opts["ssp_rcp"],
        default="hist",
    )

    ##############################################
    # In mksurfdata.pl these options are -l --dinlc
    # But the group decided --raw_dir is more descriptive.
    # If everyone agrees, the commented out line should be removed.
    # parser.add_argument('-l','--dinlc',  #--raw-dir or --rawdata-dir

    parser.add_argument(
        "--raw-dir",
        "--rawdata-dir",
        help="""
            /path/of/root/of/input/data',
            [default: %(default)s]
            """,
        action="store",
        dest="input_path",
        default="/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/",
    )
    parser.add_argument(
        "--vic",
        help="""
            Flag for adding the fields required for the VIC model.
            """,
        action="store_true",
        dest="vic_flag",
        default=False,
    )
    parser.add_argument(
        "--glc",
        help="""
            Flag for adding the optional 3D glacier fields for verification of the glacier model.
            """,
        action="store_true",
        dest="glc_flag",
        default=False,
    )
    parser.add_argument(
        "--hirespft",
        help="""
            If you want to use the high-resolution pft dataset rather
            than the default lower resolution dataset.
            (Low resolution is at quarter-degree, high resolution at 3-minute)
            [Note: hires only available for 1850 and 2005.]
            """,
        action="store_true",
        dest="hres_flag",
        default=False,
    )
    parser.add_argument(
        "--nocrop",
        help="""
            Create datasets with the extensive list of prognostic crop types.
            """,
        action="store_false",
        dest="crop_flag",
        default=True,
    )
    parser.add_argument(
        "-f",
        "--fast",
        help="Toggle fast mode which does not user the large mapping file",
        action="store_true",
        dest="fast_flag",
        default=False,
    )
    parser.add_argument(
        "-r",
        "--res",
        help="""
            Resolution is the supported resolution(s) to use for files.
            [default: %(default)s]
            """,
        action="store",
        dest="res",
        choices=valid_opts["res"],
        required=False,
        default="4x5",
    )
    return parser


# -- types for this parser


def glc_nec_type(glc):
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
    glc = int(glc)
    if (glc <= 0) or (glc >= 100):
        raise argparse.ArgumentTypeError("ERROR: glc_nec must be between 1 and 99.")
    return glc.__str__()


def start_year_type(x):
    """
    Function for defining acceptable start_year input.

    Args:
        x (str) : start_year string from command line args.

    Raises:
        Error if value of start_year is not in the range
        of 850-2105.

    Returns:
        x (int) : Acceptable start_year value.
    """
    x = int(x)
    if (x < 850) or (x > 2105):
        raise argparse.ArgumentTypeError(
            "ERROR: Simulation start year should be between 850 and 2105."
        )
    return x


def main():
    # -- add logging flags from ctsm_logging
    setup_logging_pre_config()
    parser = get_parser()
    add_logging_args(parser)

    args = parser.parse_args()
    process_logging_args(args)

    res = args.res
    glc_nec = args.glc_nec
    input_path = args.input_path
    ssp_rcp = args.ssp_rcp
    crop_flag = args.crop_flag
    vic_flag = args.vic_flag
    glc_flag = args.glc_flag
    hres_flag = args.hres_flag

    start_year = args.start_year
    end_year = args.end_year

    # -- determine end_year if not given as an argument:
    if not end_year:
        end_year = start_year

    # -- check if the input path exist
    if not os.path.exists(input_path):
        sys.exit(
            "ERROR: \n"
            + "\t raw_dir does not exist on this machine. \n"
            + "\t Please point to the correct raw_dir using --raw-dir"
            + "or --rawdata-dir flags."
        )

    ctsm_case = CtsmCase(
        res,
        glc_nec,
        ssp_rcp,
        crop_flag,
        input_path,
        vic_flag,
        glc_flag,
        start_year,
        end_year,
    )

    logger.info("--------------------------")
    logger.info(" ctsm case : %s", ctsm_case)
    logger.info("--------------------------")

    ctsm_case.create_namelist_file()


if __name__ == "__main__":
    main()
