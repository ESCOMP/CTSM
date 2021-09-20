"""
Run this code by using the following wrapper script:
tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import sys
from getpass import getuser
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from ctsm.modify_fsurdat.modify_fsurdat import ModifyFsurdat


def get_parser():

    """
    Get parser object for this script.

    Command-line inputs
    - fsurdat_in: input file (str)
    - fsurdat_out: output file (str)
    - variable: variable to modify (str)
    - value: value assigned to the variable to be modified (float)
    """

    parser = ArgumentParser(description=__doc__,
                           formatter_class=RawDescriptionHelpFormatter)

    parser.print_usage = parser.print_help

    parser.add_argument('--dom_pft',
                help='Dominant PFT if overwrite_single_pft = .true. [default: %(default)s] ',
                action="store",
                dest="dom_pft",
                type=int)
    parser.add_argument('--uni_snow',
                help='Turn on the flag to create uniform snowpack. [default: %(default)s]',
                action="store_true",
                dest="uni_snowpack")
    parser.add_argument('--overwrite_single_pft',
                help='Turn on flag to make whole grid 100% single PFT. [default: %(default)s]',
                action="store_true",
                dest="overwrite_single_pft")
    parser.add_argument('--zero_nonveg',
                help='Set all non-vegetation landunits to zero. [default: %(default)s]',
                action="store_true",
                dest="zero_nonveg_lu")
    parser.add_argument('--no_saturation_excess',
                help='Turn off saturation excess. [default: %(default)s]',
                action="store_true",
                dest="no_saturation_excess")
    parser.add_argument('--fsurdat_in',
                help='Input surface dataset. [default: %(default)s]',
                action="store",
                dest="fsurdat_in",
                type=str,
                # TODO: Require user to enter this and remove the default
                default="/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/release-clm5.0.18/surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr2000_c190214.nc")
    parser.add_argument('--fsurdat_out',
                help='Output surface dataset. [default: %(default)s]',
                action="store",
                dest="fsurdat_out",
                type=str,
                # TODO: Require user to enter this and remove the default
                default="/glade/scratch/" + getuser() +
                          "/surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr2000_c190214_modified.nc")

    return parser


def main ():
    """
    Description
    -----------
    Calls various functions that modify an fsurdat (surface dataset)
    """

    # Parse arguments from the command line
    args = get_parser().parse_args()

    # create ModifyFsurdat object
    modify_fsurdat = ModifyFsurdat(args.fsurdat_in)

    # modify surface data properties
    if args.overwrite_single_pft:
        modify_fsurdat.overwrite_single_pft(args.dom_pft)
    if args.zero_nonveg_lu:
        modify_fsurdat.zero_nonveg_lu()
    if args.uni_snowpack:
        modify_fsurdat.uni_snowpack()
    if args.no_saturation_excess:
        modify_fsurdat.no_saturation_excess()

    # output CTSM surface data file
    modify_fsurdat.write_output(args.fsurdat_in, args.fsurdat_out)

    sys.exit('SUCCESS')
