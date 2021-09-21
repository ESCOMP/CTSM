"""
Run this code by using the following wrapper script:
tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import sys
from getpass import getuser
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from ctsm.utils import abort
from ctsm.modify_fsurdat.modify_fsurdat import ModifyFsurdat


def get_parser():

    """
    Description
    -----------
    Get parser object for this script.
    """

    parser = ArgumentParser(description=__doc__,
                           formatter_class=RawDescriptionHelpFormatter)

    parser.print_usage = parser.print_help

    parser.add_argument('--dom_pft',
                help='Enter PFT that you wish to set to 100% everywhere [default: %(default)s] ',
                action="store",
                dest="dom_pft",
                type=int,
                default=-999)
    parser.add_argument('--uni_snow',
                help='Create uniform snowpack. [default: %(default)s]',
                action="store_true",
                dest="uni_snowpack")
    parser.add_argument('--zero_nonveg',
                help='Set all non-vegetation landunits to zero. [default: %(default)s]',
                action="store_true",
                dest="zero_nonveg_lu")
    parser.add_argument('--max_sat_area',
                help='Enter maximum fractional saturated area (FMAX from ' + \
                     '0 to 1). [default: %(default)s]',
                action="store",
                dest="max_sat_area",
                type=float,
                default=-999)
    parser.add_argument('--fsurdat_in',
                help='Input surface dataset. [default: %(default)s]',
                action="store",
                dest="fsurdat_in",
                type=str,
                # TODO: Require user to enter this and remove the default
                default="/glade/p/cesmdata/cseg/inputdata/lnd/clm2/" + \
                        "surfdata_map/release-clm5.0.18/" \
                        "surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr2000_c190214.nc")
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

    # set dom_pft to 100% everywhere
    min_dom_pft = int(min(modify_fsurdat.file.lsmpft))
    max_dom_pft = int(max(modify_fsurdat.file.lsmpft))
    if args.dom_pft >= min_dom_pft and args.dom_pft <= max_dom_pft:
        modify_fsurdat.dom_pft(args.dom_pft)
    else:
        errmsg = 'Argument --dom_pft can only be set to values from ' + \
                 str(min_dom_pft) + ' to ' + str(max_dom_pft)
        abort(errmsg)

    # set all non-vegetation landunits to zero
    if args.zero_nonveg_lu:
        modify_fsurdat.zero_nonveg_lu()

    # create uniform snowpack
    if args.uni_snowpack:
        modify_fsurdat.uni_snowpack()

    # set max_sat_area to a constant everywhere
    min_max_sat_area = 0
    max_max_sat_area = 1
    if args.max_sat_area >= min_max_sat_area and \
       args.max_sat_area <= max_max_sat_area:
        modify_fsurdat.max_sat_area(args.max_sat_area)
    else:
        errmsg = 'Argument --max_sat_area can only be set to values from ' + \
                 str(min_max_sat_area) + ' to ' + str(max_max_sat_area)
        abort(errmsg)

    # output CTSM surface data file
    modify_fsurdat.write_output(args.fsurdat_in, args.fsurdat_out)

    sys.exit('SUCCESS')
