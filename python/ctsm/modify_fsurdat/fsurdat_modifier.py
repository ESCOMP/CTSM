"""
Run this code by using the following wrapper script:
tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter

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

    # Accepted values
    # TODO Read control/namelist file rather than command-line options/arguments

    parser.add_argument('--dom_nat_pft',
                help='Optional non-crop PFT to be set to 100% everywhere ' \
                     '(integer from 0 to natpft; see natpft in your ' \
                     'fsurdat_in file). [default: %(default)s]',
                action="store",
                dest="dom_nat_pft",
                type=int,
                choices=range(0, 15),  # maximum natveg = 14
                default=-999)
    parser.add_argument('--lai',
                help='Optional non-crop PFT leaf area index value to be ' \
                     'set when setting --dom_nat_pft to non-bare soil '
                     'only when also setting the lon/lat boundaries '
                     '(float). [default: %(default)s]',
                action="store",
                dest="lai",
                type=float,
                default=0)
    parser.add_argument('--sai',
                help='Optional non-crop PFT stem area index value to be ' \
                     'set when setting --dom_nat_pft to non-bare soil '
                     'only when also setting the lon/lat boundaries '
                     '(float). [default: %(default)s]',
                action="store",
                dest="sai",
                type=float,
                default=0)
    parser.add_argument('--hgt_top',
                help='Optional non-crop PFT top of canopy height to be ' \
                     'set when setting --dom_nat_pft to non-bare soil '
                     'only when also setting the lon/lat boundaries '
                     '(float, meters). [default: %(default)s]',
                action="store",
                dest="hgt_top",
                type=float,
                default=0)
    parser.add_argument('--hgt_bot',
                help='Optional non-crop PFT bottop of canopy height to be ' \
                     'set when setting --dom_nat_pft to non-bare soil '
                     'only when also setting the lon/lat boundaries '
                     '(float, meters). [default: %(default)s]',
                action="store",
                dest="hgt_bot",
                type=float,
                default=0)
    parser.add_argument('--std_elev',
                help='Optional STD_ELEV value to specify uniform snowpack ' \
                     'everywhere (integer from 0 to 100 m). [default: %(default)s]',
                action="store",
                dest="std_elev",
                type=int,
                choices=range(0, 101),  # TODO reasonable range?
                default=-999)
    parser.add_argument('--zero_nonveg',
                help='Optional: sets non-vegetation landunits to 0 globally. ' \
                     'Redundant if defining new land mask using --lnd_lat_1, ' \
                     '--lnd_lat_2, --lnd_lon_1, and --lnd_lon_2 [default: %(default)s]',
                action="store_true",
                dest="zero_nonveg")
    parser.add_argument('--max_sat_area',
                help='Optional maximum fractional saturated area, aka. FMAX ' \
                     'to set to everywhere (from 0 to 1). [default: %(default)s]',
                action="store",
                dest="max_sat_area",
                type=float,
                default=-999)
    parser.add_argument('--lnd_lat_1',
                help='Optional southernmost latitude for rectangle (integer ' \
                     '-90 to 90). If lnd_lat_1 > lnd_lat_2 the code creates ' \
                     'two rectangles, one in the north and one in the ' \
                     'south. Required for lnd_lat_2, lnd_lon_1, and ' \
                     'lnd_lon_2 to work together [default: %(default)s]',
                dest='lnd_lat_1',
                action="store",
                type=int,
                choices=range(-90, 91),
                default=-999)
    parser.add_argument('--lnd_lat_2',
                help='Optional northernmost latitude for rectangle (integer ' \
                     '-90 to 90). If lnd_lat_1 > lnd_lat_2 the code creates ' \
                     'two rectangles, one in the north and one in the ' \
                     'south. Required for lnd_lat_1, lnd_lon_1, and ' \
                     'lnd_lon_2 to work together [default: %(default)s]',
                action="store",
                dest="lnd_lat_2",
                type=int,
                choices=range(-90, 91),
                default=-999)
    parser.add_argument('--lnd_lon_1',
                help='Optional minimum longitude for rectangle (integer 0 ' \
                     'to 360). If lnd_lon_1 > lnd_lon_2 the rectangle wraps ' \
                     'around the 0-degree meridian. Required for lnd_lon_2, ' \
                     'lnd_lat_1, and lnd_lat_2 to work together [default: %(default)s]',
                action="store",
                dest="lnd_lon_1",
                type=int,
                choices=range(0, 361),
                default=-999)
    parser.add_argument('--lnd_lon_2',
                help='Optional maximum longitude for rectangle (integer 0 ' \
                     'to 360). If lnd_lon_1 > lnd_lon_2 the rectangle wraps ' \
                     'around the 0-degree meridian. Required for lnd_lon_1, ' \
                     'lnd_lat_1, and lnd_lat_2 to work together [default: %(default)s]',
                action="store",
                dest="lnd_lon_2",
                type=int,
                choices=range(0, 361),
                default=-999)
    parser.add_argument('--fsurdat_in',
                help='Path and name of input surface dataset [default: %(default)s]',
                action="store",
                dest="fsurdat_in",
                type=str,
                # TODO: Require user to enter this and remove the default
                default="/glade/p/cesmdata/cseg/inputdata" \
                        "/lnd/clm2/surfdata_map/" \
                        "surfdata_10x15_78pfts_CMIP6_simyr1850_c170824.nc")
    parser.add_argument('--fsurdat_out',
                help='Required path and name of output surface dataset [default: %(default)s]',
                action="store",
                dest="fsurdat_out",
                type=str,
                required=True)  # TODO Doesn't return 'this is required' error

    return parser


def main ():
    """
    Description
    -----------
    Calls various functions that modify an fsurdat (surface dataset)
    """

    # Parse arguments from the command line
    args = get_parser().parse_args()

    # Create ModifyFsurdat object
    modify_fsurdat = ModifyFsurdat(args.fsurdat_in)

    # ------------------------------
    # modify surface data properties
    # ------------------------------

    # 1) Set fsurdat variables in a rectangle
    # Note that the land/ocean mask gets specified in the domain file for
    # MCT or the ocean mesh files for NUOPC. The function set_a_rectangle
    # can specify fsurdat variables inside a box but it cannot
    # change which points will run as land and which as ocean.
    if args.lnd_lon_1 == -999 or args.lnd_lat_1 == -999 or \
       args.lnd_lon_2 == -999 or args.lnd_lat_2 == -999:
        warning_msg = 'Warning: One or more of the optional arguments ' \
                      'lnd_lon_1, lnd_lon_2, lnd_lat_1, lnd_lat_2 were not ' \
                      'set, so all four are ignored.'
        print(warning_msg)  # TODO Use logging for this statement

        # 1a) If not setting in a rectangle, can set user-selected
        # dom_nat_pft to 100% globally.
        if args.dom_nat_pft != -999:
            modify_fsurdat.dom_nat_pft(args.dom_nat_pft)
        # 1b) If not setting in a rectangle, can set non-vegetation
        # landunits to zero globally.
        if args.zero_nonveg:
            modify_fsurdat.zero_nonveg()
        # 1c) If not modifying a rectangle, can create uniform snowpack
        # by setting STD_ELEV to a constant everywhere
        if args.std_elev != -999:
            modify_fsurdat.std_elev(args.std_elev)
        # 1d) If not modifying a rectangle, can set max_sat_area (FMAX)
        # to a constant everywhere
        if args.max_sat_area != -999:
            modify_fsurdat.max_sat_area(args.max_sat_area)
    else:
        modify_fsurdat.set_a_rectangle(lon_in_1=args.lnd_lon_1,
                                       lon_in_2=args.lnd_lon_2,
                                       lat_in_1=args.lnd_lat_1,
                                       lat_in_2=args.lnd_lat_2,
                                       dom_nat_pft=args.dom_nat_pft,
                                       lai=args.lai,
                                       sai=args.sai,
                                       hgt_top=args.hgt_top,
                                       hgt_bot=args.hgt_bot)

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    modify_fsurdat.write_output(args.fsurdat_in, args.fsurdat_out)

    sys.exit('SUCCESS')
