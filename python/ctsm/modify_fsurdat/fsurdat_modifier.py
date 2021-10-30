"""
Run this code by using the following wrapper script:
tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

import sys
import argparse
from configparser import ConfigParser
from ctsm.utils import get_config_value
from ctsm.utils import CONFIG_UNSET
from ctsm.modify_fsurdat.modify_fsurdat import ModifyFsurdat


def main ():
    """
    Description
    -----------
    Calls function that modifies an fsurdat (surface dataset)
    """

    # read the command line argument to obtain the path to the .cfg file
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg_path',
                        help='/path/name.cfg of input file, eg ./modify.cfg')
    args = parser.parse_args()

    # read the .cfg (config) file
    config = ConfigParser()
    config.read(args.cfg_path)
    section = config.sections()[0]  # name of the first section

    # required: user must set these in the .cfg file
    fsurdat_in = get_config_value(config=config, section=section,
        item='fsurdat_in', file_path=args.cfg_path)
    fsurdat_out = get_config_value(config=config, section=section,
        item='fsurdat_out', file_path=args.cfg_path)

    # required but fallback values available for variables omitted
    # entirely from the .cfg file
    idealized = config.getboolean(section, 'idealized', fallback=False)
    zero_nonveg = config.getboolean(section, 'zero_nonveg',
                                    fallback=False)

    lnd_lat_1 = config.getfloat(section, 'lnd_lat_1', fallback=-90)
    lnd_lat_2 = config.getfloat(section, 'lnd_lat_2', fallback=90)
    lnd_lon_1 = config.getfloat(section, 'lnd_lon_1', fallback=0)
    lnd_lon_2 = config.getfloat(section, 'lnd_lon_2', fallback=360)

    # not required: user may set these in the .cfg file
    dom_nat_pft = get_config_value(config=config, section=section,
        item='dom_nat_pft', file_path=args.cfg_path,
        allowed_values=['0','1','2','3','4','5','6','7','8','9','10',
        '11','12','13','14',CONFIG_UNSET], convert_to_type=int)

    lai = get_config_value(config=config, section=section, item='lai',
        file_path=args.cfg_path, is_list=True, convert_to_type=float)
    sai = get_config_value(config=config, section=section, item='sai',
        file_path=args.cfg_path, is_list=True, convert_to_type=float)
    hgt_top = get_config_value(config=config, section=section,
        item='hgt_top', file_path=args.cfg_path, is_list=True,
        convert_to_type=float)
    hgt_bot = get_config_value(config=config, section=section,
        item='hgt_bot', file_path=args.cfg_path, is_list=True,
        convert_to_type=float)

    soil_color = get_config_value(config=config, section=section,
        item='soil_color', file_path=args.cfg_path,
        allowed_values=['1','2','3','4','5','6','7','8','9','10','11',
        '12','13','14','15','16','17','18','19','20',CONFIG_UNSET],
        convert_to_type=int)

    std_elev = get_config_value(config=config, section=section,
        item='std_elev', file_path=args.cfg_path, convert_to_type=float)
    max_sat_area = get_config_value(config=config, section=section,
        item='max_sat_area', file_path=args.cfg_path,
        convert_to_type=float)

    # Create ModifyFsurdat object
    modify_fsurdat = ModifyFsurdat(fsurdat_in, lon_1=lnd_lon_1,
        lon_2=lnd_lon_2, lat_1=lnd_lat_1, lat_2=lnd_lat_2)

    # ------------------------------
    # modify surface data properties
    # ------------------------------

    # Set fsurdat variables in a rectangle that could be global (default).
    # Note that the land/ocean mask gets specified in the domain file for
    # MCT or the ocean mesh files for NUOPC. Here the user may specify
    # fsurdat variables inside a box but cannot change which points will
    # run as land and which as ocean.
    if idealized:
        modify_fsurdat.set_idealized()  # set 2D variables
        # set 3D and 4D variables pertaining to natural vegetation
        modify_fsurdat.set_dom_nat_pft(dom_nat_pft=0, lai=[], sai=[],
                                       hgt_top=[], hgt_bot=[])

    if dom_nat_pft is not None:  # overwrite "idealized" value
        modify_fsurdat.set_dom_nat_pft(dom_nat_pft=dom_nat_pft,
                                       lai=lai, sai=sai,
                                       hgt_top=hgt_top, hgt_bot=hgt_bot)

    if max_sat_area is not None:  # overwrite "idealized" value
        modify_fsurdat.file['FMAX'] = \
         modify_fsurdat.file['FMAX'].where(modify_fsurdat.not_rectangle, other=max_sat_area)

    if std_elev is not None:  # overwrite "idealized" value
        modify_fsurdat.file['STD_ELEV'] = \
         modify_fsurdat.file['STD_ELEV'].where(modify_fsurdat.not_rectangle, other=std_elev)

    if soil_color is not None:  # overwrite "idealized" value
        modify_fsurdat.file['SOIL_COLOR'] = \
         modify_fsurdat.file['SOIL_COLOR'].where(modify_fsurdat.not_rectangle, other=soil_color)

    if zero_nonveg:
        modify_fsurdat.zero_nonveg()

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    modify_fsurdat.write_output(fsurdat_in, fsurdat_out)

    sys.exit('SUCCESS')
