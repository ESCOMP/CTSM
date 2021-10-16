"""
Run this code by using the following wrapper script:
tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import sys
from configparser import ConfigParser
from ctsm.utils import get_config_value, select_value
from ctsm.modify_fsurdat.modify_fsurdat import ModifyFsurdat


def main ():
    """
    Description
    -----------
    Calls various functions that modify an fsurdat (surface dataset)
    """

    # read the config file
    config = ConfigParser()
    cfg_path = './modify.cfg'
    config.read(cfg_path)

    # required: user must set these in ./modify.cfg
    fsurdat_in = get_config_value(config, 'modify_input', 'fsurdat_in',
                                  cfg_path)
    fsurdat_out = get_config_value(config, 'modify_input', 'fsurdat_out',
                                   cfg_path)

    # not required: user may set these in ./modify.cfg
    # TODO Error checking not complete for these. Should it be done here
    # while reading or later when using them?
    temp = get_config_value(config, 'modify_input', 'idealized', cfg_path,
                            allowed_values=['True','False','UNSET'])
    idealized = select_value(var=temp, default=False, type_of_var=bool)

    temp = get_config_value(config, 'modify_input', 'lnd_lat_1', cfg_path)
    lnd_lat_1 = select_value(var=temp, default=-90, type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'lnd_lat_2', cfg_path)
    lnd_lat_2 = select_value(var=temp, default=90, type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'lnd_lon_1', cfg_path)
    lnd_lon_1 = select_value(var=temp, default=0, type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'lnd_lon_2', cfg_path)
    lnd_lon_2 = select_value(var=temp, default=360, type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'dom_nat_pft', cfg_path,
        allowed_values=['0','1','2','3','4','5','6','7','8','9','10','11',
                        '12','13','14','UNSET'])
    dom_nat_pft = select_value(var=temp, default=None, type_of_var=int)

    temp = get_config_value(config, 'modify_input', 'lai', cfg_path)
    lai = select_value(var=temp, default=[], type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'sai', cfg_path)
    sai = select_value(var=temp, default=[], type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'hgt_top', cfg_path)
    hgt_top = select_value(var=temp, default=[], type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'hgt_bot', cfg_path)
    hgt_bot = select_value(var=temp, default=[], type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'soil_color', cfg_path,
        allowed_values=['1','2','3','4','5','6','7','8','9','10','11','12',
                        '13','14','15','16','17','18','19','20','UNSET'])
    soil_color = select_value(var=temp, default=None, type_of_var=int)

    temp = get_config_value(config, 'modify_input', 'std_elev', cfg_path)
    std_elev = select_value(var=temp, default=None, type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'max_sat_area', cfg_path)
    max_sat_area = select_value(var=temp, default=None, type_of_var=float)

    temp = get_config_value(config, 'modify_input', 'zero_nonveg', cfg_path,
                            allowed_values=['True','False','UNSET'])
    zero_nonveg = select_value(var=temp, default=False, type_of_var=bool)

    # Create ModifyFsurdat object
    modify_fsurdat = ModifyFsurdat(fsurdat_in)

    # ------------------------------
    # modify surface data properties
    # ------------------------------

    # Set fsurdat variables in a rectangle that could be global (default).
    # Note that the land/ocean mask gets specified in the domain file for
    # MCT or the ocean mesh files for NUOPC. The function set_in_rectangle
    # can specify fsurdat variables inside a box but it cannot
    # change which points will run as land and which as ocean.
    modify_fsurdat.set_in_rectangle(_idealized=idealized,
                                    _lon_in_1=lnd_lon_1,
                                    _lon_in_2=lnd_lon_2,
                                    _lat_in_1=lnd_lat_1,
                                    _lat_in_2=lnd_lat_2,
                                    _dom_nat_pft=dom_nat_pft,
                                    _lai=lai,
                                    _sai=sai,
                                    _hgt_top=hgt_top,
                                    _hgt_bot=hgt_bot,
                                    _zero_nonveg=zero_nonveg,
                                    _std_elev=std_elev,
                                    _soil_color=soil_color,
                                    _max_sat_area=max_sat_area)

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    modify_fsurdat.write_output(fsurdat_in, fsurdat_out)

    sys.exit('SUCCESS')
