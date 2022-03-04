"""
Run this code by using the following wrapper script:
tools/modify_mesh_mask/mesh_mask_modifier

The wrapper script includes a full description and instructions.
"""

import logging
import argparse
from configparser import ConfigParser
from ctsm.utils import get_config_value
from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args
from ctsm.modify_mesh_mask.modify_mesh_mask import ModifyMeshMask

logger = logging.getLogger(__name__)

def main ():
    """
    Description
    -----------
    Calls function that modifies mesh mask
    """

    # set up logging allowing user control
    setup_logging_pre_config()

    # read the command line argument to obtain the path to the .cfg file
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg_path',
                        help='/path/name.cfg of input file, eg ./modify.cfg')
    add_logging_args(parser)
    args = parser.parse_args()
    process_logging_args(args)
    mesh_mask_modifier(args.cfg_path)

def mesh_mask_modifier(cfg_path):
    """Implementation of mesh_mask_modifier command"""
    # read the .cfg (config) file
    config = ConfigParser()
    config.read(cfg_path)
    section = config.sections()[0]  # name of the first section

    # required: user must set these in the .cfg file
    mesh_mask_in = get_config_value(config=config, section=section,
        item='mesh_mask_in', file_path=cfg_path)
    mesh_mask_out = get_config_value(config=config, section=section,
        item='mesh_mask_out', file_path=cfg_path)

    # required but fallback values available for variables omitted
    # entirely from the .cfg file
    lnd_lat_1 = get_config_value(config=config, section=section,
        item='lnd_lat_1', file_path=cfg_path, convert_to_type=float)
    lnd_lat_2 = get_config_value(config=config, section=section,
        item='lnd_lat_2', file_path=cfg_path, convert_to_type=float)
    lnd_lon_1 = get_config_value(config=config, section=section,
        item='lnd_lon_1', file_path=cfg_path, convert_to_type=float)
    lnd_lon_2 = get_config_value(config=config, section=section,
        item='lnd_lon_2', file_path=cfg_path, convert_to_type=float)

    # not required: user may set these in the .cfg file
    landmask_file = get_config_value(config=config, section=section,
        item='landmask_file', file_path=cfg_path, can_be_unset=True)

    # Create ModifyMeshMask object
    modify_mesh_mask = ModifyMeshMask.init_from_file(mesh_mask_in,
        lnd_lon_1, lnd_lon_2, lnd_lat_1, lnd_lat_2, landmask_file)

    # ----------------
    # modify mesh mask
    # ----------------

    # Modify mesh mask in a rectangle that could be global (default)
    # TODO Use one of these as a template. Look for mesh_mask files
    # that are gridded rather than gx1/gx3/etc.
    modify_mesh_mask.set_idealized()  # set 2D variables
    modify_mesh_mask.setvar_lev0('FMAX', max_sat_area)

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    modify_mesh_mask.write_output(mesh_mask_in, mesh_mask_out)
