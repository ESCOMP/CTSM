"""
Run this code by using the following wrapper script:
tools/modify_input_files/mesh_mask_modifier

The wrapper script includes a full description and instructions.
"""

import os
import logging
import argparse
from configparser import ConfigParser

from ctsm.utils import abort, write_output
from ctsm.config_utils import get_config_value
from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args
from ctsm.modify_input_files.modify_mesh_mask import ModifyMeshMask

logger = logging.getLogger(__name__)


def main():
    """
    Description
    -----------
    Calls function that modifies mesh mask
    """

    # set up logging allowing user control
    setup_logging_pre_config()

    # read the command line argument to obtain the path to the .cfg file
    parser = argparse.ArgumentParser()
    parser.add_argument("cfg_path", help="/path/name.cfg of input file, eg ./modify.cfg")
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
    mesh_mask_in = get_config_value(
        config=config, section=section, item="mesh_mask_in", file_path=cfg_path
    )
    mesh_mask_out = get_config_value(
        config=config, section=section, item="mesh_mask_out", file_path=cfg_path
    )
    landmask_file = get_config_value(
        config=config, section=section, item="landmask_file", file_path=cfg_path
    )
    lat_dimname = get_config_value(
        config=config, section=section, item="lat_dimname", file_path=cfg_path
    )
    lon_dimname = get_config_value(
        config=config, section=section, item="lon_dimname", file_path=cfg_path
    )
    lat_varname = get_config_value(
        config=config, section=section, item="lat_varname", file_path=cfg_path
    )
    lon_varname = get_config_value(
        config=config, section=section, item="lon_varname", file_path=cfg_path
    )
    lon_type = get_config_value(config=config, section=section, item="lon_type", file_path=cfg_path)

    # Create ModifyMeshMask object
    modify_mesh_mask = ModifyMeshMask.init_from_file(
        file_in=mesh_mask_in,
        landmask_file=landmask_file,
        lat_dimname=lat_dimname,
        lon_dimname=lon_dimname,
        lat_varname=lat_varname,
        lon_varname=lon_varname,
        lon_type=lon_type,
    )

    # If output file exists, abort before starting work
    if os.path.exists(mesh_mask_out):
        errmsg = "Output file already exists: " + mesh_mask_out
        abort(errmsg)

    # ----------------
    # modify mesh mask
    # ----------------

    # Modify mesh mask
    modify_mesh_mask.set_mesh_mask("elementMask")

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    write_output(modify_mesh_mask.file, mesh_mask_in, mesh_mask_out, "mesh")
