"""
Run this code by using the following wrapper script:
tools/modify_input_files/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

import os
import logging
import argparse
from configparser import ConfigParser

from ctsm.utils import abort, write_output
from ctsm.config_utils import get_config_value, get_config_value_or_array
from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)
from ctsm.modify_input_files.modify_fsurdat import ModifyFsurdat

logger = logging.getLogger(__name__)


def main():
    """
    Description
    -----------
    Calls function that modifies an fsurdat (surface dataset)
    """

    args = fsurdat_modifier_arg_process()
    fsurdat_modifier(args)


def fsurdat_modifier_arg_process():
    """Argument processing for fsurdat_modifier script"""
    # set up logging allowing user control
    setup_logging_pre_config()

    # read the command line argument to obtain the path to the .cfg file
    parser = argparse.ArgumentParser()
    parser.add_argument("cfg_path", help="/path/name.cfg of input file, eg ./modify.cfg")
    parser.add_argument(
        "-i",
        "--fsurdat_in",
        default="UNSET",
        required=False,
        type=str,
        help="The input surface dataset to modify. ",
    )
    parser.add_argument(
        "-o",
        "--fsurdat_out",
        required=False,
        default="UNSET",
        type=str,
        help="The output surface dataset with the modifications. ",
    )
    parser.add_argument(
        "--overwrite",
        required=False,
        default=False,
        action="store_true",
        help="Overwrite the output file if it already exists. ",
    )
    add_logging_args(parser)
    args = parser.parse_args()
    process_logging_args(args)
    # Error checking of arguments
    if not os.path.exists(args.cfg_path):
        abort("Config file does NOT exist: " + str(args.cfg_path))

    return args


def check_no_subgrid_section(config):
    """Check that there isn't a subgrid section when it's processing is turned off"""
    section = "modify_fsurdat_subgrid_fractions"
    if config.has_section(section):
        abort(
            "Config file does have a section: "
            + section
            + " that should NOT be there since it is turned off"
        )


def check_no_varlist_section(config):
    """Check that there isn't a var list section when it's processing is turned off"""
    section = "modify_fsurdat_variable_list"
    if config.has_section(section):
        abort(
            "Config file does have a section: "
            + section
            + " that should NOT be there since it is turned off"
        )


def check_range(var, section, value, minval, maxval):
    """Check that the value is within range"""
    if value < minval or value > maxval:
        abort("Variable " + var + " in " + section + " is out of range of 0 to 100 = " + str(value))


def read_cfg_subgrid(config, cfg_path, numurbl=3):
    """Read the subgrid fraction section from the config file"""
    section = "modify_fsurdat_subgrid_fractions"
    if not config.has_section(section):
        abort("Config file does not have the expected section: " + section)

    subgrid_settings = {}
    var_list = config.options(section)
    valid_list = [
        "pct_natveg",
        "pct_crop",
        "pct_lake",
        "pct_glacier",
        "pct_wetland",
        "pct_urban",
        "pct_ocean",
    ]
    varsum = 0
    for var in var_list:
        if valid_list.count(var) == 0:
            abort(
                "Variable "
                + var
                + " in "
                + section
                + " is not a valid variable name. Valid vars ="
                + str(valid_list)
            )
        # Urban is multidimensional
        if var == "pct_urban":
            vallist = get_config_value(
                config=config,
                section=section,
                item=var,
                file_path=cfg_path,
                is_list=True,
                convert_to_type=float,
            )
            if len(vallist) != numurbl:
                abort("PCT_URBAN is not a list of the expected size of " + str(numurbl))
            # so if a scalar value, must be multiplied # by the density dimension
            for val in vallist:
                check_range(var, section, val, 0.0, 100.0)
                varsum += val
            value = vallist
        else:
            value = get_config_value(
                config=config, section=section, item=var, file_path=cfg_path, convert_to_type=float
            )
            check_range(var, section, value, 0.0, 100.0)
            varsum += value

        subgrid_settings[var.upper()] = value

    if varsum != 100.0:
        abort(
            "PCT fractions in subgrid section do NOT sum to a hundred as they should. Sum = "
            + str(varsum)
        )

    return subgrid_settings


def read_cfg_var_list(config, idealized=True):
    """Read the variable list section from the config file"""
    section = "modify_fsurdat_variable_list"
    if not config.has_section(section):
        abort("Config file does not have the expected section: " + section)

    varlist_settings = {}
    var_list = config.options(section)
    ideal_list = [
        "soil_color",
        "pct_sand",
        "pct_clay",
        "organic",
        "pct_cft",
        "pct_nat_pft",
        "fmax",
        "std_elev",
    ]
    subgrid_list = ["pct_natveg", "pct_crop", "pct_lake", "pct_glacier", "pct_wetland", "pct_urban"]
    # List of variables that should be excluded because they are changed elsewhere,
    # or they shouldn't be changed # Ds, Dsmax, and Ws are excluded because they
    # are of mixed case and we only search for varaibles in lowercase
    # or uppercase and not mixed case.
    monthly_list = [
        "monthly_lai",
        "monthly_sai",
        "monthly_height_top",
        "monthly_height_bot",
        "ds",
        "mxsoil_color",
        "natpft",
        "cft",
        "time",
        "longxy",
        "latixy",
        "dsmax",
        "area",
        "ws",
    ]
    for var in var_list:
        if idealized and ideal_list.count(var) != 0:
            abort(
                var
                + " is a special variable handled in the idealized section."
                + " This should NOT be handled in the variable list section."
                + " Special idealized vars ="
                + str(ideal_list)
            )
        if subgrid_list.count(var) != 0:
            abort(
                var
                + " is a variable handled in the subgrid section."
                + " This should NOT be handled in the variable list section."
                + " Subgrid vars ="
                + str(subgrid_list)
            )
        if monthly_list.count(var) != 0:
            abort(
                var
                + " is a variable handled as part of the dom_pft handling."
                + " This should NOT be handled in the variable list section."
                + " Monthly vars handled this way ="
                + str(monthly_list)
            )
        value = get_config_value_or_array(
            config=config, section=section, item=var, convert_to_type=float
        )
        varlist_settings[var] = value

    return varlist_settings


def modify_optional(
    *,
    modify_fsurdat,
    idealized,
    include_nonveg,
    max_sat_area,
    std_elev,
    soil_color,
    dom_pft,
    evenly_split_cropland,
    lai,
    sai,
    hgt_top,
    hgt_bot,
):
    """Modify the dataset according to the optional settings"""

    # Set fsurdat variables in a rectangle that could be global (default).
    # Note that the land/ocean mask gets specified in
    # the ocean mesh files. Here the user may specify
    # fsurdat variables inside a box but cannot change which points will
    # run as land and which as ocean.
    if idealized:
        modify_fsurdat.set_idealized()  # set 2D variables
        # set 3D and 4D variables pertaining to natural vegetation
        # to default values here; allow override values with the later call
        # to set_dom_pft
        modify_fsurdat.set_dom_pft(dom_pft=0, lai=[], sai=[], hgt_top=[], hgt_bot=[])
        logger.info("idealized complete")

    if max_sat_area is not None:  # overwrite "idealized" value
        modify_fsurdat.setvar_lev0("FMAX", max_sat_area)
        logger.info("max_sat_area complete")

    if std_elev is not None:  # overwrite "idealized" value
        modify_fsurdat.setvar_lev0("STD_ELEV", std_elev)
        logger.info("std_elev complete")

    if soil_color is not None:  # overwrite "idealized" value
        modify_fsurdat.setvar_lev0("SOIL_COLOR", soil_color)
        logger.info("soil_color complete")

    if not include_nonveg:
        modify_fsurdat.zero_nonveg()
        logger.info("zero_nonveg complete")

    # set_dom_pft follows idealized and zero_nonveg because it modifies
    # PCT_NATVEG and PCT_CROP in the user-defined rectangle
    if dom_pft is not None:
        modify_fsurdat.set_dom_pft(
            dom_pft=dom_pft, lai=lai, sai=sai, hgt_top=hgt_top, hgt_bot=hgt_bot
        )
        logger.info("dom_pft complete")

    if evenly_split_cropland:
        modify_fsurdat.evenly_split_cropland()
        logger.info("evenly_split_cropland complete")


def read_cfg_optional_basic_opts(modify_fsurdat, config, cfg_path, section):
    """Read the optional parts of the main section of the config file.
    The main section is called modify_fsurdat_basic_options.
    Users may set these optional parts but are not required to do so."""

    lai = get_config_value(
        config=config,
        section=section,
        item="lai",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    sai = get_config_value(
        config=config,
        section=section,
        item="sai",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    hgt_top = get_config_value(
        config=config,
        section=section,
        item="hgt_top",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    hgt_bot = get_config_value(
        config=config,
        section=section,
        item="hgt_bot",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )

    max_soil_color = int(modify_fsurdat.file.mxsoil_color)
    soil_color = get_config_value(
        config=config,
        section=section,
        item="soil_color",
        file_path=cfg_path,
        allowed_values=range(1, max_soil_color + 1),  # 1 to max_soil_color
        convert_to_type=int,
        can_be_unset=True,
    )

    std_elev = get_config_value(
        config=config,
        section=section,
        item="std_elev",
        file_path=cfg_path,
        convert_to_type=float,
        can_be_unset=True,
    )
    max_sat_area = get_config_value(
        config=config,
        section=section,
        item="max_sat_area",
        file_path=cfg_path,
        convert_to_type=float,
        can_be_unset=True,
    )
    return (
        max_sat_area,
        std_elev,
        soil_color,
        lai,
        sai,
        hgt_top,
        hgt_bot,
    )


def read_cfg_option_control(
    modify_fsurdat,
    config,
    section,
    cfg_path,
):
    """Read the option control section"""
    # required but fallback values available for variables omitted
    # entirely from the .cfg file
    idealized = get_config_value(
        config=config,
        section=section,
        item="idealized",
        file_path=cfg_path,
        convert_to_type=bool,
    )
    if idealized:
        logger.info("idealized option is on")
    else:
        logger.info("idealized option is off")
    process_subgrid = get_config_value(
        config=config,
        section=section,
        item="process_subgrid_section",
        file_path=cfg_path,
        convert_to_type=bool,
    )
    if process_subgrid:
        logger.info("process_subgrid_section option is on")
    else:
        logger.info("process_subgrid_section option is off")
    process_var_list = get_config_value(
        config=config,
        section=section,
        item="process_var_list_section",
        file_path=cfg_path,
        convert_to_type=bool,
    )
    if process_var_list:
        logger.info("process_var_list_section option is on")
    else:
        logger.info("process_var_list_section option is off")
    include_nonveg = get_config_value(
        config=config,
        section=section,
        item="include_nonveg",
        file_path=cfg_path,
        convert_to_type=bool,
    )
    if include_nonveg:
        logger.info("include_nonveg option is on")
    else:
        logger.info("include_nonveg option is off")
    max_pft = int(max(modify_fsurdat.file.lsmpft))
    dom_pft = get_config_value(
        config=config,
        section=section,
        item="dom_pft",
        file_path=cfg_path,
        allowed_values=range(max_pft + 1),  # integers from 0 to max_pft
        convert_to_type=int,
        can_be_unset=True,
    )
    if dom_pft:
        logger.info("dom_pft option is on and = %s", str(dom_pft))
    else:
        logger.info("dom_pft option is off")
    evenly_split_cropland = get_config_value(
        config=config,
        section=section,
        item="evenly_split_cropland",
        file_path=cfg_path,
        convert_to_type=bool,
    )
    if (
        evenly_split_cropland
        and dom_pft is not None
        and dom_pft > int(max(modify_fsurdat.file.natpft.values))
    ):
        abort("dom_pft must not be set to a crop PFT when evenly_split_cropland is True")
    if process_subgrid and idealized:
        abort("idealized AND process_subgrid_section can NOT both be on, pick one or the other")

    return (
        idealized,
        process_subgrid,
        process_var_list,
        include_nonveg,
        dom_pft,
        evenly_split_cropland,
    )


def read_cfg_required_basic_opts(config, section, cfg_path):
    """Read the required part of the control section"""
    lnd_lat_1 = get_config_value(
        config=config,
        section=section,
        item="lnd_lat_1",
        file_path=cfg_path,
        convert_to_type=float,
    )
    lnd_lat_2 = get_config_value(
        config=config,
        section=section,
        item="lnd_lat_2",
        file_path=cfg_path,
        convert_to_type=float,
    )
    lnd_lon_1 = get_config_value(
        config=config,
        section=section,
        item="lnd_lon_1",
        file_path=cfg_path,
        convert_to_type=float,
    )
    lnd_lon_2 = get_config_value(
        config=config,
        section=section,
        item="lnd_lon_2",
        file_path=cfg_path,
        convert_to_type=float,
    )
    lon_type = get_config_value(
        config=config,
        section=section,
        item="lon_type",
        file_path=cfg_path,
        convert_to_type=int,
    )

    landmask_file = get_config_value(
        config=config,
        section=section,
        item="landmask_file",
        file_path=cfg_path,
        can_be_unset=True,
    )

    lat_dimname = get_config_value(
        config=config, section=section, item="lat_dimname", file_path=cfg_path, can_be_unset=True
    )
    lon_dimname = get_config_value(
        config=config, section=section, item="lon_dimname", file_path=cfg_path, can_be_unset=True
    )
    return (
        lnd_lat_1,
        lnd_lat_2,
        lnd_lon_1,
        lnd_lon_2,
        landmask_file,
        lat_dimname,
        lon_dimname,
        lon_type,
    )


def fsurdat_modifier(parser):
    """Implementation of fsurdat_modifier command"""
    # read the .cfg (config) file
    cfg_path = str(parser.cfg_path)
    config = ConfigParser()
    config.read(cfg_path)
    section = "modify_fsurdat_basic_options"
    if not config.has_section(section):
        abort("Config file does not have the expected section: " + section)

    if parser.fsurdat_in == "UNSET":
        # required: user must set these in the .cfg file
        fsurdat_in = get_config_value(
            config=config, section=section, item="fsurdat_in", file_path=cfg_path
        )
    else:
        if config.has_option(section=section, option="fsurdat_in"):
            abort("fsurdat_in is specified in both the command line and the config file, pick one")
        fsurdat_in = str(parser.fsurdat_in)

    # Error checking of input file
    if not os.path.exists(fsurdat_in):
        abort("Input fsurdat_in file does NOT exist: " + str(fsurdat_in))

    if parser.fsurdat_out == "UNSET":
        fsurdat_out = get_config_value(
            config=config, section=section, item="fsurdat_out", file_path=cfg_path
        )
    else:
        if config.has_option(section=section, option="fsurdat_out"):
            abort("fsurdat_out is specified in both the command line and the config file, pick one")
        fsurdat_out = str(parser.fsurdat_out)

    # If output file exists, abort before starting work
    if os.path.exists(fsurdat_out):
        if not parser.overwrite:
            errmsg = "Output file already exists: " + fsurdat_out
            abort(errmsg)
        else:
            warnmsg = (
                "Output file already exists"
                + ", but the overwrite option was selected so the file will be overwritten."
            )
            logger.warning(warnmsg)
    (
        lnd_lat_1,
        lnd_lat_2,
        lnd_lon_1,
        lnd_lon_2,
        landmask_file,
        lat_dimname,
        lon_dimname,
        lon_type,
    ) = read_cfg_required_basic_opts(config, section, cfg_path)
    # Create ModifyFsurdat object
    modify_fsurdat = ModifyFsurdat.init_from_file(
        fsurdat_in=fsurdat_in,
        lon_1=lnd_lon_1,
        lon_2=lnd_lon_2,
        lat_1=lnd_lat_1,
        lat_2=lnd_lat_2,
        landmask_file=landmask_file,
        lat_dimname=lat_dimname,
        lon_dimname=lon_dimname,
        lon_type=lon_type,
    )

    # Read control information about the optional sections
    (
        idealized,
        process_subgrid,
        process_var_list,
        include_nonveg,
        dom_pft,
        evenly_split_cropland,
    ) = read_cfg_option_control(
        modify_fsurdat,
        config,
        section,
        cfg_path,
    )

    # Read parts that are optional
    (
        max_sat_area,
        std_elev,
        soil_color,
        lai,
        sai,
        hgt_top,
        hgt_bot,
    ) = read_cfg_optional_basic_opts(modify_fsurdat, config, cfg_path, section)
    # ------------------------------
    # modify surface data properties
    # ------------------------------

    modify_optional(
        modify_fsurdat=modify_fsurdat,
        idealized=idealized,
        include_nonveg=include_nonveg,
        max_sat_area=max_sat_area,
        std_elev=std_elev,
        soil_color=soil_color,
        dom_pft=dom_pft,
        evenly_split_cropland=evenly_split_cropland,
        lai=lai,
        sai=sai,
        hgt_top=hgt_top,
        hgt_bot=hgt_bot,
    )
    #
    # Handle optional sections
    #
    if process_subgrid:
        subgrid = read_cfg_subgrid(config, cfg_path, numurbl=modify_fsurdat.get_urb_dens())
        modify_fsurdat.set_varlist(subgrid, cfg_path)
        logger.info("process_subgrid is complete")
    else:
        check_no_subgrid_section(config)

    if process_var_list:
        varlist = read_cfg_var_list(config, idealized=idealized)
        update_list = modify_fsurdat.check_varlist(
            varlist, allow_uppercase_vars=True, source="Config file: " + cfg_path
        )
        modify_fsurdat.set_varlist(update_list, cfg_path)
        logger.info("process_var_list is complete")
    else:
        check_no_varlist_section(config)

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    write_output(modify_fsurdat.file, fsurdat_in, fsurdat_out, "fsurdat")
