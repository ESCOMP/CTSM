"""Functions implementing LILAC's make_runtime_inputs command"""

import os
import subprocess
import argparse
import logging

from configparser import ConfigParser

from CIME.buildnml import create_namelist_infile  # pylint: disable=import-error

from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)
from ctsm.path_utils import path_to_ctsm_root
from ctsm.utils import abort
from ctsm.config_utils import get_config_value

logger = logging.getLogger(__name__)

# ========================================================================
# Define some constants
# ========================================================================

_CONFIG_CACHE_TEMPLATE = """
<?xml version="1.0"?>
<config_definition>
<commandline></commandline>
<entry id="phys" value="{clm_phys}" list="" valid_values="clm4_5,clm5_0,clm6_0">Specifies ctsm physics</entry>
</config_definition>
"""

# Note the following is needed in env_lilac.xml otherwise the following error appears in
# the call to build_namelist

# err=ERROR : CLM build-namelist::CLMBuildNamelist::logical_to_fortran() :
# Unexpected value in logical_to_fortran:

_ENV_LILAC_TEMPLATE = """
<?xml version="1.0"?>
<file id="env_lilac.xml" version="2.0">
  <group id="run_glc">
    <entry id="GLC_TWO_WAY_COUPLING" value="FALSE">
      <type>logical</type>
      <valid_values>TRUE,FALSE</valid_values>
    </entry>
  </group>
  <group id="run_cpl">
    <entry id="LND_SETS_DUST_EMIS_DRV_FLDS" value="TRUE">
      <type>logical</type>
      <valid_values>TRUE,FALSE</valid_values>
    </entry>
  </group>
</file>
"""

# ========================================================================
# Fake case class that can be used to satisfy the interface of CIME functions that need a
# case object
# ========================================================================


class CaseFake:
    """Fake case class to satisfy interface of CIME functions that need a case object"""

    # pylint: disable=too-few-public-methods

    def __init__(self):
        pass

    @staticmethod
    def get_resolved_value(value):
        """Make sure get_resolved_value doesn't get called

        (since we don't have a real case object to resolve values with)
        """
        abort("Cannot resolve value with a '$' variable: {}".format(value))


###############################################################################
def parse_command_line():
    ###############################################################################

    """Parse the command line, return object holding arguments"""

    description = """
Script to create runtime inputs when running CTSM via LILAC
"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, description=description
    )

    parser.add_argument(
        "--rundir",
        type=str,
        default=os.getcwd(),
        help="Full path of the run directory (containing ctsm.cfg & user_nl_ctsm)",
    )

    add_logging_args(parser)

    arguments = parser.parse_args()

    # Perform some error checking on arguments

    if not os.path.isdir(arguments.rundir):
        abort("rundir {} does not exist".format(arguments.rundir))

    return arguments


###############################################################################
def determine_bldnml_opts(bgc_mode, crop, vichydro):
    ###############################################################################
    """Return a string giving bldnml options, given some other inputs"""
    bldnml_opts = ""
    bldnml_opts += " -bgc {}".format(bgc_mode)
    if bgc_mode == "fates":
        # BUG(wjs, 2020-06-12, ESCOMP/CTSM#115) For now, FATES is incompatible with MEGAN
        bldnml_opts += " -no-megan"

    if crop == "on":
        if bgc_mode not in ["bgc", "cn"]:
            abort("Error: setting crop to 'on' is only compatible with bgc_mode of 'bgc' or 'cn'")
        bldnml_opts += " -crop"

    if vichydro == "on":
        if bgc_mode != "sp":
            abort("Error: setting vichydro to 'on' is only compatible with bgc_mode of 'sp'")
        bldnml_opts += " -vichydro"

    return bldnml_opts


###############################################################################
def buildnml(cime_path, rundir):
    ###############################################################################

    """Build the ctsm namelist"""

    # pylint: disable=too-many-locals
    # pylint: disable=too-many-statements

    ctsm_cfg_path = os.path.join(rundir, "ctsm.cfg")

    # read the config file
    config = ConfigParser()
    config.read(ctsm_cfg_path)

    lnd_domain_file = get_config_value(config, "buildnml_input", "lnd_domain_file", ctsm_cfg_path)
    fsurdat = get_config_value(
        config, "buildnml_input", "fsurdat", ctsm_cfg_path, can_be_unset=True
    )
    finidat = get_config_value(
        config, "buildnml_input", "finidat", ctsm_cfg_path, can_be_unset=True
    )

    ctsm_phys = get_config_value(
        config,
        "buildnml_input",
        "ctsm_phys",
        ctsm_cfg_path,
        allowed_values=["clm4_5", "clm5_0", "clm6_0"],
    )
    configuration = get_config_value(
        config,
        "buildnml_input",
        "configuration",
        ctsm_cfg_path,
        allowed_values=["nwp", "clm"],
    )
    structure = get_config_value(
        config,
        "buildnml_input",
        "structure",
        ctsm_cfg_path,
        allowed_values=["fast", "standard"],
    )
    bgc_mode = get_config_value(
        config,
        "buildnml_input",
        "bgc_mode",
        ctsm_cfg_path,
        allowed_values=["sp", "bgc", "cn", "fates"],
    )
    crop = get_config_value(
        config, "buildnml_input", "crop", ctsm_cfg_path, allowed_values=["off", "on"]
    )
    vichydro = get_config_value(
        config,
        "buildnml_input",
        "vichydro",
        ctsm_cfg_path,
        allowed_values=["off", "on"],
    )

    bldnml_opts = determine_bldnml_opts(bgc_mode=bgc_mode, crop=crop, vichydro=vichydro)

    co2_ppmv = get_config_value(config, "buildnml_input", "co2_ppmv", ctsm_cfg_path)
    use_case = get_config_value(config, "buildnml_input", "use_case", ctsm_cfg_path)
    lnd_tuning_mode = get_config_value(config, "buildnml_input", "lnd_tuning_mode", ctsm_cfg_path)
    spinup = get_config_value(
        config, "buildnml_input", "spinup", ctsm_cfg_path, allowed_values=["off", "on"]
    )

    inputdata_path = get_config_value(config, "buildnml_input", "inputdata_path", ctsm_cfg_path)

    # Parse the user_nl_ctsm file
    infile = os.path.join(rundir, ".namelist")
    create_namelist_infile(
        case=CaseFake(),
        user_nl_file=os.path.join(rundir, "user_nl_ctsm"),
        namelist_infile=infile,
    )

    # create config_cache.xml file
    # Note that build-namelist utilizes the contents of the config_cache.xml file in
    # the namelist_defaults.xml file to obtain namelist variables
    config_cache = os.path.join(rundir, "config_cache.xml")
    config_cache_text = _CONFIG_CACHE_TEMPLATE.format(clm_phys=ctsm_phys)
    with open(config_cache, "w") as tempfile:
        tempfile.write(config_cache_text)

    # create temporary env_lilac.xml
    env_lilac = os.path.join(rundir, "env_lilac.xml")
    env_lilac_text = _ENV_LILAC_TEMPLATE.format()
    with open(env_lilac, "w") as tempfile:
        tempfile.write(env_lilac_text)

    # remove any existing clm.input_data_list file
    inputdatalist_path = os.path.join(rundir, "ctsm.input_data_list")
    if os.path.exists(inputdatalist_path):
        os.remove(inputdatalist_path)

    # determine if fsurdat and/or finidat should appear in the -namelist option
    extra_namelist_opts = ""
    if fsurdat is not None:
        # NOTE(wjs, 2020-06-30) With the current logic, fsurdat should never be UNSET
        # (ie None here) but it's possible that this will change in the future.
        extra_namelist_opts = extra_namelist_opts + " fsurdat = '{}' ".format(fsurdat)
    if finidat is not None:
        extra_namelist_opts = extra_namelist_opts + " finidat = '{}' ".format(finidat)

    # call build-namelist
    cmd = os.path.abspath(os.path.join(path_to_ctsm_root(), "bld", "build-namelist"))
    command = [
        cmd,
        "-driver",
        "nuopc",
        "-cimeroot",
        cime_path,
        "-infile",
        infile,
        "-csmdata",
        inputdata_path,
        "-inputdata",
        inputdatalist_path,
        # Hard-code start_ymd of year-2000. This is used to set the run type (for
        # which a setting of 2000 gives 'startup', which is what we want) and pick
        # the initial conditions file (which is pretty much irrelevant when running
        # with lilac).
        "-namelist",
        "&clm_inparm  start_ymd=20000101 {} /".format(extra_namelist_opts),
        "-use_case",
        use_case,
        # For now, we assume ignore_ic_year, not ignore_ic_date
        "-ignore_ic_year",
        # -clm_start_type seems unimportant (see discussion in
        # https://github.com/ESCOMP/CTSM/issues/876)
        "-clm_start_type",
        "default",
        "-configuration",
        configuration,
        "-structure",
        structure,
        "-lilac",
        "-lnd_frac",
        lnd_domain_file,
        "-glc_nec",
        str(10),
        "-co2_ppmv",
        co2_ppmv,
        "-co2_type",
        "constant",
        "-clm_accelerated_spinup",
        spinup,
        "-lnd_tuning_mode",
        lnd_tuning_mode,
        # Eventually make -no-megan dynamic (see
        # https://github.com/ESCOMP/CTSM/issues/926)
        "-no-megan",
        "-config",
        os.path.join(rundir, "config_cache.xml"),
        "-envxml_dir",
        rundir,
    ]
    # NOTE(wjs, 2020-06-16) Note that we do NOT use the -mask argument; it's possible that
    # we should be using it in some circumstances (I haven't looked into how it's used).
    command.extend(["-res", "lilac", "-clm_usr_name", "lilac"])
    command.extend(bldnml_opts.split())

    subprocess.check_call(command, universal_newlines=True)

    # remove temporary files in rundir
    os.remove(os.path.join(rundir, "config_cache.xml"))
    os.remove(os.path.join(rundir, "env_lilac.xml"))
    os.remove(infile)


###############################################################################
def main(cime_path):
    """Main function

    Args:
    cime_path (str): path to the cime that we're using (this is passed in explicitly
        rather than relying on calling path_to_cime so that we can be absolutely sure that
        the scripts called here are coming from the same cime as the cime library we're
        using).
    """
    setup_logging_pre_config()
    args = parse_command_line()
    process_logging_args(args)

    buildnml(cime_path=cime_path, rundir=args.rundir)


###############################################################################
