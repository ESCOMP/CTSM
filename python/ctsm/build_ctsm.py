"""Functions implementing build_ctsm command"""

import argparse
import logging
import os
import string
import subprocess
import sys

from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args
from ctsm.os_utils import run_cmd_output_on_error
from ctsm.path_utils import path_to_ctsm_root

logger = logging.getLogger(__name__)

# ========================================================================
# Define some constants
# ========================================================================

# this matches the machine name in config_machines_template.xml
_MACH_NAME = 'ctsm_build'

# these are arbitrary, since we only use the case for its build, not any of the runtime
# settings; they just need to be valid
_COMPSET = 'I2000Ctsm50NwpSpNldasRsGs'
_RES = 'nldas2_rnldas2_mnldas2'

_MACHINE_CONFIG_DIRNAME = 'machine_configuration'
_INPUTDATA_DIRNAME = 'inputdata'
_GPTL_NANOTIMERS_CPPDEFS = '-DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY' # pylint: disable=line-too-long

# ========================================================================
# Public functions
# ========================================================================

def main(cime_path):
    """Main function called when build_ctsm is run from the command-line

    Args:
    cime_path (str): path to the cime that we're using (this is passed in explicitly
        rather than relying on calling path_to_cime so that we can be absolutely sure that
        the scripts called here are coming from the same cime as the cime library we're
        using).
    """
    setup_logging_pre_config()
    args = _commandline_args()
    process_logging_args(args)

    if args.rebuild:
        sys.exit('ERROR: --rebuild not yet implemented')
    else:
        build_ctsm(cime_path=cime_path,
                   build_dir=args.build_dir,
                   os_type=args.os,
                   compiler=args.compiler,
                   netcdf_path=args.netcdf_path,
                   esmf_lib_path=args.esmf_lib_path,
                   gmake=args.gmake,
                   gmake_j=args.gmake_j,
                   pnetcdf_path=args.pnetcdf_path,
                   pio_filesystem_hints=args.pio_filesystem_hints,
                   gptl_nano_timers=args.gptl_nano_timers,
                   extra_fflags=args.extra_fflags,
                   extra_cflags=args.extra_cflags)

def build_ctsm(cime_path,
               build_dir,
               os_type,
               compiler,
               netcdf_path,
               esmf_lib_path,
               gmake,
               gmake_j,
               pnetcdf_path=None,
               pio_filesystem_hints=None,
               gptl_nano_timers=False,
               extra_fflags='',
               extra_cflags=''):
    """Implementation of build_ctsm command

    Args:
    cime_path (str): path to root of cime
    build_dir (str): path to build directory
    os_type (str): operating system type; one of linux, aix, darwin or cnl
    compiler (str): compiler type
    netcdf_path (str): path to NetCDF installation
    esmf_lib_path (str): path to ESMF library directory
    gmake (str): name of GNU make tool
    gmake_j (int): number of threads to use when building
    pnetcdf_path (str): path to PNetCDF installation, if present (or None)
    pio_filesystem_hints (str): if present (not None), enable filesystem hints for the
        given filesystem type
    gptl_nano_timers (bool): if True, enable timers in build of the GPTL timing library
    extra_fflags (str): any extra flags to include when compiling Fortran files
    extra_cflags (str): any extra flags to include when compiling C files
    """

    os_type = _check_and_transform_os(os_type)
    _create_build_dir(build_dir)
    _fill_out_machine_files(build_dir=build_dir,
                            os_type=os_type,
                            compiler=compiler,
                            netcdf_path=netcdf_path,
                            esmf_lib_path=esmf_lib_path,
                            gmake=gmake,
                            gmake_j=gmake_j,
                            pnetcdf_path=pnetcdf_path,
                            pio_filesystem_hints=pio_filesystem_hints,
                            gptl_nano_timers=gptl_nano_timers,
                            extra_fflags=extra_fflags,
                            extra_cflags=extra_cflags)
    _create_and_build_case(cime_path=cime_path,
                           build_dir=build_dir)

# ========================================================================
# Private functions
# ========================================================================

def _commandline_args(args_to_parse=None):
    """Parse and return command-line arguments

    Args:
    args_to_parse: list of strings or None: Generally only used for unit testing; if None,
        reads args from sys.argv
    """

    description = """
Script to build CTSM library and its dependencies

Typical usage:

    For a fresh build:

        build_ctsm /path/to/nonexistent/directory --os OS --compiler COMPILER --netcdf-path NETCDF_PATH --esmf-lib-path ESMF_LIB_PATH

        (Other optional arguments are also allowed in this usage.)

    For rebuilding:

        build_ctsm /path/to/existing/directory --rebuild

        (No other arguments are allowed in this usage.)
"""

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('build_dir',
                        help='Path to build directory\n'
                        'If --rebuild is given, this should be the path to an existing build,\n'
                        'otherwise this directory must not already exist.')

    parser.add_argument('--rebuild', action='store_true',
                        help='Rebuild in an existing build directory\n'
                        'If given, none of the build-related optional arguments should be given.\n')

    non_rebuild_required = parser.add_argument_group(
        title='required arguments without --rebuild; not allowed with --rebuild')
    non_rebuild_required_list = []

    non_rebuild_required.add_argument('--os', type=str.lower,
                                      choices=['linux', 'aix', 'darwin', 'cnl'],
                                      help='Operating system type')
    non_rebuild_required_list.append('os')

    # For now, only support the compilers that we regularly test with, even though cime
    # supports many other options
    non_rebuild_required.add_argument('--compiler', type=str.lower,
                                      choices=['gnu', 'intel', 'nag', 'pgi'],
                                      help='Compiler type')
    non_rebuild_required_list.append('compiler')

    non_rebuild_required.add_argument('--netcdf-path',
                                      help='Path to NetCDF installation\n'
                                      '(path to top-level directory, containing subdirectories\n'
                                      'named lib, include, etc.)')
    non_rebuild_required_list.append('netcdf-path')

    non_rebuild_required.add_argument('--esmf-lib-path',
                                      help='Path to ESMF library directory\n'
                                      'This directory should include an esmf.mk file')
    non_rebuild_required_list.append('esmf-lib-path')

    non_rebuild_optional = parser.add_argument_group(
        title='optional arguments without --rebuild; not allowed with --rebuild')
    non_rebuild_optional_list = []

    non_rebuild_optional.add_argument('--gmake', default='gmake',
                                      help='Name of GNU Make tool on your system\n'
                                      'Default: gmake')
    non_rebuild_optional_list.append('gmake')

    non_rebuild_optional.add_argument('--gmake-j', default=8, type=int,
                                      help='Number of threads to use when building\n'
                                      'Default: 8')
    non_rebuild_optional_list.append('gmake-j')

    non_rebuild_optional.add_argument('--pnetcdf-path',
                                      help='Path to PNetCDF installation, if present\n')
    non_rebuild_optional_list.append('pnetcdf-path')

    non_rebuild_optional.add_argument('--pio-filesystem-hints', type=str.lower,
                                      choices=['gpfs', 'lustre'],
                                      help='Enable filesystem hints for the given filesystem type\n'
                                      'when building the Parallel IO library')
    non_rebuild_optional_list.append('pio-filesystem-hints')

    non_rebuild_optional.add_argument('--gptl-nano-timers', action='store_true',
                                      help='Enable nano timers in build of the GPTL timing library')
    non_rebuild_optional_list.append('gptl-nano-timers')

    non_rebuild_optional.add_argument('--extra-fflags', default='',
                                      help='Any extra, non-standard flags to include\n'
                                      'when compiling Fortran files\n'
                                      'Tip: to allow a dash at the start of these flags,\n'
                                      'use a quoted string with an initial space, as in:\n'
                                      '    --extra-fflags " -flag1 -flag2"')
    non_rebuild_optional_list.append('extra-fflags')

    non_rebuild_optional.add_argument('--extra-cflags', default='',
                                      help='Any extra, non-standard flags to include\n'
                                      'when compiling C files\n'
                                      'Tip: to allow a dash at the start of these flags,\n'
                                      'use a quoted string with an initial space, as in:\n'
                                      '    --extra-cflags " -flag1 -flag2"')
    non_rebuild_optional_list.append('extra-cflags')

    add_logging_args(parser)

    args = parser.parse_args(args_to_parse)
    if args.rebuild:
        _check_args_rebuild(parser, args, non_rebuild_required_list+non_rebuild_optional_list)
    else:
        _check_args_non_rebuild(parser, args, non_rebuild_required_list)

    return args

def _check_args_rebuild(parser, args, args_not_allowed_in_rebuild):
    """Checks if any arguments not allowed with --rebuild are set

    Calls parser.error if there are problems

    Args:
    parser: ArgumentParser
    args: list of parsed arguments
    args_not_allowed_in_rebuild: list of strings - argument names in this category
    """
    for arg in args_not_allowed_in_rebuild:
        arg_no_dashes = arg.replace('-', '_')
        # To determine whether the user specified an argument, we look at whether it's
        # value differs from its default value. This won't catch the case where the user
        # explicitly set an argument to its default value, but it's not a big deal if we
        # miss printing an error in that case.
        if vars(args)[arg_no_dashes] != parser.get_default(arg_no_dashes):
            parser.error('--{} cannot be provided if --rebuild is set'.format(arg))

def _check_args_non_rebuild(parser, args, non_rebuild_required_list):
    """Checks if any arguments required without --rebuild are absent

    Calls parser.error if there are problems

    Args:
    parser: ArgumentParser
    args: list of parsed arguments
    non_rebuild_required_list: list of strings - argument names in this category
    """
    for arg in non_rebuild_required_list:
        arg_no_dashes = arg.replace('-', '_')
        if vars(args)[arg_no_dashes] is None:
            parser.error('--{} must be provided if --rebuild is not set'.format(arg))

def _check_and_transform_os(os_type):
    """Check validity of os_type argument and transform it to proper case

    os_type should be a lowercase string; returns a transformed string
    """
    transforms = {'linux': 'LINUX',
                  'aix': 'AIX',
                  'darwin': 'Darwin',
                  'cnl': 'CNL'}
    try:
        os_type_transformed = transforms[os_type]
    except KeyError:
        raise ValueError("Unknown OS: {}".format(os_type))
    return os_type_transformed

def _create_build_dir(build_dir):
    """Create the given build directory and any necessary sub-directories

    Args:
    build_dir (str): path to build directory; this directory shouldn't exist yet!
    """
    if os.path.exists(build_dir):
        sys.exit('ERROR: When running without --rebuild, the build directory must not exist yet\n'
                 '(<{}> already exists)'.format(build_dir))
    os.makedirs(build_dir)
    os.makedirs(os.path.join(build_dir, _INPUTDATA_DIRNAME))
    os.makedirs(os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME))

def _fill_out_machine_files(build_dir,
                            os_type,
                            compiler,
                            netcdf_path,
                            esmf_lib_path,
                            gmake,
                            gmake_j,
                            pnetcdf_path=None,
                            pio_filesystem_hints=None,
                            gptl_nano_timers=False,
                            extra_fflags='',
                            extra_cflags=''):
    """Fill out the machine porting templates for this machine / compiler

    For documentation of args, see the documentation in the build_ctsm function
    """
    path_to_templates = os.path.join(path_to_ctsm_root(),
                                     'lilac_config',
                                     'build_templates')

    # ------------------------------------------------------------------------
    # Fill in config_machines.xml
    # ------------------------------------------------------------------------

    with open(os.path.join(path_to_templates, 'config_machines_template.xml')) as cm_template_file:
        cm_template_file_contents = cm_template_file.read()
    config_machines_template = string.Template(cm_template_file_contents)
    config_machines = config_machines_template.substitute(
        OS=os_type,
        COMPILER=compiler,
        CIME_OUTPUT_ROOT=build_dir,
        GMAKE=gmake,
        GMAKE_J=gmake_j)
    with open(os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME, 'config_machines.xml'),
              'w') as cm_file:
        cm_file.write(config_machines)

    # ------------------------------------------------------------------------
    # Fill in config_compilers.xml
    # ------------------------------------------------------------------------

    if gptl_nano_timers:
        gptl_cppdefs = _GPTL_NANOTIMERS_CPPDEFS
    else:
        gptl_cppdefs = ''

    if pio_filesystem_hints:
        pio_filesystem_hints_tag = '<PIO_FILESYSTEM_HINTS>{}</PIO_FILESYSTEM_HINTS>'.format(
            pio_filesystem_hints)
    else:
        pio_filesystem_hints_tag = ''

    if pnetcdf_path:
        pnetcdf_path_tag = '<PNETCDF_PATH>{}</PNETCDF_PATH>'.format(
            pnetcdf_path)
    else:
        pnetcdf_path_tag = ''

    with open(os.path.join(path_to_templates, 'config_compilers_template.xml')) as cc_template_file:
        cc_template_file_contents = cc_template_file.read()
    config_compilers_template = string.Template(cc_template_file_contents)
    config_compilers = config_compilers_template.substitute(
        COMPILER=compiler,
        GPTL_CPPDEFS=gptl_cppdefs,
        NETCDF_PATH=netcdf_path,
        PIO_FILESYSTEM_HINTS=pio_filesystem_hints_tag,
        PNETCDF_PATH=pnetcdf_path_tag,
        ESMF_LIBDIR=esmf_lib_path,
        EXTRA_CFLAGS=extra_cflags,
        EXTRA_FFLAGS=extra_fflags)
    with open(os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME, 'config_compilers.xml'),
              'w') as cc_file:
        cc_file.write(config_compilers)

def _create_and_build_case(cime_path, build_dir):
    """Create a case and build the CTSM library and its dependencies

    Args:
    cime_path (str): path to root of cime
    build_dir (str): path to build directory
    """
    casedir = os.path.join(build_dir, 'case')

    # Note that, for some commands, we want to suppress output, only showing the output if
    # the command fails; for these we use run_cmd_output_on_error. For other commands, we
    # want to always show output (or there should be no output in general); for these, we
    # directly use subprocess.check_call or similar.

    run_cmd_output_on_error(
        ['create_newcase',
         '--case', casedir,
         '--compset', _COMPSET,
         '--res', _RES,
         '--machine', _MACH_NAME,
         '--driver', 'nuopc',
         '--extra-machines-dir', os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME),
         '--run-unsupported'],
        errmsg='Problem creating CTSM case directory',
        cwd=os.path.join(cime_path, 'scripts'))

    run_cmd_output_on_error(['case.setup'],
                            errmsg='Problem setting up CTSM case directory',
                            cwd=casedir)

    subprocess.check_call(['xmlchange', 'LILAC_MODE=on'], cwd=casedir)

    try:
        subprocess.check_call(
            ['case.build',
             '--sharedlib-only'],
            cwd=casedir)
    except subprocess.CalledProcessError:
        sys.exit('ERROR building CTSM or its dependencies - see above for details')
