"""Functions implementing build_ctsm command"""

import argparse
import logging
import os
import subprocess

from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args
from ctsm.os_utils import run_cmd_output_on_error, make_link
from ctsm.path_utils import path_to_ctsm_root
from ctsm.utils import abort, fill_template_file

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

_PATH_TO_TEMPLATES = os.path.join(path_to_ctsm_root(),
                                  'lilac',
                                  'bld_templates')

_MACHINE_CONFIG_DIRNAME = 'machine_configuration'
_INPUTDATA_DIRNAME = 'inputdata'
_RUNTIME_INPUTS_DIRNAME = 'runtime_inputs'

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
    build_dir = os.path.abspath(args.build_dir)

    if args.rebuild:
        rebuild_ctsm(build_dir=build_dir)
    else:
        build_ctsm(cime_path=cime_path,
                   build_dir=build_dir,
                   compiler=args.compiler,
                   no_build=args.no_build,
                   machine=args.machine,
                   os_type=args.os,
                   netcdf_path=args.netcdf_path,
                   esmf_lib_path=args.esmf_lib_path,
                   gmake=args.gmake,
                   gmake_j=args.gmake_j,
                   pnetcdf_path=args.pnetcdf_path,
                   pio_filesystem_hints=args.pio_filesystem_hints,
                   gptl_nano_timers=args.gptl_nano_timers,
                   extra_fflags=args.extra_fflags,
                   extra_cflags=args.extra_cflags,
                   build_debug=args.build_debug,
                   build_without_openmp=args.build_without_openmp)

def build_ctsm(cime_path,
               build_dir,
               compiler,
               no_build=False,
               machine=None,
               os_type=None,
               netcdf_path=None,
               esmf_lib_path=None,
               gmake=None,
               gmake_j=None,
               pnetcdf_path=None,
               pio_filesystem_hints=None,
               gptl_nano_timers=False,
               extra_fflags='',
               extra_cflags='',
               build_debug=False,
               build_without_openmp=False):
    """Implementation of build_ctsm command

    Args:
    cime_path (str): path to root of cime
    build_dir (str): path to build directory
    compiler (str): compiler type
    no_build (bool): If True, set things up, but skip doing the actual build
    machine (str or None): machine name (a machine known to cime)
    os_type (str or None): operating system type; one of linux, aix, darwin or cnl
        Must be given if machine isn't given; ignored if machine is given
    netcdf_path (str or None): path to NetCDF installation
        Must be given if machine isn't given; ignored if machine is given
    esmf_lib_path (str or None): path to ESMF library directory
        Must be given if machine isn't given; ignored if machine is given
    gmake (str or None): name of GNU make tool
        Must be given if machine isn't given; ignored if machine is given
    gmake_j (int or None): number of threads to use when building
        Must be given if machine isn't given; ignored if machine is given
    pnetcdf_path (str or None): path to PNetCDF installation, if present (or None)
        Ignored if machine is given
    pio_filesystem_hints (str or None): if present (not None), enable filesystem hints for the
        given filesystem type
        Ignored if machine is given
    gptl_nano_timers (bool): if True, enable timers in build of the GPTL timing library
        Ignored if machine is given
    extra_fflags (str): any extra flags to include when compiling Fortran files
        Ignored if machine is given
    extra_cflags (str): any extra flags to include when compiling C files
        Ignored if machine is given
    build_debug (bool): if True, build with flags for debugging
    build_without_openmp (bool): if True, build without OpenMP support
    """

    existing_machine = machine is not None
    _create_build_dir(build_dir=build_dir,
                      existing_machine=existing_machine)

    if machine is None:
        assert os_type is not None, 'with machine absent, os_type must be given'
        assert netcdf_path is not None, 'with machine absent, netcdf_path must be given'
        assert esmf_lib_path is not None, 'with machine absent, esmf_lib_path must be given'
        os_type = _check_and_transform_os(os_type)
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

    _create_case(cime_path=cime_path,
                 build_dir=build_dir,
                 compiler=compiler,
                 machine=machine,
                 build_debug=build_debug,
                 build_without_openmp=build_without_openmp)

    if existing_machine:
        # For a user-defined machine, we create an inputdata directory for this case. For
        # an existing cime-ported machine, we still want an inputdata directory alongside
        # the other directories, but now it will just be a link to the real inputdata
        # space on that machine. (Note that, for a user-defined machine, it's important
        # that we have created this directory before creating the case, whereas for an
        # existing machine, we need to wait until after we have created the case to know
        # where to make the sym link point to.)
        _link_to_inputdata(build_dir=build_dir)

    _stage_runtime_inputs(build_dir=build_dir)
    if not no_build:
        _build_case(build_dir=build_dir)

def rebuild_ctsm(build_dir):
    """Re-run the build in an existing directory

    Args:
    build_dir (str): path to build directory
    """
    if not os.path.exists(build_dir):
        abort('When running with --rebuild, the build directory must already exist\n'
              '(<{}> does not exist)'.format(build_dir))

    case_dir = _get_case_dir(build_dir)
    if not os.path.exists(case_dir):
        abort('It appears there was a problem setting up the initial build in\n'
              '<{}>\n'
              'You should start over with a fresh build directory.'.format(build_dir))

    try:
        subprocess.check_call(
            [os.path.join(case_dir, 'case.build'),
             '--clean-depends',
             'lnd'],
            cwd=case_dir)
    except subprocess.CalledProcessError:
        abort('ERROR resetting build for CTSM in order to rebuild - see above for details')

    _build_case(build_dir)

# ========================================================================
# Private functions
# ========================================================================

def _commandline_args(args_to_parse=None):
    """Parse and return command-line arguments

    Args:
    args_to_parse: list of strings or None: Generally only used for unit testing; if None,
        reads args from sys.argv
    """
    # pylint: disable=line-too-long
    # pylint: disable=too-many-statements

    description = """
Script to build CTSM library and its dependencies

Typical usage:

    For a fresh build with a machine that has been ported to cime
    (http://esmci.github.io/cime/versions/master/html/users_guide/porting-cime.html):

        build_ctsm /path/to/nonexistent/directory --machine MACHINE --compiler COMPILER

        (Some other optional arguments are also allowed in this usage, but many are not.)

    For a fresh build with a machine that has NOT been ported to cime:

        build_ctsm /path/to/nonexistent/directory --compiler COMPILER --os OS --netcdf-path NETCDF_PATH --esmf-lib-path ESMF_LIB_PATH

        (Other optional arguments are also allowed in this usage.)

    For rebuilding:

        build_ctsm /path/to/existing/directory --rebuild

        (Most other arguments are NOT allowed in this usage.)
"""

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('build_dir',
                        help='Path to build directory\n'
                        'If --rebuild is given, this should be the path to an existing build,\n'
                        'otherwise this directory must not already exist.')

    main_opts = parser.add_mutually_exclusive_group()

    main_opts.add_argument('--machine',
                           help='Name of machine; this must be a machine that has been ported to cime\n'
                           '(http://esmci.github.io/cime/versions/master/html/users_guide/porting-cime.html)\n'
                           'If given, then none of the machine-definition optional arguments should be given.\n')

    main_opts.add_argument('--rebuild', action='store_true',
                           help='Rebuild in an existing build directory\n'
                           'If given, none of the machine-definition or build-related optional arguments\n'
                           'should be given.\n')

    non_rebuild_required = parser.add_argument_group(
        title='required arguments when not rebuilding',
        description='These arguments are required if --rebuild is not given; '
        'they are not allowed with --rebuild:')
    non_rebuild_required_list = []

    # For now, only support the compilers that we regularly test with, even though cime
    # supports many other options
    non_rebuild_required.add_argument('--compiler', type=str.lower,
                                      choices=['gnu', 'intel', 'nag', 'pgi'],
                                      help='Compiler type')
    non_rebuild_required_list.append('compiler')

    non_rebuild_optional = parser.add_argument_group(
        title='optional arguments when not rebuilding',
        description='These arguments are optional if --rebuild is not given; '
        'they are not allowed with --rebuild:')
    non_rebuild_optional_list = []

    non_rebuild_optional.add_argument('--build-debug', action='store_true',
                                      help='Build with flags for debugging rather than production runs')
    non_rebuild_optional_list.append('build-debug')

    non_rebuild_optional.add_argument('--build-without-openmp', action='store_true',
                                      help='By default, CTSM is built with support for OpenMP threading;\n'
                                      'if this flag is set, then CTSM is built without this support.\n'
                                      'This is mainly useful if your machine/compiler does not support OpenMP.')
    non_rebuild_optional_list.append('build-without-openmp')

    non_rebuild_optional.add_argument('--no-build', action='store_true',
                                      help='Do the pre-build setup, but do not actually build CTSM\n'
                                      '(This is useful for testing, or for expert use.)')
    non_rebuild_optional_list.append('no-build')

    new_machine_required = parser.add_argument_group(
        title='required arguments for a user-defined machine',
        description='These arguments are required if neither --machine nor --rebuild are given; '
        'they are not allowed with either of those arguments:')
    new_machine_required_list = []

    new_machine_required.add_argument('--os', type=str.lower,
                                      choices=['linux', 'aix', 'darwin', 'cnl'],
                                      help='Operating system type')
    new_machine_required_list.append('os')

    new_machine_required.add_argument('--netcdf-path',
                                      help='Path to NetCDF installation\n'
                                      '(path to top-level directory, containing subdirectories\n'
                                      'named lib, include, etc.)')
    new_machine_required_list.append('netcdf-path')

    new_machine_required.add_argument('--esmf-lib-path',
                                      help='Path to ESMF library directory\n'
                                      'This directory should include an esmf.mk file')
    new_machine_required_list.append('esmf-lib-path')

    new_machine_optional = parser.add_argument_group(
        title='optional arguments for a user-defined machine',
        description='These arguments are optional if neither --machine nor --rebuild are given; '
        'they are not allowed with either of those arguments:')
    new_machine_optional_list = []

    new_machine_optional.add_argument('--gmake', default='gmake',
                                      help='Name of GNU Make tool on your system\n'
                                      'Default: gmake')
    new_machine_optional_list.append('gmake')

    new_machine_optional.add_argument('--gmake-j', default=8, type=int,
                                      help='Number of threads to use when building\n'
                                      'Default: 8')
    new_machine_optional_list.append('gmake-j')

    new_machine_optional.add_argument('--pnetcdf-path',
                                      help='Path to PNetCDF installation, if present\n')
    new_machine_optional_list.append('pnetcdf-path')

    new_machine_optional.add_argument('--pio-filesystem-hints', type=str.lower,
                                      choices=['gpfs', 'lustre'],
                                      help='Enable filesystem hints for the given filesystem type\n'
                                      'when building the Parallel IO library')
    new_machine_optional_list.append('pio-filesystem-hints')

    new_machine_optional.add_argument('--gptl-nano-timers', action='store_true',
                                      help='Enable nano timers in build of the GPTL timing library')
    new_machine_optional_list.append('gptl-nano-timers')

    new_machine_optional.add_argument('--extra-fflags', default='',
                                      help='Any extra, non-standard flags to include\n'
                                      'when compiling Fortran files\n'
                                      'Tip: to allow a dash at the start of these flags,\n'
                                      'use a quoted string with an initial space, as in:\n'
                                      '    --extra-fflags " -flag1 -flag2"')
    new_machine_optional_list.append('extra-fflags')

    new_machine_optional.add_argument('--extra-cflags', default='',
                                      help='Any extra, non-standard flags to include\n'
                                      'when compiling C files\n'
                                      'Tip: to allow a dash at the start of these flags,\n'
                                      'use a quoted string with an initial space, as in:\n'
                                      '    --extra-cflags " -flag1 -flag2"')
    new_machine_optional_list.append('extra-cflags')

    add_logging_args(parser)

    args = parser.parse_args(args_to_parse)
    if args.rebuild:
        _confirm_args_absent(parser, args, "cannot be provided if --rebuild is set",
                             (non_rebuild_required_list + non_rebuild_optional_list +
                              new_machine_required_list + new_machine_optional_list))
    else:
        _confirm_args_present(parser, args, "must be provided if --rebuild is not set",
                              non_rebuild_required_list)
        if args.machine:
            _confirm_args_absent(parser, args, "cannot be provided if --machine is set",
                                 new_machine_required_list + new_machine_optional_list)
        else:
            _confirm_args_present(parser, args, "must be provided if neither --machine nor --rebuild are set",
                                  new_machine_required_list)

    return args

def _confirm_args_absent(parser, args, errmsg, args_not_allowed):
    """Confirms that all args not allowed in this usage are absent

    Calls parser.error if there are problems

    Args:
    parser: ArgumentParser
    args: list of parsed arguments
    errmsg: string - message printed if there is a problem
    args_not_allowed: list of strings - argument names in this category
    """
    for arg in args_not_allowed:
        arg_no_dashes = arg.replace('-', '_')
        # To determine whether the user specified an argument, we look at whether it's
        # value differs from its default value. This won't catch the case where the user
        # explicitly set an argument to its default value, but it's not a big deal if we
        # miss printing an error in that case.
        if vars(args)[arg_no_dashes] != parser.get_default(arg_no_dashes):
            parser.error('--{} {}'.format(arg, errmsg))

def _confirm_args_present(parser, args, errmsg, args_required):
    """Confirms that all args required in this usage are present

    Calls parser.error if there are problems

    Args:
    parser: ArgumentParser
    args: list of parsed arguments
    errmsg: string - message printed if there is a problem
    args_required: list of strings - argument names in this category
    """
    for arg in args_required:
        arg_no_dashes = arg.replace('-', '_')
        if vars(args)[arg_no_dashes] is None:
            parser.error('--{} {}'.format(arg, errmsg))

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

def _get_case_dir(build_dir):
    """Given the path to build_dir, return the path to the case directory"""
    return os.path.join(build_dir, 'case')

def _create_build_dir(build_dir, existing_machine):
    """Create the given build directory and any necessary sub-directories

    Args:
    build_dir (str): path to build directory; this directory shouldn't exist yet!
    existing_machine (bool): whether this build is for a machine known to cime
        (as opposed to an on-the-fly machine port)
    """
    if os.path.exists(build_dir):
        abort('When running without --rebuild, the build directory must not exist yet\n'
              '(<{}> already exists)'.format(build_dir))
    os.makedirs(build_dir)
    if not existing_machine:
        os.makedirs(os.path.join(build_dir, _INPUTDATA_DIRNAME))

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
    os.makedirs(os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME))

    # ------------------------------------------------------------------------
    # Fill in config_machines.xml
    # ------------------------------------------------------------------------

    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES, 'config_machines_template.xml'),
        path_to_final=os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME, 'config_machines.xml'),
        substitutions={'OS':os_type,
                       'COMPILER':compiler,
                       'CIME_OUTPUT_ROOT':build_dir,
                       'GMAKE':gmake,
                       'GMAKE_J':gmake_j})

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

    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES,
                                      'config_compilers_template.xml'),
        path_to_final=os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME, 'config_compilers.xml'),
        substitutions={'COMPILER':compiler,
                       'GPTL_CPPDEFS':gptl_cppdefs,
                       'NETCDF_PATH':netcdf_path,
                       'PIO_FILESYSTEM_HINTS':pio_filesystem_hints_tag,
                       'PNETCDF_PATH':pnetcdf_path_tag,
                       'ESMF_LIBDIR':esmf_lib_path,
                       'EXTRA_CFLAGS':extra_cflags,
                       'EXTRA_FFLAGS':extra_fflags})


def _create_case(cime_path, build_dir, compiler,
                 machine=None, build_debug=False, build_without_openmp=False):
    """Create a case that can later be used to build the CTSM library and its dependencies

    Args:
    cime_path (str): path to root of cime
    build_dir (str): path to build directory
    compiler (str): compiler to use
    machine (str or None): name of machine or None
        If None, we assume we're using an on-the-fly machine port
        Otherwise, machine should be the name of a machine known to cime
    build_debug (bool): if True, build with flags for debugging
    build_without_openmp (bool): if True, build without OpenMP support
    """
    # Note that, for some commands, we want to suppress output, only showing the output if
    # the command fails; for these we use run_cmd_output_on_error. For other commands,
    # there should be no output in general; for these, we directly use
    # subprocess.check_call or similar.

    # Also note that, for commands executed from the case directory, we specify the path
    # to the case directory both in the command itself and in the cwd argument. We do the
    # former in case dot isn't in the user's path; we do the latter in case the commands
    # require you to be in the case directory when you execute them.

    case_dir = _get_case_dir(build_dir)
    xmlchange = os.path.join(case_dir, 'xmlchange')

    if machine is None:
        machine_args = ['--machine', _MACH_NAME,
                        '--extra-machines-dir', os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME)]
    else:
        machine_args = ['--machine', machine]

    create_newcase_cmd = [os.path.join(cime_path, 'scripts', 'create_newcase'),
                          '--output-root', build_dir,
                          '--case', case_dir,
                          '--compset', _COMPSET,
                          '--res', _RES,
                          '--compiler', compiler,
                          '--driver', 'nuopc',
                          '--run-unsupported']
    create_newcase_cmd.extend(machine_args)
    run_cmd_output_on_error(create_newcase_cmd,
                            errmsg='Problem creating CTSM case directory')

    run_cmd_output_on_error([os.path.join(case_dir, 'case.setup')],
                            errmsg='Problem setting up CTSM case directory',
                            cwd=case_dir)

    subprocess.check_call([xmlchange, 'LILAC_MODE=on'], cwd=case_dir)
    if build_debug:
        subprocess.check_call([xmlchange, 'DEBUG=TRUE'], cwd=case_dir)
    if not build_without_openmp:
        subprocess.check_call([xmlchange, 'FORCE_BUILD_SMP=TRUE'], cwd=case_dir)

    make_link(os.path.join(case_dir, 'bld'),
              os.path.join(build_dir, 'bld'))
    if machine is not None:
        # For a pre-existing machine, the .env_mach_specific files are likely useful to
        # the user. Make sym links to these with more intuitive names.
        for extension in ('sh', 'csh'):
            make_link(os.path.join(case_dir, '.env_mach_specific.{}'.format(extension)),
                      os.path.join(build_dir, 'ctsm_build_environment.{}'.format(extension)))

def _link_to_inputdata(build_dir):
    """For an existing machine, make a sym link to the inputdata directory

    Args:
    build_dir (str): path to build directory
    """
    case_dir = _get_case_dir(build_dir)
    xmlquery = os.path.join(case_dir, 'xmlquery')

    inputdata_dir = subprocess.check_output([xmlquery, '--value', 'DIN_LOC_ROOT'],
                                            cwd=case_dir,
                                            universal_newlines=True)
    make_link(inputdata_dir,
              os.path.join(build_dir, _INPUTDATA_DIRNAME))

def _stage_runtime_inputs(build_dir):
    """Stage CTSM and LILAC runtime inputs

    Args:
    build_dir (str): path to build directory
    """
    os.makedirs(os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME))

    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES, 'ctsm_template.cfg'),
        path_to_final=os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME, 'ctsm.cfg'),
        substitutions={'INPUTDATA':os.path.join(build_dir, _INPUTDATA_DIRNAME)})

def _build_case(build_dir):
    """Build the CTSM library and its dependencies

    Args:
    build_dir (str): path to build directory
    """
    # We want user to see output from the build command, so we use subprocess.check_call
    # rather than run_cmd_output_on_error.

    # See comment in _create_case for why we use case_dir in both the path to the command
    # and in the cwd argument to check_call.
    case_dir = _get_case_dir(build_dir)
    try:
        subprocess.check_call(
            [os.path.join(case_dir, 'case.build'),
             '--sharedlib-only'],
            cwd=case_dir)
    except subprocess.CalledProcessError:
        abort('ERROR building CTSM or its dependencies - see above for details')

    make_link(os.path.join(case_dir, 'bld', 'ctsm.mk'),
              os.path.join(build_dir, 'ctsm.mk'))
