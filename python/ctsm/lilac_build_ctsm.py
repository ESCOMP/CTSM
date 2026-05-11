"""Functions implementing LILAC's build_ctsm command"""

import argparse
import logging
import os
import shutil
import subprocess

from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)
from ctsm.os_utils import run_cmd_output_on_error, make_link
from ctsm.path_utils import path_to_ctsm_root
from ctsm.utils import abort, fill_template_file

logger = logging.getLogger(__name__)

# ========================================================================
# Define some constants
# ========================================================================

# this matches the machine name in config_machines_template.xml
_MACH_NAME = "ctsm-build"

# these are arbitrary, since we only use the case for its build, not any of the runtime
# settings; they just need to be valid
_COMPSET = "I2000Ctsm50NwpSpAsRs"
_RES = "f10_f10_mg37"

_PATH_TO_TEMPLATES = os.path.join(path_to_ctsm_root(), "lilac", "bld_templates")

_PATH_TO_MAKE_RUNTIME_INPUTS = os.path.join(path_to_ctsm_root(), "lilac", "make_runtime_inputs")

_PATH_TO_DOWNLOAD_INPUT_DATA = os.path.join(path_to_ctsm_root(), "lilac", "download_input_data")

_MACHINE_CONFIG_DIRNAME = "machine_configuration"
_INPUTDATA_DIRNAME = "inputdata"
_RUNTIME_INPUTS_DIRNAME = "runtime_inputs"

_GPTL_NANOTIMERS_CPPDEFS = "-DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY"  # pylint: disable=line-too-long

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
        build_ctsm(
            cime_path=cime_path,
            build_dir=build_dir,
            compiler=args.compiler,
            no_build=args.no_build,
            machine=args.machine,
            os_type=args.os,
            netcdf_path=args.netcdf_path,
            esmf_mkfile_path=args.esmf_mkfile_path,
            max_mpitasks_per_node=args.max_mpitasks_per_node,
            gmake=args.gmake,
            gmake_j=args.gmake_j,
            pnetcdf_path=args.pnetcdf_path,
            pio_filesystem_hints=args.pio_filesystem_hints,
            gptl_nano_timers=args.gptl_nano_timers,
            extra_fflags=args.extra_fflags,
            extra_cflags=args.extra_cflags,
            no_pnetcdf=args.no_pnetcdf,
            build_debug=args.build_debug,
            build_with_openmp=args.build_with_openmp,
            inputdata_path=args.inputdata_path,
        )


def build_ctsm(
    cime_path,
    build_dir,
    compiler,
    *,
    no_build=False,
    machine=None,
    os_type=None,
    netcdf_path=None,
    esmf_mkfile_path=None,
    max_mpitasks_per_node=None,
    gmake=None,
    gmake_j=None,
    pnetcdf_path=None,
    pio_filesystem_hints=None,
    gptl_nano_timers=False,
    extra_fflags="",
    extra_cflags="",
    no_pnetcdf=False,
    build_debug=False,
    build_with_openmp=False,
    inputdata_path=None,
):
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
    esmf_mkfile_path (str or None): path to esmf.mk file (typically within ESMF library directory)
        Must be given if machine isn't given; ignored if machine is given
    max_mpitasks_per_node (int or None): number of physical processors per shared-memory node
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
    no_pnetcdf (bool): if True, use netcdf rather than pnetcdf
    build_debug (bool): if True, build with flags for debugging
    build_with_openmp (bool): if True, build with OpenMP support
    inputdata_path (str or None): path to existing inputdata directory on this machine
        If None, an inputdata directory will be created for this build
        (If machine is given, then we use the machine's inputdata directory by default;
        but if inputdata_path is given, it overrides the machine's inputdata directory.)
    """

    existing_machine = machine is not None
    existing_inputdata = existing_machine or inputdata_path is not None
    _create_build_dir(build_dir=build_dir, existing_inputdata=existing_inputdata)

    # Some error checking
    if inputdata_path is not None:
        if not os.path.isdir(inputdata_path):
            abort("Input inputdata_path directory does NOT exist = " + inputdata_path)

    if not os.path.isdir(build_dir):
        abort("Input build_dir directory does NOT exist = " + build_dir)

    if machine is None:
        assert os_type is not None, "with machine absent, os_type must be given"
        assert netcdf_path is not None, "with machine absent, netcdf_path must be given"
        assert esmf_mkfile_path is not None, "with machine absent, esmf_mkfile_path must be given"
        assert max_mpitasks_per_node is not None, (
            "with machine absent " "max_mpitasks_per_node must be given"
        )
        os_type = _check_and_transform_os(os_type)
        _fill_out_machine_files(
            build_dir=build_dir,
            os_type=os_type,
            compiler=compiler,
            netcdf_path=netcdf_path,
            esmf_mkfile_path=esmf_mkfile_path,
            max_mpitasks_per_node=max_mpitasks_per_node,
            gmake=gmake,
            gmake_j=gmake_j,
            pnetcdf_path=pnetcdf_path,
            pio_filesystem_hints=pio_filesystem_hints,
            gptl_nano_timers=gptl_nano_timers,
            extra_fflags=extra_fflags,
            extra_cflags=extra_cflags,
        )
    assert os.path.isdir(cime_path), "cime_path must be a directory"

    _create_case(
        cime_path=cime_path,
        build_dir=build_dir,
        compiler=compiler,
        machine=machine,
        build_debug=build_debug,
        build_with_openmp=build_with_openmp,
        inputdata_path=inputdata_path,
    )

    _stage_runtime_inputs(build_dir=build_dir, no_pnetcdf=no_pnetcdf)

    print(
        "Initial setup complete; it is now safe to work with the runtime inputs in\n"
        "{}\n".format(os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME))
    )

    if not no_build:
        _build_case(build_dir=build_dir)


def rebuild_ctsm(build_dir):
    """Re-run the build in an existing directory

    Args:
    build_dir (str): path to build directory
    """
    if not os.path.exists(build_dir):
        abort(
            "When running with --rebuild, the build directory must already exist\n"
            "(<{}> does not exist)".format(build_dir)
        )

    case_dir = _get_case_dir(build_dir)
    if not os.path.exists(case_dir):
        abort(
            "It appears there was a problem setting up the initial build in\n"
            "<{}>\n"
            "You should start over with a fresh build directory.".format(build_dir)
        )

    try:
        subprocess.check_call(
            [os.path.join(case_dir, "case.build"), "--clean-depends", "lnd"],
            cwd=case_dir,
        )
    except subprocess.CalledProcessError:
        abort("ERROR resetting build for CTSM in order to rebuild - see above for details")

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

        build_ctsm /path/to/nonexistent/directory --os OS --compiler COMPILER --netcdf-path NETCDF_PATH --esmf-mkfile-path ESMF_MKFILE_PATH --max-mpitasks-per-node MAX_MPITASKS_PER_NODE --pnetcdf-path PNETCDF_PATH

        If PNetCDF is not available, set --no-pnetcdf instead of --pnetcdf-path.

        (Other optional arguments are also allowed in this usage.)

    For rebuilding:

        build_ctsm /path/to/existing/directory --rebuild

        (Most other arguments are NOT allowed in this usage.)
"""

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "build_dir",
        help="Path to build directory\n"
        "If --rebuild is given, this should be the path to an existing build,\n"
        "otherwise this directory must not already exist.",
    )

    main_opts = parser.add_mutually_exclusive_group()

    main_opts.add_argument(
        "--machine",
        help="Name of machine; this must be a machine that has been ported to cime\n"
        "(http://esmci.github.io/cime/versions/master/html/users_guide/porting-cime.html)\n"
        "If given, then none of the machine-definition optional arguments should be given.\n",
    )

    main_opts.add_argument(
        "--rebuild",
        action="store_true",
        help="Rebuild in an existing build directory\n"
        "If given, none of the machine-definition or build-related optional arguments\n"
        "should be given.\n",
    )

    non_rebuild_required = parser.add_argument_group(
        title="required arguments when not rebuilding",
        description="These arguments are required if --rebuild is not given; "
        "they are not allowed with --rebuild:",
    )
    non_rebuild_required_list = []

    # For now, only support the compilers that we regularly test with, even though cime
    # supports many other options
    non_rebuild_required.add_argument(
        "--compiler",
        type=str.lower,
        choices=["gnu", "intel", "nag", "pgi"],
        help="Compiler type",
    )
    non_rebuild_required_list.append("compiler")

    non_rebuild_optional = parser.add_argument_group(
        title="optional arguments when not rebuilding",
        description="These arguments are optional if --rebuild is not given; "
        "they are not allowed with --rebuild:",
    )
    non_rebuild_optional_list = []

    non_rebuild_optional.add_argument(
        "--no-pnetcdf",
        action="store_true",
        help="Use NetCDF instead of PNetCDF for CTSM I/O.\n"
        "On a user-defined machine, you must either set this flag\n"
        "or set --pnetcdf-path. On a cime-ported machine,\n"
        "this flag must be set if PNetCDF is not available\n"
        "for this machine/compiler.",
    )
    non_rebuild_optional_list.append("no-pnetcdf")

    non_rebuild_optional.add_argument(
        "--build-debug",
        action="store_true",
        help="Build with flags for debugging rather than production runs",
    )
    non_rebuild_optional_list.append("build-debug")

    non_rebuild_optional.add_argument(
        "--build-with-openmp",
        action="store_true",
        help="By default, CTSM is built WITHOUT support for OpenMP threading;\n"
        "if this flag is set, then CTSM is built WITH this support.\n"
        "This is important for performance if you will be running with\n"
        "OpenMP threading-based parallelization, or hybrid MPI/OpenMP.",
    )
    non_rebuild_optional_list.append("build-with-openmp")

    non_rebuild_optional.add_argument(
        "--inputdata-path",
        help="Path to directory containing CTSM's NetCDF inputs.\n"
        "For a machine that has been ported to cime, the default is to\n"
        "use this machine's standard inputdata location; this argument\n"
        "can be used to override this default.\n"
        "For a user-defined machine, the default is to create an inputdata\n"
        "directory in the build directory; again, this argument can be\n"
        "used to override this default.",
    )
    non_rebuild_optional_list.append("inputdata-path")

    non_rebuild_optional.add_argument(
        "--no-build",
        action="store_true",
        help="Do the pre-build setup, but do not actually build CTSM\n"
        "(This is useful for testing, or for expert use.)",
    )
    non_rebuild_optional_list.append("no-build")

    new_machine_required = parser.add_argument_group(
        title="required arguments for a user-defined machine",
        description="These arguments are required if neither --machine nor --rebuild are given; "
        "they are not allowed with either of those arguments:",
    )
    new_machine_required_list = []

    new_machine_required.add_argument(
        "--os",
        type=str.lower,
        choices=["linux", "aix", "darwin", "cnl"],
        help="Operating system type",
    )
    new_machine_required_list.append("os")

    new_machine_required.add_argument(
        "--netcdf-path",
        help="Path to NetCDF installation\n"
        "(path to top-level directory, containing subdirectories\n"
        "named lib, include, etc.)",
    )
    new_machine_required_list.append("netcdf-path")

    new_machine_required.add_argument(
        "--esmf-mkfile-path",
        help="Path to esmf.mk file\n" "(typically within ESMF library directory)",
    )
    new_machine_required_list.append("esmf-mkfile-path")

    new_machine_required.add_argument(
        "--max-mpitasks-per-node",
        type=int,
        help="Number of physical processors per shared-memory node\n" "on this machine",
    )
    new_machine_required_list.append("max-mpitasks-per-node")

    new_machine_optional = parser.add_argument_group(
        title="optional arguments for a user-defined machine",
        description="These arguments are optional if neither --machine nor --rebuild are given; "
        "they are not allowed with either of those arguments:",
    )
    new_machine_optional_list = []

    new_machine_optional.add_argument(
        "--gmake",
        default="gmake",
        help="Name of GNU Make tool on your system\n" "Default: gmake",
    )
    new_machine_optional_list.append("gmake")

    new_machine_optional.add_argument(
        "--gmake-j",
        default=8,
        type=int,
        help="Number of threads to use when building\n" "Default: 8",
    )
    new_machine_optional_list.append("gmake-j")

    new_machine_optional.add_argument(
        "--pnetcdf-path",
        help="Path to PNetCDF installation, if present\n"
        "You must either specify this or set --no-pnetcdf",
    )
    new_machine_optional_list.append("pnetcdf-path")

    new_machine_optional.add_argument(
        "--pio-filesystem-hints",
        type=str.lower,
        choices=["gpfs", "lustre"],
        help="Enable filesystem hints for the given filesystem type\n"
        "when building the Parallel IO library",
    )
    new_machine_optional_list.append("pio-filesystem-hints")

    new_machine_optional.add_argument(
        "--gptl-nano-timers",
        action="store_true",
        help="Enable nano timers in build of the GPTL timing library",
    )
    new_machine_optional_list.append("gptl-nano-timers")

    new_machine_optional.add_argument(
        "--extra-fflags",
        default="",
        help="Any extra, non-standard flags to include\n"
        "when compiling Fortran files\n"
        "Tip: to allow a dash at the start of these flags,\n"
        "use a quoted string with an initial space, as in:\n"
        '    --extra-fflags " -flag1 -flag2"',
    )
    new_machine_optional_list.append("extra-fflags")

    new_machine_optional.add_argument(
        "--extra-cflags",
        default="",
        help="Any extra, non-standard flags to include\n"
        "when compiling C files\n"
        "Tip: to allow a dash at the start of these flags,\n"
        "use a quoted string with an initial space, as in:\n"
        '    --extra-cflags " -flag1 -flag2"',
    )
    new_machine_optional_list.append("extra-cflags")

    add_logging_args(parser)

    args = parser.parse_args(args_to_parse)
    if args.rebuild:
        _confirm_args_absent(
            parser,
            args,
            "cannot be provided if --rebuild is set",
            (
                non_rebuild_required_list
                + non_rebuild_optional_list
                + new_machine_required_list
                + new_machine_optional_list
            ),
        )
    else:
        _confirm_args_present(
            parser,
            args,
            "must be provided if --rebuild is not set",
            non_rebuild_required_list,
        )
        if args.machine:
            _confirm_args_absent(
                parser,
                args,
                "cannot be provided if --machine is set",
                new_machine_required_list + new_machine_optional_list,
            )
        else:
            _confirm_args_present(
                parser,
                args,
                "must be provided if neither --machine nor --rebuild are set",
                new_machine_required_list,
            )
            if not args.no_pnetcdf and args.pnetcdf_path is None:
                parser.error(
                    "For a user-defined machine, need to specify either --no-pnetcdf or --pnetcdf-path"
                )
            if args.no_pnetcdf and args.pnetcdf_path is not None:
                parser.error("--no-pnetcdf cannot be given if you set --pnetcdf-path")

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
        arg_no_dashes = arg.replace("-", "_")
        # To determine whether the user specified an argument, we look at whether its
        # value differs from its default value. This won't catch the case where the user
        # explicitly set an argument to its default value, but it's not a big deal if we
        # miss printing an error in that case.
        if vars(args)[arg_no_dashes] != parser.get_default(arg_no_dashes):
            parser.error("--{} {}".format(arg, errmsg))


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
        arg_no_dashes = arg.replace("-", "_")
        if vars(args)[arg_no_dashes] is None:
            parser.error("--{} {}".format(arg, errmsg))


def _check_and_transform_os(os_type):
    """Check validity of os_type argument and transform it to proper case

    os_type should be a lowercase string; returns a transformed string
    """
    transforms = {"linux": "LINUX", "aix": "AIX", "darwin": "Darwin", "cnl": "CNL"}
    try:
        os_type_transformed = transforms[os_type]
    except KeyError as exc:
        raise ValueError("Unknown OS: {}".format(os_type)) from exc
    return os_type_transformed


def _get_case_dir(build_dir):
    """Given the path to build_dir, return the path to the case directory"""
    return os.path.join(build_dir, "case")


def _create_build_dir(build_dir, existing_inputdata):
    """Create the given build directory and any necessary sub-directories

    Args:
    build_dir (str): path to build directory; this directory shouldn't exist yet!
    existing_inputdata (bool): whether the inputdata directory already exists on this machine
    """
    if os.path.exists(build_dir):
        abort(
            "When running without --rebuild, the build directory must not exist yet\n"
            "(<{}> already exists)".format(build_dir)
        )
    os.makedirs(build_dir)
    if not existing_inputdata:
        os.makedirs(os.path.join(build_dir, _INPUTDATA_DIRNAME))


def _fill_out_machine_files(
    *,
    build_dir,
    os_type,
    compiler,
    netcdf_path,
    esmf_mkfile_path,
    max_mpitasks_per_node,
    gmake,
    gmake_j,
    pnetcdf_path=None,
    pio_filesystem_hints=None,
    gptl_nano_timers=False,
    extra_fflags="",
    extra_cflags="",
):
    """Fill out the machine porting templates for this machine / compiler

    For documentation of args, see the documentation in the build_ctsm function
    """
    os.makedirs(os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME, "cmake_macros"))

    # ------------------------------------------------------------------------
    # Fill in config_machines.xml
    # ------------------------------------------------------------------------

    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES, "config_machines_template.xml"),
        path_to_final=os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME, "config_machines.xml"),
        substitutions={
            "OS": os_type,
            "COMPILER": compiler,
            "CIME_OUTPUT_ROOT": build_dir,
            "GMAKE": gmake,
            "GMAKE_J": gmake_j,
            "MAX_TASKS_PER_NODE": max_mpitasks_per_node,
            "MAX_MPITASKS_PER_NODE": max_mpitasks_per_node,
            "ESMF_MKFILE_PATH": esmf_mkfile_path,
        },
    )

    # ------------------------------------------------------------------------
    # Fill in ctsm-build_template.cmake
    # ------------------------------------------------------------------------

    if gptl_nano_timers:
        gptl_cppdefs = _GPTL_NANOTIMERS_CPPDEFS
    else:
        gptl_cppdefs = ""

    if pio_filesystem_hints:
        pio_filesystem_hints_addition = 'set(PIO_FILESYSTEM_HINTS "{}")'.format(
            pio_filesystem_hints
        )
    else:
        pio_filesystem_hints_addition = ""

    if pnetcdf_path:
        pnetcdf_path_addition = 'set(PNETCDF_PATH "{}")'.format(pnetcdf_path)
    else:
        pnetcdf_path_addition = ""

    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES, "ctsm-build_template.cmake"),
        path_to_final=os.path.join(
            build_dir,
            _MACHINE_CONFIG_DIRNAME,
            "cmake_macros",
            "{}_{}.cmake".format(compiler, _MACH_NAME),
        ),
        substitutions={
            "GPTL_CPPDEFS": gptl_cppdefs,
            "NETCDF_PATH": netcdf_path,
            "PIO_FILESYSTEM_HINTS": pio_filesystem_hints_addition,
            "PNETCDF_PATH": pnetcdf_path_addition,
            "EXTRA_CFLAGS": extra_cflags,
            "EXTRA_FFLAGS": extra_fflags,
        },
    )


def _create_case(
    cime_path,
    build_dir,
    compiler,
    *,
    machine=None,
    build_debug=False,
    build_with_openmp=False,
    inputdata_path=None,
):
    """Create a case that can later be used to build the CTSM library and its dependencies

    Args:
    cime_path (str): path to root of cime
    build_dir (str): path to build directory
    compiler (str): compiler to use
    machine (str or None): name of machine or None
        If None, we assume we're using an on-the-fly machine port
        Otherwise, machine should be the name of a machine known to cime
    build_debug (bool): if True, build with flags for debugging
    build_with_openmp (bool): if True, build with OpenMP support
    inputdata_path (str or None): path to existing inputdata directory on this machine
        If None, we use the machine's default DIN_LOC_ROOT
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
    xmlchange = os.path.join(case_dir, "xmlchange")

    if machine is None:
        machine_args = [
            "--machine",
            _MACH_NAME,
            "--extra-machines-dir",
            os.path.join(build_dir, _MACHINE_CONFIG_DIRNAME),
        ]
    else:
        machine_args = ["--machine", machine]

    cmd = os.path.join(cime_path, "scripts", "create_newcase")
    if not os.path.exists(cmd):
        abort(
            "The create_newcase command doesn't exist as expected <{}> does not exist)".format(cmd)
        )
    create_newcase_cmd = [
        cmd,
        "--output-root",
        build_dir,
        "--case",
        case_dir,
        "--compset",
        _COMPSET,
        "--res",
        _RES,
        "--compiler",
        compiler,
        "--driver",
        "nuopc",
        # Project isn't used for anything in the LILAC workflow, but it
        # still needs to be specified on machines that expect it.
        "--project",
        "UNSET",
        "--run-unsupported",
    ]
    create_newcase_cmd.extend(machine_args)
    if inputdata_path:
        create_newcase_cmd.extend(["--input-dir", inputdata_path])
        if not os.path.isdir(inputdata_path):
            abort("inputdata_path directory (<{}> does not exist)".format(inputdata_path))
    run_cmd_output_on_error(
        create_newcase_cmd,
        errmsg="Problem running create_newcase to create the CTSM case directory",
    )

    subprocess.check_call([xmlchange, "LILAC_MODE=on"], cwd=case_dir)
    if build_debug:
        subprocess.check_call([xmlchange, "DEBUG=TRUE"], cwd=case_dir)
    if build_with_openmp:
        subprocess.check_call([xmlchange, "FORCE_BUILD_SMP=TRUE"], cwd=case_dir)

    run_cmd_output_on_error(
        [os.path.join(case_dir, "case.setup")],
        errmsg="Problem setting up CTSM case directory",
        cwd=case_dir,
    )

    make_link(os.path.join(case_dir, "bld"), os.path.join(build_dir, "bld"))
    if machine is not None:
        # For a pre-existing machine, the .env_mach_specific files are likely useful to
        # the user. Make sym links to these with more intuitive names.
        for extension in ("sh", "csh"):
            make_link(
                os.path.join(case_dir, ".env_mach_specific.{}".format(extension)),
                os.path.join(build_dir, "ctsm_build_environment.{}".format(extension)),
            )


def _stage_runtime_inputs(build_dir, no_pnetcdf):
    """Stage CTSM and LILAC runtime inputs

    Args:
    build_dir (str): path to build directory
    no_pnetcdf (bool): if True, use netcdf rather than pnetcdf
    """
    os.makedirs(os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME))

    inputdata_dir = _xmlquery("DIN_LOC_ROOT", build_dir)
    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES, "ctsm_template.cfg"),
        path_to_final=os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME, "ctsm.cfg"),
        substitutions={"INPUTDATA": inputdata_dir},
    )

    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES, "lilac_in_template"),
        path_to_final=os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME, "lilac_in"),
        substitutions={"INPUTDATA": inputdata_dir},
    )

    pio_stride = _xmlquery("MAX_MPITASKS_PER_NODE", build_dir)
    if no_pnetcdf:
        pio_typename = "netcdf"
        # pio_rearranger = 1 is generally more efficient with netcdf (see
        # https://github.com/ESMCI/cime/pull/3732#discussion_r508954806 and the following
        # discussion)
        pio_rearranger = 1
    else:
        pio_typename = "pnetcdf"
        pio_rearranger = 2
    fill_template_file(
        path_to_template=os.path.join(_PATH_TO_TEMPLATES, "lnd_modelio_template.nml"),
        path_to_final=os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME, "lnd_modelio.nml"),
        substitutions={
            "PIO_REARRANGER": pio_rearranger,
            "PIO_STRIDE": pio_stride,
            "PIO_TYPENAME": pio_typename,
        },
    )

    shutil.copyfile(
        src=os.path.join(
            path_to_ctsm_root(), "cime_config", "usermods_dirs", "clm", "lilac", "user_nl_ctsm"
        ),
        dst=os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME, "user_nl_ctsm"),
    )

    make_link(
        _PATH_TO_MAKE_RUNTIME_INPUTS,
        os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME, "make_runtime_inputs"),
    )

    make_link(
        _PATH_TO_DOWNLOAD_INPUT_DATA,
        os.path.join(build_dir, _RUNTIME_INPUTS_DIRNAME, "download_input_data"),
    )


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
            [os.path.join(case_dir, "case.build"), "--sharedlib-only"], cwd=case_dir
        )
    except subprocess.CalledProcessError:
        abort("ERROR building CTSM or its dependencies - see above for details")

    make_link(os.path.join(case_dir, "bld", "ctsm.mk"), os.path.join(build_dir, "ctsm.mk"))


def _xmlquery(varname, build_dir):
    """Run xmlquery from the case in build_dir and return the value of the given variable"""
    case_dir = _get_case_dir(build_dir)
    xmlquery_path = os.path.join(case_dir, "xmlquery")
    value = subprocess.check_output(
        [xmlquery_path, "--value", varname], cwd=case_dir, universal_newlines=True
    )
    return value
