"""Functions implementing run_sys_tests command"""

import argparse
import logging
import os
import sys
import subprocess
from datetime import datetime

from CIME.test_utils import get_tests_from_xml  # pylint: disable=import-error
from CIME.test_utils import test_to_string  # pylint: disable=import-error
from CIME.cs_status_creator import create_cs_status  # pylint: disable=import-error

from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)
from ctsm.machine_utils import get_machine_name
from ctsm.machine import (
    create_machine,
    get_possibly_overridden_mach_value,
    CREATE_TEST_QUEUE_UNSPECIFIED,
)
from ctsm.machine_defaults import MACHINE_DEFAULTS
from ctsm.os_utils import make_link
from ctsm.path_utils import path_to_ctsm_root
from ctsm.joblauncher.job_launcher_factory import JOB_LAUNCHER_NOBATCH

logger = logging.getLogger(__name__)

# Number of initial characters from the compiler name to use in a testid
_NUM_COMPILER_CHARS = 3

# For job launchers that use 'nice', the level of niceness we should use
_NICE_LEVEL = 19

# Extra arguments for the cs.status.fails command
_CS_STATUS_FAILS_EXTRA_ARGS = "--fails-only --count-performance-fails"

# ========================================================================
# Public functions
# ========================================================================


def main(cime_path):
    """Main function called when run_sys_tests is run from the command-line

    Args:
    cime_path (str): path to the cime that we're using (this is passed in explicitly
        rather than relying on calling path_to_cime so that we can be absolutely sure that
        the scripts called here are coming from the same cime as the cime library we're
        using).
    """
    setup_logging_pre_config()
    args = _commandline_args()
    process_logging_args(args)
    logger.info("Running on machine: %s", args.machine_name)
    if args.job_launcher_nobatch:
        job_launcher_type = JOB_LAUNCHER_NOBATCH
    else:
        job_launcher_type = None
    machine = create_machine(
        machine_name=args.machine_name,
        job_launcher_type=job_launcher_type,
        defaults=MACHINE_DEFAULTS,
        account=args.account,
        job_launcher_queue=args.job_launcher_queue,
        job_launcher_walltime=args.job_launcher_walltime,
        job_launcher_nice_level=_NICE_LEVEL,
        job_launcher_extra_args=args.job_launcher_extra_args,
    )
    logger.debug("Machine info: %s", machine)

    run_sys_tests(
        machine=machine,
        cime_path=cime_path,
        skip_testroot_creation=args.skip_testroot_creation,
        skip_git_status=args.skip_git_status,
        dry_run=args.dry_run,
        suite_name=args.suite_name,
        testfile=args.testfile,
        testlist=args.testname,
        suite_compilers=args.suite_compiler,
        testid_base=args.testid_base,
        testroot_base=args.testroot_base,
        rerun_existing_failures=args.rerun_existing_failures,
        compare_name=args.compare,
        generate_name=args.generate,
        baseline_root=args.baseline_root,
        walltime=args.walltime,
        queue=args.queue,
        retry=args.retry,
        extra_create_test_args=args.extra_create_test_args,
    )


def run_sys_tests(
    machine,
    cime_path,
    *,
    skip_testroot_creation=False,
    skip_git_status=False,
    dry_run=False,
    suite_name=None,
    testfile=None,
    testlist=None,
    suite_compilers=None,
    testid_base=None,
    testroot_base=None,
    rerun_existing_failures=False,
    compare_name=None,
    generate_name=None,
    baseline_root=None,
    walltime=None,
    queue=None,
    retry=None,
    extra_create_test_args="",
):
    """Implementation of run_sys_tests command

    Exactly one of suite_name, testfile or testlist should be provided

    Args:
    machine: Machine object, as defined in machine.py
    cime_path (str): path to root of cime
    skip_testroot_creation (bool): if True, assume the testroot directory has already been
        created, so don't try to recreate it or re-make the link to it
    skip_git_status (bool): if True, skip printing git and git-fleximod status
    dry_run (bool): if True, print commands to be run but don't run them
    suite_name (str): name of test suite/category to run
    testfile (str): path to file containing list of tests to run
    testlist (list of strings): names of tests to run
    suite_compilers (list of strings): compilers to use in the test suite; only applicable
        with suite_name; if not specified, use all compilers that are defined for this
        test suite
    testid_base (str): test id, or start of the test id in the case of a test suite (if
        not provided, will be generated automatically)
    testroot_base (str): path to the directory that will contain the testroot (if not
        provided, will be determined based on machine defaults)
    rerun_existing_failures (bool): if True, add the '--use-existing' option to create_test
        If specified, then --testid-base should also be specified. This also implies
        skip_testroot_creation.
    compare_name (str): if not None, baseline name to compare against
    generate_name (str): if not None, baseline name to generate
    baseline_root (str): path in which baselines should be compared and generated (if not
        provided, this will be obtained from the default value in the machine object; if
        that is None, then the test suite will determine it automatically)
    walltime (str): walltime to use for each test (if not provided, the test suite will
        determine it automatically)
    queue (str): queue to use for each test (if not provided, will use the default for
        this machine based on the passed-in machine object; if that is unspecified, then
        the test suite will determine it automatically)
    retry (int): retry value to pass to create_test (if not provided, will use the default
        for this machine)
    extra_create_test_args (str): any extra arguments to create_test, as a single,
        space-delimited string
    testlist: list of strings giving test names to run

    """
    num_provided_options = (
        (suite_name is not None)
        + (testfile is not None)
        + (testlist is not None and len(testlist) > 0)
    )
    if num_provided_options != 1:
        raise RuntimeError("Exactly one of suite_name, testfile or testlist must be provided")

    if testid_base is None:
        testid_base = _get_testid_base(machine.name)
    if testroot_base is None:
        testroot_base = _get_testroot_base(machine)
    testroot = _get_testroot(testroot_base, testid_base)
    if not (skip_testroot_creation or rerun_existing_failures):
        _make_testroot(testroot, testid_base, dry_run)
    else:
        if not os.path.exists(testroot):
            raise RuntimeError(
                "testroot directory does NOT exist as expected when a rerun"
                + " option is used: directory expected = "
                + testroot
            )
    print("Testroot: {}\n".format(testroot))
    retry_final = get_possibly_overridden_mach_value(
        machine, varname="create_test_retry", value=retry
    )
    # Note the distinction between a queue of None and a queue of
    # CREATE_TEST_QUEUE_UNSPECIFIED in the following: If queue is None (meaning that the
    # user hasn't specified a '--queue' argument to run_sys_tests), then we'll use the
    # queue specified in the machine object; if queue is CREATE_TEST_QUEUE_UNSPECIFIED,
    # then we'll force queue_final to be None, which means we won't add a '--queue'
    # argument to create_test, regardless of what is specified in the machine object.
    # (It's also possible for the machine object to specify a queue of
    # CREATE_TEST_QUEUE_UNSPECIFIED, which means that we won't use a '--queue' argument to
    # create_test unless the user specifies a '--queue' argument to run_sys_tests.)
    queue_final = get_possibly_overridden_mach_value(
        machine, varname="create_test_queue", value=queue
    )
    if queue_final == CREATE_TEST_QUEUE_UNSPECIFIED:
        queue_final = None
    if not skip_git_status:
        _record_git_status(testroot, retry_final, dry_run)

    baseline_root_final = get_possibly_overridden_mach_value(
        machine, varname="baseline_dir", value=baseline_root
    )
    create_test_args = _get_create_test_args(
        compare_name=compare_name,
        generate_name=generate_name,
        baseline_root=baseline_root_final,
        account=machine.account,
        walltime=walltime,
        queue=queue_final,
        retry=retry_final,
        rerun_existing_failures=rerun_existing_failures,
        extra_create_test_args=extra_create_test_args,
    )

    running_ctsm_py_tests = (
        testfile == "/path/to/testfile"
        or testlist in [["test1", "test2"], ["foo"]]
        or suite_name == "my_suite"
    )

    if suite_name:
        if not dry_run:
            _make_cs_status_for_suite(testroot, testid_base)
        _run_test_suite(
            cime_path=cime_path,
            suite_name=suite_name,
            suite_compilers=suite_compilers,
            machine=machine,
            testid_base=testid_base,
            testroot=testroot,
            create_test_args=create_test_args,
            dry_run=dry_run,
            running_ctsm_py_tests=running_ctsm_py_tests,
        )
    else:
        if not dry_run:
            _make_cs_status_non_suite(testroot, testid_base)
        running_ctsm_py_tests = testfile == "/path/to/testfile"
        testname_list = None
        if testfile:
            test_args = ["--testfile", os.path.abspath(testfile)]
            if not running_ctsm_py_tests:
                with open(test_args[1], "r") as testfile_abspath:
                    testname_list = testfile_abspath.readlines()
        elif testlist:
            test_args = testlist
            testname_list = testlist
        else:
            raise RuntimeError("None of suite_name, testfile or testlist were provided")
        if testname_list:
            _check_py_env(testname_list)
        _run_create_test(
            cime_path=cime_path,
            test_args=test_args,
            machine=machine,
            testid=testid_base,
            testroot=testroot,
            create_test_args=create_test_args,
            dry_run=dry_run,
        )


# ========================================================================
# Private functions
# ========================================================================


def _commandline_args():
    """Parse and return command-line arguments"""

    description = """
Driver for running CTSM system tests

Typical usage:

./run_sys_tests -s aux_clm -c COMPARE_NAME -g GENERATE_NAME

    This automatically detects the machine and launches the appropriate components of the
    aux_clm test suite on that machine. This script also implements other aspects of the
    typical CTSM system testing workflow, such as running create_test via qsub on
    cheyenne, and setting up a directory to hold all of the tests in the test suite. A
    symbolic link will be created in the current directory pointing to the testroot
    directory containing all of the test directories in the test suite.

    Note that the -c/--compare and -g/--generate arguments are required, unless you specify
    --skip-compare and/or --skip-generate.

    Any other test suite can be given as well: clm_short, aux_glc, etc.

This can also be used to run tests listed in a text file (via the -f/--testfile argument),
or tests listed individually on the command line (via the -t/--testname argument).
"""

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )

    machine_name = get_machine_name()

    default_machine = create_machine(
        machine_name, defaults=MACHINE_DEFAULTS, allow_missing_entries=True
    )

    tests_to_run = parser.add_mutually_exclusive_group(required=True)

    tests_to_run.add_argument("-s", "--suite-name", help="Name of test suite to run")

    tests_to_run.add_argument("-f", "--testfile", help="Path to file listing tests to run")

    tests_to_run.add_argument(
        "-t",
        "--testname",
        "--testnames",
        nargs="+",
        help="One or more test names to run (space-delimited)",
    )

    compare = parser.add_mutually_exclusive_group(required=True)

    compare.add_argument(
        "-c",
        "--compare",
        metavar="COMPARE_NAME",
        help="Baseline name (often tag) to compare against\n"
        "(required unless --skip-compare is given)",
    )

    compare.add_argument(
        "--skip-compare", action="store_true", help="Do not compare against baselines"
    )

    generate = parser.add_mutually_exclusive_group(required=True)

    generate.add_argument(
        "-g",
        "--generate",
        metavar="GENERATE_NAME",
        help="Baseline name (often tag) to generate\n" "(required unless --skip-generate is given)",
    )

    generate.add_argument("--skip-generate", action="store_true", help="Do not generate baselines")

    parser.add_argument(
        "--suite-compiler",
        "--suite-compilers",
        nargs="+",
        help="Compiler(s) from the given test suite for which tests are run\n"
        "Only valid in conjunction with -s/--suite-name;\n"
        "if not specified, use all compilers defined for this suite and machine\n",
    )

    parser.add_argument(
        "--account",
        help="Account number to use for job submission.\n"
        "This is needed on some machines; if not provided explicitly,\n"
        "the script will attempt to guess it using the same rules as in CIME.\n"
        "Default for this machine: {}".format(default_machine.account),
    )

    parser.add_argument(
        "--testid-base",
        help="Base string used for the test id.\n"
        "Default is to auto-generate this with a date and time stamp.",
    )

    parser.add_argument(
        "--testroot-base",
        help="Path in which testroot should be put.\n"
        "For supported machines, this can be left out;\n"
        "for non-supported machines, it must be provided.\n"
        "Default for this machine: {}".format(default_machine.scratch_dir),
    )

    parser.add_argument(
        "--rerun-existing-failures",
        action="store_true",
        help="Rerun failed tests from the last PEND or FAIL state.\n"
        "This triggers the --use-existing option to create_test.\n"
        "To use this option, provide the same options to run_sys_tests\n"
        "as in the initial run, but also adding --testid-base\n"
        "corresponding to the base testid used initially.\n"
        "(However, many of the arguments to create_test are ignored,\n"
        "so it is not important for all of the options to exactly match\n"
        "those in the initial run.)\n"
        "This option implies --skip-testroot-creation (that option does not\n"
        "need to be specified separately if using --rerun-existing-failures).",
    )

    if default_machine.baseline_dir:
        baseline_root_default_msg = "Default for this machine: {}".format(
            default_machine.baseline_dir
        )
    else:
        baseline_root_default_msg = "Default for this machine: use cime's default"
    parser.add_argument(
        "--baseline-root",
        help="Path in which baselines should be compared and generated.\n"
        + baseline_root_default_msg,
    )

    parser.add_argument(
        "--walltime",
        help="Walltime for each test.\n"
        "If running a test suite, you can generally leave this unset,\n"
        "because it is set in the file defining the test suite.\n"
        "For other uses, providing this helps decrease the time spent\n"
        "waiting in the queue.",
    )

    parser.add_argument(
        "--queue",
        help="Queue to which tests are submitted.\n"
        'The special value "{}" means do not add a --queue option to create_test,\n'
        "instead allowing CIME to pick an appropriate queue for each test\n"
        "using its standard mechanisms.\n"
        "Default for this machine: {}".format(
            CREATE_TEST_QUEUE_UNSPECIFIED, default_machine.create_test_queue
        ),
    )

    parser.add_argument(
        "--retry",
        type=int,
        help="Argument to create_test: Number of times to retry failed tests.\n"
        "Default for this machine: {}".format(default_machine.create_test_retry),
    )

    parser.add_argument(
        "--extra-create-test-args",
        default="",
        help="String giving extra arguments to pass to create_test\n"
        "(To allow the argument parsing to accept this, enclose the string\n"
        'in quotes, with a leading space, as in " --my-arg foo".)',
    )

    parser.add_argument(
        "--job-launcher-nobatch",
        action="store_true",
        help="Run create_test on the login node, even if this machine\n"
        "is set up to submit create_test to a compute node by default.",
    )

    parser.add_argument(
        "--job-launcher-queue",
        help="Queue to which the create_test command is submitted.\n"
        "Only applies on machines for which we submit the create_test command\n"
        "rather than running it on the login node.\n"
        "Default for this machine: {}".format(default_machine.job_launcher.get_queue()),
    )

    parser.add_argument(
        "--job-launcher-walltime",
        help="Walltime for the create_test command.\n"
        "Only applies on machines for which we submit the create_test command\n"
        "rather than running it on the login node.\n"
        "Default for this machine: {}".format(default_machine.job_launcher.get_walltime()),
    )

    parser.add_argument(
        "--job-launcher-extra-args",
        help="Extra arguments for the command that launches the\n"
        "create_test command.\n"
        "(To allow the argument parsing to accept this, enclose the string\n"
        'in quotes, with a leading space, as in " --my-arg foo".)\n'
        "Default for this machine: {}".format(default_machine.job_launcher.get_extra_args()),
    )

    parser.add_argument(
        "--skip-testroot-creation",
        action="store_true",
        help="Do not create the directory that will hold the tests.\n"
        "This should be used if the desired testroot directory already exists.",
    )

    parser.add_argument(
        "--skip-git-status",
        action="store_true",
        help="Skip printing git and git-fleximod status,\n"
        "both to screen and to the SRCROOT_GIT_STATUS file in TESTROOT.\n"
        "This printing can often be helpful, but this option can be used to\n"
        "avoid extraneous output, to reduce the time needed to run this script,\n"
        "or if git or git-fleximod are currently broken in your sandbox.\n",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would happen, but do not run any commands.\n"
        "(Generally should be run with --verbose.)\n",
    )

    parser.add_argument(
        "--machine-name",
        default=machine_name,
        help="Name of machine for which create_test is run.\n"
        "This typically is not needed, but can be provided\n"
        "for the sake of testing this script.\n"
        "Defaults to current machine: {}".format(machine_name),
    )

    add_logging_args(parser)

    args = parser.parse_args()

    _check_arg_validity(args)

    return args


def _check_arg_validity(args):
    if args.suite_compiler and not args.suite_name:
        raise RuntimeError("--suite-compiler can only be specified if using --suite-name")
    if args.rerun_existing_failures and not args.testid_base:
        raise RuntimeError("With --rerun-existing-failures, must also specify --testid-base")


def _get_testid_base(machine_name):
    """Returns a base testid based on the current date and time and the machine name"""
    now = datetime.now()
    now_str = now.strftime("%m%d-%H%M%S")
    machine_start = machine_name[0:2]
    return "{}{}".format(now_str, machine_start)


def _get_testroot_base(machine):
    scratch_dir = machine.scratch_dir
    if scratch_dir is None:
        raise RuntimeError(
            "For a machine without a default specified for the scratch directory, "
            "must specify --testroot-base"
        )
    return scratch_dir


def _get_testroot(testroot_base, testid_base):
    """Get the path to the test root, given a base test id"""
    return os.path.join(testroot_base, _get_testdir_name(testid_base))


def _get_testdir_name(testid_base):
    return "tests_{}".format(testid_base)


def _make_testroot(testroot, testid_base, dry_run):
    """Make the testroot directory at the given location, as well as a link in the current
    directory
    """
    if os.path.exists(testroot):
        raise RuntimeError("{} already exists".format(testroot))
    logger.info("Making directory: %s", testroot)
    if not dry_run:
        os.makedirs(testroot)
        make_link(testroot, _get_testdir_name(testid_base))


def _record_git_status(testroot, retry, dry_run):
    """Record git status and related information to stdout and a file"""
    output = ""
    ctsm_root = path_to_ctsm_root()

    output += "create_test --retry: {}\n\n".format(retry)

    current_hash = subprocess.check_output(
        ["git", "show", "--no-patch", "--format=format:%h (%an, %ad) %s\n", "HEAD"],
        cwd=ctsm_root,
        universal_newlines=True,
    )
    output += "Current hash: {}".format(current_hash)
    git_status = subprocess.check_output(
        ["git", "-c", "color.ui=always", "status", "--short", "--branch"],
        cwd=ctsm_root,
        universal_newlines=True,
    )
    output += git_status
    if git_status.count("\n") == 1:
        # Only line in git status is the branch info
        output += "(clean sandbox)\n"
    fleximod = os.path.join("bin", "git-fleximod")
    fleximod_status = subprocess.check_output(
        [fleximod, "status"], cwd=ctsm_root, universal_newlines=True
    )
    output += 72 * "-" + "\n" + "git-fleximod status:" + "\n"
    output += fleximod_status
    output += 72 * "-" + "\n"

    print(output)

    if not dry_run:
        git_status_filepath = os.path.join(testroot, "SRCROOT_GIT_STATUS")
        if os.path.exists(git_status_filepath):
            # If we're reusing an existing directory, it could happen that
            # SRCROOT_GIT_STATUS already exists. It's still helpful to record the current
            # SRCROOT_GIT_STATUS information, but we don't want to clobber the old. So
            # make a new file with a date/time-stamp.
            now = datetime.now()
            now_str = now.strftime("%m%d-%H%M%S")
            git_status_filepath = git_status_filepath + "_" + now_str
        with open(git_status_filepath, "w") as git_status_file:
            git_status_file.write(" ".join(sys.argv) + "\n\n")
            git_status_file.write("SRCROOT: {}\n".format(ctsm_root))
            git_status_file.write(output)


def _get_create_test_args(
    *,
    compare_name,
    generate_name,
    baseline_root,
    account,
    walltime,
    queue,
    retry,
    rerun_existing_failures,
    extra_create_test_args,
):
    args = []
    if compare_name:
        args.extend(["--compare", compare_name])
    if generate_name:
        args.extend(["--generate", generate_name])
    if baseline_root:
        args.extend(["--baseline-root", baseline_root])
    if account:
        args.extend(["--project", account])
    if walltime:
        args.extend(["--walltime", walltime])
    if queue:
        args.extend(["--queue", queue])
    args.extend(["--retry", str(retry)])
    if rerun_existing_failures:
        # In addition to --use-existing, we also need --allow-baseline-overwrite in this
        # case; otherwise, create_test throws an error saying that the baseline
        # directories already exist.
        args.extend(["--use-existing", "--allow-baseline-overwrite"])
    args.extend(extra_create_test_args.split())
    return args


def _make_cs_status_for_suite(testroot, testid_base):
    """Makes a cs.status file that can be run for the entire test suite"""
    testid_pattern = testid_base + "_" + _NUM_COMPILER_CHARS * "?"
    # The basic cs.status just aggregates results from all of the individual create_tests
    create_cs_status(
        test_root=testroot,
        test_id=testid_pattern,
        extra_args=_cs_status_xfail_arg(),
        filename="cs.status",
    )
    # cs.status.fails additionally filters the results so that only failures are shown
    create_cs_status(
        test_root=testroot,
        test_id=testid_pattern,
        extra_args=(_CS_STATUS_FAILS_EXTRA_ARGS + " " + _cs_status_xfail_arg()),
        filename="cs.status.fails",
    )


def _make_cs_status_non_suite(testroot, testid_base):
    """Makes a cs.status file for a single run of create_test - not a whole test suite"""
    create_cs_status(
        test_root=testroot,
        test_id=testid_base,
        extra_args=(_CS_STATUS_FAILS_EXTRA_ARGS + " " + _cs_status_xfail_arg()),
        filename="cs.status.fails",
    )


def _cs_status_xfail_arg():
    """Returns a string giving the argument to cs_status that will point to CTSM's
    expected fails xml file
    """
    ctsm_root = path_to_ctsm_root()
    xfail_path = os.path.join(ctsm_root, "cime_config", "testdefs", "ExpectedTestFails.xml")
    return "--expected-fails-file {}".format(xfail_path)


def _run_test_suite(
    *,
    cime_path,
    suite_name,
    suite_compilers,
    machine,
    testid_base,
    testroot,
    create_test_args,
    dry_run,
    running_ctsm_py_tests,
):
    if not suite_compilers:
        suite_compilers = _get_compilers_for_suite(suite_name, machine.name, running_ctsm_py_tests)
    for compiler in suite_compilers:
        test_args = [
            "--xml-category",
            suite_name,
            "--xml-machine",
            machine.name,
            "--xml-compiler",
            compiler,
        ]
        testid = testid_base + "_" + compiler[0:_NUM_COMPILER_CHARS]
        _run_create_test(
            cime_path=cime_path,
            test_args=test_args,
            machine=machine,
            testid=testid,
            testroot=testroot,
            create_test_args=create_test_args,
            dry_run=dry_run,
        )


def _get_testmod_list(test_attributes, unique=False):
    # Isolate testmods, producing a list like
    # ["clm-test1mod1", "clm-test2mod1", "clm-test2mod2", ...]
    # Handles test attributes passed in from run_sys_tests calls using -t, -f, or -s

    testmods = []
    for test_attribute in test_attributes:
        for dot_split in test_attribute.split("."):
            slash_replaced = dot_split.replace("/", "-")
            for ddash_split in slash_replaced.split("--"):
                if "clm-" in ddash_split and (ddash_split not in testmods or not unique):
                    testmods.append(ddash_split)

    return testmods


def _check_py_env(test_attributes):
    err_msg = " can't be loaded. Do you need to activate the ctsm_pylib conda environment?"
    # Suppress pylint import-outside-toplevel warning because (a) we only want to import
    # this when certain tests are requested, and (b) the import needs to be in a try-except
    # block to produce a nice error message.
    # pylint: disable=import-outside-toplevel disable
    # Suppress pylint unused-import warning because the import itself IS the use.
    # pylint: disable=unused-import disable
    # Suppress pylint import-error warning because the whole point here is to check
    # whether import is possible.
    # pylint: disable=import-error disable

    # Check requirements for using modify_fsurdat Python module, if needed
    modify_fsurdat_users = ["FSURDATMODIFYCTSM", "RXCROPMATURITY"]
    if any(any(u in t for u in modify_fsurdat_users) for t in test_attributes):
        try:
            import ctsm.modify_input_files.modify_fsurdat
        except ModuleNotFoundError as err:
            raise ModuleNotFoundError("modify_fsurdat" + err_msg) from err

    # Check requirements for RXCROPMATURITY, if needed
    if any("RXCROPMATURITY" in t for t in test_attributes):
        try:
            import ctsm.crop_calendars.check_rxboth_run
        except ModuleNotFoundError as err:
            raise ModuleNotFoundError("check_rxboth_run" + err_msg) from err
        try:
            import ctsm.crop_calendars.generate_gdds
        except ModuleNotFoundError as err:
            raise ModuleNotFoundError("generate_gdds" + err_msg) from err
        try:
            import ctsm.crop_calendars.interpolate_gdds
        except ModuleNotFoundError as err:
            raise ModuleNotFoundError("interpolate_gdds" + err_msg) from err

    # Check that list for any testmods that use modify_fates_paramfile.py
    testmods_to_check = ["clm-FatesColdTwoStream", "clm-FatesColdTwoStreamNoCompFixedBioGeo"]
    testmods = _get_testmod_list(test_attributes)
    if any(t in testmods_to_check for t in testmods):
        # This bit is needed because it's outside the top-level python/ directory.
        fates_dir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "src", "fates"
        )
        sys.path.insert(1, fates_dir)
        try:
            import tools.modify_fates_paramfile
        except ModuleNotFoundError as err:
            raise ModuleNotFoundError("modify_fates_paramfile" + err_msg) from err


def _get_compilers_for_suite(suite_name, machine_name, running_ctsm_py_tests):
    test_data = get_tests_from_xml(xml_machine=machine_name, xml_category=suite_name)
    if not test_data:
        raise RuntimeError(
            "No tests found for suite {} on machine {}".format(suite_name, machine_name)
        )
    if logger.getEffectiveLevel() <= logging.INFO:
        logger.info("Tests:")
        for test in test_data:
            test_string = test_to_string(test).split(" ")[1]
            logger.info("   %s", test_string)
    if not running_ctsm_py_tests:
        _check_py_env([t["testname"] for t in test_data])
        _check_py_env([t["testmods"] for t in test_data if "testmods" in t.keys()])
    compilers = sorted({one_test["compiler"] for one_test in test_data})
    logger.info("Running with compilers: %s", compilers)
    return compilers


def _run_create_test(*, cime_path, test_args, machine, testid, testroot, create_test_args, dry_run):
    create_test_cmd = _build_create_test_cmd(
        cime_path=cime_path,
        test_args=test_args,
        testid=testid,
        testroot=testroot,
        create_test_args=create_test_args,
    )
    stdout_path = os.path.join(testroot, "STDOUT.{}".format(testid))
    stderr_path = os.path.join(testroot, "STDERR.{}".format(testid))
    machine.job_launcher.run_command(
        create_test_cmd,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        dry_run=dry_run,
    )


def _build_create_test_cmd(cime_path, test_args, testid, testroot, create_test_args):
    """Builds and returns the create_test command

    This is a list, where each element of the list is one argument
    """
    command = [
        os.path.join(cime_path, "scripts", "create_test"),
        "--test-id",
        testid,
        "--output-root",
        testroot,
    ]
    command.extend(test_args)
    command.extend(create_test_args)
    return command
