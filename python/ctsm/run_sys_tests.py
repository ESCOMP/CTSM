"""Functions implementing run_sys_tests command"""

from __future__ import print_function
import argparse
import logging
from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args
from ctsm.machine_utils import get_machine_name, get_user
from ctsm.machine import create_machine
from ctsm.machine_defaults import MACHINE_DEFAULTS

logger = logging.getLogger(__name__)

# ========================================================================
# Public functions
# ========================================================================

def main(description):
    """Main function called when run_sys_tests is run from the command-line

    Args:
    description (str): description printed to usage message
    """
    setup_logging_pre_config()
    args = _commandline_args(description)
    process_logging_args(args)
    logger.info('Running on machine: {}'.format(args.machine_name))
    machine = create_machine(machine_name=args.machine_name,
                             defaults=MACHINE_DEFAULTS,
                             scratch_dir=args.testroot_base,
                             account=args.account,
                             job_launcher_queue=args.job_launcher_queue,
                             job_launcher_walltime=args.job_launcher_walltime,
                             job_launcher_extra_args=args.job_launcher_extra_args)
    logger.debug("Machine info: {}".format(machine))
    # FIXME(wjs, 2018-08-24) Fill this in with real arguments
    run_sys_tests(machine=machine, suite_name='aux_clm')

def run_sys_tests(machine, suite_name=None, testfile=None, testlist=None,
                  testid_base=None, test_root_base=None,
                  compare_name=None, generate_name=None, baselineroot=None,
                  account=None, walltime=None, queue=None,
                  extra_create_test_args=''):
    """Implementation of run_sys_tests command

    Args:
    testlist: list of strings giving test names to run
    """
    num_provided_options = ((suite_name is not None) +
                            (testfile is not None) +
                            (testlist is not None and len(testlist) > 0))
    if num_provided_options != 1:
        raise RuntimeError("Exactly one of suite_name, testfile or testlist must be provided")

    if testid_base is None:
        testid_base = _get_testid_base(machine.name)
    if test_root_base is None:
        test_root_base = _get_test_root_base(machine)
    test_root = _get_test_root(test_root_base, testid_base)
    _make_test_root(test_root)

    create_test_args = _get_create_test_args(compare_name=compare_name,
                                             generate_name=generate_name,
                                             baselineroot=baselineroot,
                                             account=account,
                                             walltime=walltime,
                                             queue=queue,
                                             extra_create_test_args=extra_create_test_args)
    if suite_name:
        _run_test_suite(suite_name=suite_name, machine=machine,
                        testid_base=testid_base, test_root=test_root,
                        create_test_args=create_test_args)
    else:
        if testfile:
            test_args = ['--testfile', testfile]
        elif testlist:
            test_args = testlist
        else:
            raise RuntimeError("None of suite_name, testfile or testlist were provided")
        _run_create_test(test_args=test_args, machine=machine,
                         testid=testid_base, test_root=test_root,
                         create_test_args=create_test_args)

# ========================================================================
# Private functions
# ========================================================================

def _commandline_args(description):
    """Parse and return command-line arguments
    """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter)

    machine_name = get_machine_name()

    default_machine = create_machine(machine_name,
                                     defaults=MACHINE_DEFAULTS,
                                     allow_missing_entries=True)

    tests_to_run = parser.add_mutually_exclusive_group(required=True)

    tests_to_run.add_argument('-s', '--suite-name',
                              help='Name of test suite to run')

    tests_to_run.add_argument('-f', '--testfile',
                              help='Path to file listing tests to run')

    tests_to_run.add_argument('-t', '--testname', nargs='+',
                              help='One or more test names to run (space-delimited)')

    compare = parser.add_mutually_exclusive_group(required=True)

    compare.add_argument('-c', '--compare', metavar='COMPARE_NAME',
                        help='Baseline name (often tag) to compare against\n'
                        '(required unless --skip-compare is given)')

    compare.add_argument('--skip-compare', action='store_true',
                        help='Do not compare against baselines')

    generate = parser.add_mutually_exclusive_group(required=True)

    generate.add_argument('-g', '--generate', metavar='GENERATE_NAME',
                        help='Baseline name (often tag) to generate\n'
                        '(required unless --skip-generate is given)')

    generate.add_argument('--skip-generate', action='store_true',
                        help='Do not generate baselines')

    parser.add_argument('--account',
                        help='Account number to use for job submission.\n'
                        'This is needed on some machines; if not provided explicitly,\n'
                        'the script will attempt to guess it using the same rules as in CIME.\n'
                        'Default for this machine: {}'.format(default_machine.account))

    parser.add_argument('--testid-base',
                        help='Base string used for the test id.\n'
                        'Default is to auto-generate this with a date and time stamp.')

    parser.add_argument('--testroot-base',
                        help='Path in which testroot should be put.\n'
                        'For supported machines, this can be left out; for non-supported machines, it must be provided.\n'
                        'Default for this machine: {}'.format(default_machine.scratch_dir))

    parser.add_argument('--baseline-root',
                        help='Path in which baselines should be compared and generated.\n'
                        'Default is to use the default for this machine.\n')

    parser.add_argument('--walltime',
                        help='Walltime for each test.\n'
                        'If running a test suite, you can generally leave this unset,\n'
                        'because it is set in the file defining the test suite.\n'
                        'For other uses, providing this helps decrease the time spent waiting in the queue.')

    parser.add_argument('--queue',
                        help='Queue to which tests are submitted.\n'
                        'If not provided, uses machine default.')

    parser.add_argument('--extra-create-test-args',
                        help='String giving extra arguments to pass to create_test')


    parser.add_argument('--job-launcher-queue',
                        help='Queue to which the create_test command is submitted.\n'
                        'Only applies on machines for which we submit the create_test command\n'
                        'rather than running it on the login node.\n'
                        'Default for this machine: {}'.format(default_machine.job_launcher.get_queue()))

    parser.add_argument('--job-launcher-walltime',
                        help='Walltime for the create_test command.\n'
                        'Only applies on machines for which we submit the create_test command\n'
                        'rather than running it on the login node.\n'
                        'Default for this machine: {}'.format(default_machine.job_launcher.get_walltime()))

    parser.add_argument('--job-launcher-extra-args',
                        help='Extra arguments for the command that launches the create_test command.\n'
                        'Default for this machine: {}'.format(default_machine.job_launcher.get_extra_args()))

    parser.add_argument('--machine-name', default=machine_name,
                        help='Name of machine for which create_test is run.\n'
                        'This typically is not needed, but can be provided for the sake of testing this script.\n'
                        'Defaults to current machine: {}'.format(machine_name))

    add_logging_args(parser)

    args = parser.parse_args()

    return args

def _get_testid_base(machine_name):
    """Returns a base testid based on the current date and time and the machine name"""
    # FIXME(wjs, 2018-08-24) Implement this: use the current month, day, hour & minute
    # plus the first two characters of the machine name
    return ""

def _get_test_root_base(machine):
    # FIXME(wjs, 2018-08-24) Implement this
    return ""

def _get_test_root(test_root_base, testid_base):
    """Get the path to the test root, given a base test id"""
    # FIXME(wjs, 2018-08-24) Implement this
    return ""

def _make_test_root(test_root):
    """Make the testroot directory at the given location, as well as a link in the current
    directory
    """
    # FIXME(wjs, 2018-08-24) Implement this.
    #
    # This will echo an equivalent shell command, then (if dry_run is False) execute the
    # python equivalent.
    #
    # I think this is also where we should create the cs.status scripts
    pass

def _get_create_test_args(compare_name, generate_name, baselineroot,
                          account, walltime, queue,
                          extra_create_test_args):
    args = []
    if compare_name:
        args.extend(['--compare', compare_name])
    if generate_name:
        args.extend(['--generate', generate_name])
    if baselineroot:
        args.extend(['--baseline-root', baselineroot])
    if account:
        args.extend(['--project', account])
    if walltime:
        args.extend(['--walltime', walltime])
    if queue:
        args.extend(['--queue', queue])
    args.extend(extra_create_test_args.split())
    return args

def _run_test_suite(suite_name, machine, testid_base, test_root, create_test_args):

    compilers = _get_compilers_for_suite(suite_name, machine.name)
    for compiler in compilers:
        test_args = ['--xml-category', suite_name,
                     '--xml-machine', machine.name,
                     '--xml-compiler', compiler]
        testid = testid_base + '_' + compiler[0:2]
        _run_create_test(test_args=test_args, machine=machine,
                         testid=testid, test_root=test_root,
                         create_test_args=create_test_args)

def _get_compilers_for_suite(suite_name, machine_name):
    # FIXME(wjs, 2018-08-24) Implement this. Use cime.test_status, similar to what's done
    # in query_testlists. If no tests are found for this suite and machine, raise an exception.
    return []

def _run_create_test(test_args, machine, testid, test_root, create_test_args):
    # FIXME(wjs, 2018-08-24) Finish this implementation
    create_test_cmd = _build_create_test_cmd()

def _build_create_test_cmd():
    """Builds and returns the create_test command"""
    # FIXME(wjs, 2018-08-24) Implement this
    return ""

