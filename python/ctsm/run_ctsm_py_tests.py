"""Runner for the python unit tests defined here

This is the main implementation of the run_ctsm_py_tests script contained in the
parent directory
"""

import unittest
import os
import argparse
import logging
from ctsm import unit_testing

logger = logging.getLogger(__name__)


def main(description):
    """Main function called when run_tests is run from the command-line

    Args:
    description (str): description printed to usage message
    """
    args = _commandline_args(description)
    verbosity = _get_verbosity_level(args)

    if args.pattern is not None:
        pattern = args.pattern
    elif args.unit:
        pattern = "test_unit*.py"
    elif args.sys:
        pattern = "test_sys*.py"
    else:
        pattern = "test*.py"

    # This setup_for_tests call is the main motivation for having this wrapper script to
    # run the tests rather than just using 'python -m unittest discover'
    unit_testing.setup_for_tests(enable_critical_logs=args.debug)

    mydir = os.path.dirname(os.path.abspath(__file__))
    testsuite = unittest.defaultTestLoader.discover(start_dir=mydir, pattern=pattern)
    # NOTE(wjs, 2018-08-29) We may want to change the meaning of '--debug'
    # vs. '--verbose': I could imagine having --verbose set buffer=False, and --debug
    # additionally sets the logging level to much higher - e.g., debug level.
    testrunner = unittest.TextTestRunner(buffer=(not args.debug), verbosity=verbosity)
    testrunner.run(testsuite)


def _commandline_args(description):
    """Parse and return command-line arguments
    Note that run_ctsm_py_tests is not intended to be
    used without argument specifications
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )

    output_level = parser.add_mutually_exclusive_group()

    output_level.add_argument(
        "-v", "--verbose", action="store_true", help="Run tests with more verbosity"
    )

    output_level.add_argument(
        "-d", "--debug", action="store_true", help="Run tests with even more verbosity"
    )

    test_subset = parser.add_mutually_exclusive_group()

    test_subset.add_argument("-u", "--unit", action="store_true", help="Only run unit tests")

    test_subset.add_argument("-s", "--sys", action="store_true", help="Only run system tests")

    test_subset.add_argument(
        "-p", "--pattern", help="File name pattern to match\n" "Default is test*.py"
    )

    args = parser.parse_args()

    return args


def _get_verbosity_level(args):
    if args.debug or args.verbose:
        verbosity = 2
    else:
        verbosity = 1
    return verbosity
