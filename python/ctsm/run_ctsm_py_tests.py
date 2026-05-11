"""Runner for the python unit tests defined here

This is the main implementation of the run_ctsm_py_tests script contained in the
parent directory
"""

import os
import glob
import argparse
import logging
import io
import contextlib
import re

import pytest

from ctsm import unit_testing

logger = logging.getLogger(__name__)

# Helpful message explaining the fact that our -s differs from pytest -s
SYS_TESTS_DISAMBIGUATION = "If you want to use pytest's -s option, use --capture=no instead."


def _get_files_matching_pattern(pattern):
    pattern = os.path.join("**", pattern)
    result = glob.glob(pattern, recursive=True)
    result.sort()
    result = [f for f in result if f.endswith(".py")]
    return result


def main(description):
    """Main function called when run_tests is run from the command-line

    Args:
    description (str): description printed to usage message
    """
    args, pytest_args = _commandline_args(description)

    # Get list of test files to process, if any requested
    file_list = []
    if args.pattern is not None:
        file_list += _get_files_matching_pattern(args.pattern)
    if args.unit:
        file_list += _get_files_matching_pattern("test_unit*.py")
    if args.sys:
        file_list += _get_files_matching_pattern("test_sys*.py")
    pytest_args += file_list

    # This setup_for_tests call is the main motivation for having this wrapper script to
    # run the tests rather than just using 'python -m pytest'
    unit_testing.setup_for_tests(enable_critical_logs=args.debug)

    # Run the tests
    pytest.main(pytest_args)


def _commandline_args(description):
    """Parse and return command-line arguments
    Note that run_ctsm_py_tests is not intended to be
    used without argument specifications
    """

    # Get help for pytest options we're using to overload existing run_ctsm_py_test options
    debug_help = _get_pytest_help()

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Run tests with more verbosity"
    )

    parser.add_argument(
        "-d",
        "--debug",
        nargs="?",
        type=str,
        const="",
        metavar="DEBUG_FILE_NAME",  # Same as pytest
        help=(
            "If given with no argument, uses old unittest-based run_ctsm_py_tests behavior: Run"
            " tests with maximum verbosity, equivalent to ``pylint -v --capture=no``. If given"
            f" with argument, uses pytest behavior: {debug_help}"
        ),
    )

    parser.add_argument(
        "--no-disable-warnings", action="store_true", help="Show pytest's warnings summary"
    )

    test_subset = parser.add_mutually_exclusive_group()

    test_subset.add_argument("-u", "--unit", action="store_true", help="Only run unit tests")

    test_subset.add_argument(
        "-s",
        "--sys",
        action="store_true",
        help=f"Only run system tests. {SYS_TESTS_DISAMBIGUATION}",
    )

    test_subset.add_argument(
        "-p", "--pattern", help="File name pattern to match\n" "Default is test*.py"
    )

    args, unknown = parser.parse_known_args()

    pytest_args = []

    # Pre-pytest version of run_ctsm_py_tests suppressed warnings by default
    if not args.no_disable_warnings:
        pytest_args += ["--disable-warnings"]

    # Handle -v/--verbose
    if args.verbose:
        pytest_args += ["--verbose"]

    # Handle --debug with vs. without arg. Note that --debug is no longer mutually exclusive with
    # --verbose.
    if args.debug is not None:
        # Old run_ctsm_py_tests behavior: Run tests with maximum verbosity
        if args.debug == "":
            pytest_args += ["--verbose", "--capture=no"]
        # pylint's --debug
        else:
            pytest_args += ["--debug", args.debug]

    # Warn user about ambiguous -s
    if args.sys:
        logger.info("Running system tests only. %s", SYS_TESTS_DISAMBIGUATION)

    # Pass any unknown args to pytest directly
    pytest_args += unknown

    return args, pytest_args


def _get_pytest_help():
    # Get pytest help text
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            pytest.main(["--help"])
        except SystemExit:
            # pytest may call sys.exit(); ignore it
            pass
    pytest_help_text = buf.getvalue()

    # Extract help for options we care about
    debug_help = get_pytest_help_item(pytest_help_text, "--debug")
    return debug_help


def get_pytest_help_item(pytest_help_text: str, option: str) -> str:
    """
    Extract a single pytest CLI option (e.g. "--debug") from the full ``pytest --help`` output and
    collapse its multi-line description into a single line.

    The function locates the first line starting with ``option`` and then collects all immediately
    following lines that are more indented than that header line. These indented lines are treated
    as the option's description. The first line break is replaced with a colon and the description
    lines are joined with single spaces.

    Parameters
    ----------
    pytest_help_text : str
        The full text output of ``pytest --help``.
    option : str
        A substring identifying the option to extract
        (e.g. "--debug" or "--override-ini").

    Returns
    -------
    str
        A single-line string of the form:
        "<option header>: <collapsed description>"

    Raises
    ------
    RuntimeError
        If no line starting with ``option`` is found in the help text.

    Notes
    -----
    This function relies on pytest's help formatting convention that description lines are indented
    more deeply than their corresponding option header line.
    """
    lines = pytest_help_text.splitlines()

    # Create regex pattern to match the option at the start of a line (after whitespace).
    # Use word boundary \b to ensure we match exact options (e.g., "--pdb" won't match "--pdbcls").
    option_pattern = re.compile(r"^\s*" + re.escape(option) + r"\b")

    for i, line in enumerate(lines):
        if option_pattern.match(line):
            # Found line starting with the option

            # Split the line into option header and description.
            # Description typically starts after 2+ consecutive spaces.
            match = re.match(r"^(\s*\S+(?:\s+\S+)*?)\s{2,}(.*)$", line)

            if match:
                # Line has both option header and description
                option_header = match.group(1).strip()
                first_desc = match.group(2)
            else:
                # No description on this line (or only single spaces)
                option_header = line.strip()
                first_desc = ""

            # Calculate the indentation level of the first matched ("header") line.
            # This will be used to identify continuation lines.
            header_indent = len(line) - len(line.lstrip())

            # Collect all description text, starting with any on the header line
            description_lines = []
            if first_desc:
                description_lines.append(first_desc)

            # Collect all continuation lines that are more indented than the header
            for next_line in lines[i + 1 :]:
                # Stop if we hit a blank line
                if not next_line.strip():
                    break

                # Calculate indentation of this line
                indent = len(next_line) - len(next_line.lstrip())

                # Stop if indentation is not deeper than header (indicates we've reached the next
                # option or section)
                if indent <= header_indent:
                    break

                # This is a continuation line - add it to the description
                description_lines.append(next_line.strip())

            # Join all description lines into a single string with spaces
            description = " ".join(description_lines)

            # Return formatted result: "option_header: description"
            return f"{option_header}: {description}"

    # If we get here, the option was not found
    raise RuntimeError(f"Failed to get pytest help for {option}")
