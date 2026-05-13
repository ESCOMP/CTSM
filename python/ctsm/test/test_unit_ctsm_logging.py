#!/usr/bin/env python3

"""Unit tests for functions in ctsm_logging"""

import unittest
import io
from contextlib import redirect_stdout

from ctsm import unit_testing
from ctsm.ctsm_logging import log, error
from ctsm.utils import datetime_string

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


DATETIME_STR_PATTERN = r"\d{4}\-\d{2}\-\d{2} \d{2}:\d{2}:\d{2}"


class TestLog(unittest.TestCase):
    """
    Tests of log() function

    Currently not testing ability to write to file, because unittest makes that difficult. So just
    testing write to stdout.
    """

    def test_datetime_str_pattern(self):
        """Test that our regex matches output of datetime_str"""
        self.assertRegex(datetime_string(), expected_regex=DATETIME_STR_PATTERN)

    def test_log_without_logger(self):
        """
        Tests the log() function without providing a logger for writing to file
        """
        msg = "abc123"
        f = io.StringIO()
        with redirect_stdout(f):
            log(None, msg)

        # Check that stdout matches what we expect
        stdout = f.getvalue()
        expected_regex = DATETIME_STR_PATTERN + r"\s+test_log_without_logger\s+" + msg
        self.assertRegex(stdout, expected_regex=expected_regex)


class TestError(unittest.TestCase):
    """
    Tests of error() function

    Currently not testing ability to write to file, because unittest makes that difficult. So just
    testing write to stdout and error raising.
    """

    def test_error_without_logger(self):
        """
        Tests the error() function without providing a logger for writing to file and without
        specifying a custom error type
        """
        msg = "abc123"
        f = io.StringIO()
        error_raised = None
        try:
            with redirect_stdout(f):
                error(None, msg)
        except Exception as e:  # pylint: disable=broad-exception-caught
            error_raised = e

        # Check that stdout matches what we expect
        stdout = f.getvalue()
        expected_regex = DATETIME_STR_PATTERN + r"\s+test_error_without_logger\s+" + msg
        self.assertRegex(stdout, expected_regex=expected_regex)

        # Check that error is correct
        self.assertFalse(error_raised is None)
        self.assertIsInstance(error_raised, RuntimeError)
        self.assertEqual(msg, str(error_raised))

    def test_error_without_logger_custom_err(self):
        """
        Tests the error() function without providing a logger for writing to file and
        specifying a custom error type
        """
        msg = "abc123"
        f = io.StringIO()
        error_raised = None
        error_type = ValueError
        try:
            with redirect_stdout(f):
                error(None, msg, error_type=error_type)
        except Exception as e:  # pylint: disable=broad-exception-caught
            error_raised = e

        # Check that stdout matches what we expect
        stdout = f.getvalue()
        expected_regex = DATETIME_STR_PATTERN + r"\s+test_error_without_logger_custom_err\s+" + msg
        self.assertRegex(stdout, expected_regex=expected_regex)

        # Check that error is correct
        self.assertFalse(error_raised is None)
        self.assertIsInstance(error_raised, error_type)
        self.assertEqual(msg, str(error_raised))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
