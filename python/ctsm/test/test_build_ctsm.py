#!/usr/bin/env python

"""Unit tests for build_ctsm
"""

import unittest
from unittest.mock import patch
from io import StringIO

from ctsm import unit_testing
from ctsm.build_ctsm import _commandline_args, _check_and_transform_os

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

class TestBuildCtsm(unittest.TestCase):
    """Tests of build_ctsm"""

    def test_commandlineArgs_rebuild_valid(self):
        """Test _commandline_args with --rebuild, with a valid argument list (no disallowed args)"""
        # pylint: disable=no-self-use
        _ = _commandline_args(args_to_parse=['build/directory', '--rebuild'])

    @patch('sys.stderr', new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid1(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is required for non-rebuilds, without a dash
        """
        expected_re = r"--compiler cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(args_to_parse=['build/directory',
                                                 '--rebuild',
                                                 '--compiler', 'intel'])
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch('sys.stderr', new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid2(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is required for non-rebuilds, with a dash
        """
        expected_re = r"--netcdf-path cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(args_to_parse=['build/directory',
                                                 '--rebuild',
                                                 '--netcdf-path', '/path/to/netcdf'])
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch('sys.stderr', new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid3(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is optional for non-rebuilds, which also has a default
        that isn't None
        """
        expected_re = r"--gmake cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(args_to_parse=['build/directory',
                                                 '--rebuild',
                                                 '--gmake', 'mymake'])
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    def test_commandlineArgs_noRebuild_valid(self):
        """Test _commandline_args without --rebuild, with a valid argument list

        (all required things present)
        """
        # pylint: disable=no-self-use
        _ = _commandline_args(args_to_parse=['build/directory',
                                             '--os', 'linux',
                                             '--compiler', 'intel',
                                             '--netcdf-path', '/path/to/netcdf',
                                             '--esmf-lib-path', '/path/to/esmf/lib'])

    @patch('sys.stderr', new_callable=StringIO)
    def test_commandlineArgs_noRebuild_invalid(self, mock_stderr):
        """Test _commandline_args without --rebuild, with a missing required argument"""
        expected_re = r"--compiler must be provided if --rebuild is not set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(args_to_parse=['build/directory',
                                                 '--os', 'linux',
                                                 '--netcdf-path', '/path/to/netcdf',
                                                 '--esmf-lib-path', '/path/to/esmf/lib'])
        self.assertRegex(mock_stderr.getvalue(), expected_re)


    def test_checkAndTransformOs_valid(self):
        """Test _check_and_transform_os with valid input"""
        os = _check_and_transform_os('linux')
        self.assertEqual(os, 'LINUX')

    def test_checkAndTransformOs_invalid(self):
        """Test _check_and_transform_os with invalid input"""
        with self.assertRaises(ValueError):
            _ = _check_and_transform_os('bad_os')

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
