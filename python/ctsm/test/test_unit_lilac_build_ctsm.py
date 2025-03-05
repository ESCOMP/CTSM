#!/usr/bin/env python3

"""Unit tests for lilac_build_ctsm"""

import unittest
from unittest.mock import patch
from io import StringIO

from ctsm import unit_testing
from ctsm.lilac_build_ctsm import _commandline_args, _check_and_transform_os

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

# pylint: disable=line-too-long


class TestBuildCtsm(unittest.TestCase):
    """Tests of lilac_build_ctsm"""

    # ------------------------------------------------------------------------
    # Tests of _commandline_args
    # ------------------------------------------------------------------------

    def test_commandlineArgs_rebuild_valid(self):
        """Test _commandline_args with --rebuild, with a valid argument list (no disallowed args)"""
        _ = _commandline_args(args_to_parse=["build/directory", "--rebuild"])

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid1(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is required for non-rebuilds, without a dash
        """
        expected_re = r"--compiler cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=["build/directory", "--rebuild", "--compiler", "intel"]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid2(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is required for new machines, without a dash
        """
        expected_re = r"--os cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(args_to_parse=["build/directory", "--rebuild", "--os", "linux"])
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid3(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is required for new machines, with a dash
        """
        expected_re = r"--netcdf-path cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=[
                    "build/directory",
                    "--rebuild",
                    "--netcdf-path",
                    "/path/to/netcdf",
                ]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid4(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is optional for new non-rebuild
        that isn't None
        """
        expected_re = r"--no-build cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(args_to_parse=["build/directory", "--rebuild", "--no-build"])
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_rebuild_invalid5(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is invalid with this option

        This tests an argument that is optional for new machines, which also has a default
        that isn't None
        """
        expected_re = r"--gmake cannot be provided if --rebuild is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=["build/directory", "--rebuild", "--gmake", "mymake"]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    def test_commandlineArgs_noRebuild_valid(self):
        """Test _commandline_args without --rebuild or --machine, with a valid argument list

        (all required things present)
        """
        _ = _commandline_args(
            args_to_parse=[
                "build/directory",
                "--os",
                "linux",
                "--compiler",
                "intel",
                "--netcdf-path",
                "/path/to/netcdf",
                "--esmf-mkfile-path",
                "/path/to/esmf/lib/esmf.mk",
                "--max-mpitasks-per-node",
                "16",
                "--no-pnetcdf",
            ]
        )

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_noRebuild_invalid1(self, mock_stderr):
        """Test _commandline_args without --rebuild or --machine, with a missing required argument

        This tests an argument in the non-rebuild-required list
        """
        expected_re = r"--compiler must be provided if --rebuild is not set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=[
                    "build/directory",
                    "--os",
                    "linux",
                    "--netcdf-path",
                    "/path/to/netcdf",
                    "--esmf-mkfile-path",
                    "/path/to/esmf/lib/esmf.mk",
                    "--max-mpitasks-per-node",
                    "16",
                    "--no-pnetcdf",
                ]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_noRebuild_invalid2(self, mock_stderr):
        """Test _commandline_args without --rebuild or --machine, with a missing required argument

        This tests an argument in the new-machine-required list
        """
        expected_re = r"--os must be provided if neither --machine nor --rebuild are set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=[
                    "build/directory",
                    "--compiler",
                    "intel",
                    "--netcdf-path",
                    "/path/to/netcdf",
                    "--esmf-mkfile-path",
                    "/path/to/esmf/lib/esmf.mk",
                    "--max-mpitasks-per-node",
                    "16",
                    "--no-pnetcdf",
                ]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_noRebuild_invalidNeedToDictatePnetcdf(self, mock_stderr):
        """Test _commandline_args without --rebuild or --machine: invalid b/c nothing specified for pnetcdf"""
        expected_re = (
            r"For a user-defined machine, need to specify either --no-pnetcdf or --pnetcdf-path"
        )
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=[
                    "build/directory",
                    "--os",
                    "linux",
                    "--compiler",
                    "intel",
                    "--netcdf-path",
                    "/path/to/netcdf",
                    "--esmf-mkfile-path",
                    "/path/to/esmf/lib/esmf.mk",
                    "--max-mpitasks-per-node",
                    "16",
                ]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_noRebuild_invalidConflictingPnetcdf(self, mock_stderr):
        """Test _commandline_args without --rebuild or --machine: invalid b/c of conflicting specifications for pnetcdf"""
        expected_re = r"--no-pnetcdf cannot be given if you set --pnetcdf-path"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=[
                    "build/directory",
                    "--os",
                    "linux",
                    "--compiler",
                    "intel",
                    "--netcdf-path",
                    "/path/to/netcdf",
                    "--esmf-mkfile-path",
                    "/path/to/esmf/lib/esmf.mk",
                    "--max-mpitasks-per-node",
                    "16",
                    "--no-pnetcdf",
                    "--pnetcdf-path",
                    "/path/to/pnetcdf",
                ]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    def test_commandlineArgs_machine_valid(self):
        """Test _commandline_args with --machine, with a valid argument list

        (all required things present)
        """
        _ = _commandline_args(
            args_to_parse=[
                "build/directory",
                "--machine",
                "mymachine",
                "--compiler",
                "intel",
            ]
        )

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_machine_missingRequired(self, mock_stderr):
        """Test _commandline_args with --machine, with a missing required argument"""
        expected_re = r"--compiler must be provided if --rebuild is not set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(args_to_parse=["build/directory", "--machine", "mymachine"])
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_machine_illegalArg1(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is illegal with this option

        This tests an argument that is required for new machines
        """
        expected_re = r"--os cannot be provided if --machine is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=[
                    "build/directory",
                    "--machine",
                    "mymachine",
                    "--compiler",
                    "intel",
                    "--os",
                    "linux",
                ]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    @patch("sys.stderr", new_callable=StringIO)
    def test_commandlineArgs_machine_illegalArg2(self, mock_stderr):
        """Test _commandline_args with --rebuild, with an argument that is illegal with this option

        This tests an argument that is optional for new machines
        """
        expected_re = r"--gmake cannot be provided if --machine is set"
        with self.assertRaises(SystemExit):
            _ = _commandline_args(
                args_to_parse=[
                    "build/directory",
                    "--machine",
                    "mymachine",
                    "--compiler",
                    "intel",
                    "--gmake",
                    "mymake",
                ]
            )
        self.assertRegex(mock_stderr.getvalue(), expected_re)

    # ------------------------------------------------------------------------
    # Tests of _check_and_transform_os
    # ------------------------------------------------------------------------

    def test_checkAndTransformOs_valid(self):
        """Test _check_and_transform_os with valid input"""
        os = _check_and_transform_os("linux")
        self.assertEqual(os, "LINUX")

    def test_checkAndTransformOs_invalid(self):
        """Test _check_and_transform_os with invalid input"""
        with self.assertRaises(ValueError):
            _ = _check_and_transform_os("bad_os")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
