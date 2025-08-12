#!/usr/bin/env python3

"""Unit tests for set_paramfile"""

import unittest
import os
import sys

from ctsm import unit_testing

from ctsm.param_utils import set_paramfile as sp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


class TestUnitSetParamfile(unittest.TestCase):
    """Unit tests of set_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self.orig_argv

    def test_set_paramfile_args_short(self):
        """Test that all arguments can be set correctly with shortnames"""
        output_path = "/path/to/output.nc"
        sys.argv = [
            "get_arguments",
            "-i",
            PARAMFILE,
            "-p",
            "pft1,pft2",
            "-o",
            output_path,
            "-v",
            "var1,var2",
        ]
        args = sp.get_arguments()
        self.assertEqual(PARAMFILE, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(["var1", "var2"], args.variables)
        self.assertEqual(output_path, args.output)

    def test_set_paramfile_args_long(self):
        """Test that all arguments can be set correctly with longnames"""
        output_path = "/path/to/output.nc"
        sys.argv = [
            "get_arguments",
            "--input",
            PARAMFILE,
            "--pft",
            "pft1,pft2",
            "--output",
            output_path,
            "--variables",
            "var1,var2",
        ]
        args = sp.get_arguments()
        self.assertEqual(PARAMFILE, args.input)
        self.assertEqual(["pft1", "pft2"], args.pft)
        self.assertEqual(["var1", "var2"], args.variables)
        self.assertEqual(output_path, args.output)

    def test_set_paramfile_error_missing_input(self):
        """Test that it errors if input file doesn't exist"""
        output_path = "/path/to/output.nc"
        sys.argv = [
            "get_arguments",
            "--input",
            "nwuefweirbdfdiurbe",
            "--pft",
            "pft1,pft2",
            "--output",
            output_path,
        ]
        with self.assertRaises(FileNotFoundError):
            sp.get_arguments()

    def test_set_paramfile_error_existing_output(self):
        """Test that it errors if output file already exists"""
        sys.argv = [
            "get_arguments",
            "--input",
            PARAMFILE,
            "--pft",
            "pft1,pft2",
            "--output",
            PARAMFILE,
        ]
        with self.assertRaises(FileExistsError):
            sp.get_arguments()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
