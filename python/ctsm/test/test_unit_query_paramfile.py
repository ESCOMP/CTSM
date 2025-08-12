#!/usr/bin/env python3

"""Unit tests for query_paramfile"""

import unittest
import sys

from ctsm import unit_testing

from ctsm.param_utils import query_paramfile as qp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestQueryParamfile(unittest.TestCase):
    """Tests of query_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self.orig_argv

    def test_query_paramfile_args_short(self):
        """Test that all arguments can be set correctly with shortnames"""
        input_path = "/path/to/input.nc"
        sys.argv = ["get_arguments", "-i", input_path, "-p", "pft1,pft2", "var1,var2"]
        args = qp.get_arguments()
        self.assertEqual(input_path, args.input)
        self.assertEqual("pft1,pft2", args.pft)
        self.assertEqual("var1,var2", args.variables)

    def test_query_paramfile_args_long(self):
        """Test that all arguments can be set correctly with longnames"""
        input_path = "/path/to/input.nc"
        sys.argv = ["get_arguments", "--input", input_path, "--pft", "pft1,pft2", "var1,var2"]
        args = qp.get_arguments()
        self.assertEqual(input_path, args.input)
        self.assertEqual("pft1,pft2", args.pft)
        self.assertEqual("var1,var2", args.variables)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
