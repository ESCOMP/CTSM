#!/usr/bin/env python3
"""
Unit tests for mesh_plotter

You can run this by:
    python -m unittest test_advanced_unit_mesh_plotter.py
"""

import unittest
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.mesh_plotter import get_parser, process_and_check_args, main

# pylint: disable=invalid-name


class TestMeshMaker(unittest.TestCase):
    """
    Basic class for testing mesh_plotter.py.
    """

    def test_input_file_dne(self):
        """
        Test that exits if input file does not exist
        """
        sys.argv = [
            "mesh_plotter",
            "--input",
            "zztop.nc",
            "--outdir",
            ".",
        ]
        with self.assertRaisesRegex(SystemExit, "Input file not found."):
            main()

    def test_default_outfile_as_expected(self):
        """
        Test that the default outfile is as expected
        """
        infile = "ctsm/test/testinputs/default_data.cfg"
        sys.argv = [
            "mesh_plotter",
            "--input",
            infile,
        ]
        parser = get_parser()
        args = parser.parse_args()
        args = process_and_check_args(args)
        expected_outdir = os.path.join(os.getcwd(), "meshes")
        self.assertEqual(args.out_dir, expected_outdir, "Default out_dir is not as expected")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
