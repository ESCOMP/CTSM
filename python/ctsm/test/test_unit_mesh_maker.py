#!/usr/bin/env python3
"""
Unit tests for mesh_maker

You can run this by:
    python -m unittest test_unit_mesh_maker.py
"""

import unittest
import argparse
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.mesh_maker import get_parser, main

# pylint: disable=invalid-name


class TestMeshMaker(unittest.TestCase):
    """
    Basic class for testing mesh_maker.py.
    """

    def test_input_file_dne(self):
        """
        Test that exits if input file does not exist
        """
        sys.argv = ["mesh_maker", "--input", "zztop.nc", "--lat", "lsmlat", "--lon", "lsmlon"]
        with self.assertRaisesRegex(SystemExit, "Input file not found."):
            main()

    def test_outfile_and_outdir(self):
        """
        Test that exits if both outfile and outdir are provided
        """
        infile = "ctsm/test/testinputs/default_data.cfg"
        sys.argv = ["mesh_maker", "--input", infile, "--lat", "lsmlat", "--lon", "lsmlon", "--outdir", ".", "--output", "outthing.nc"]
        with self.assertRaisesRegex(argparse.ArgumentError, "You have provided both --outdir and --output."):
            main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
