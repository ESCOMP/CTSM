#!/usr/bin/env python3
"""
Advanced System tests for mesh_plotter (requires the ctsm_pylib_wdask conda environment)

"""

import unittest
import os
import sys
import tempfile
import shutil
import glob

# pylint: disable=wrong-import-position
from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.mesh_plotter import main

# pylint: disable=invalid-name


class SysTestMeshMaker(unittest.TestCase):
    """
    Basic class for testing mesh_plotter.py.
    """

    def setUp(self):
        """Setup for all tests"""
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._testinputs_path = testinputs_path
        self._infile = os.path.join(
            testinputs_path,
            "ESMF_mesh_5x5pt_amazon_from_domain_c230308.nc",
        )
        self._tempdir = tempfile.mkdtemp()
        self.mesh_out = self._tempdir + "/mesh_out"

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_basic(self):
        """Do a simple basic test"""
        sys.argv = [
            "mesh_plotter",
            "--input",
            self._infile,
            "--output",
            self.mesh_out,
        ]
        main()
        plotfiles = glob.glob(self._tempdir + "/*.png")
        if not plotfiles:
            self.fail("plot files were NOT created as they should have")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
