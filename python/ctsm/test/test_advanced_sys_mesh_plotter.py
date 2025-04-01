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
        testinputs_path = os.path.join(
            path_to_ctsm_root(),
            "python",
            "ctsm",
            "test",
            "testinputs",
        )
        self._testinputs_path = testinputs_path
        self._infile = os.path.join(
            testinputs_path,
            "ESMF_mesh_5x5pt_amazon_from_domain_c230308.nc",
        )
        self._tempdir = tempfile.mkdtemp()
        self.mesh_out = os.path.join(self._tempdir, "mesh_out")
        self.test_basic_argv = [
            "mesh_plotter",
            "--input",
            self._infile,
            "--output",
            self.mesh_out,
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_basic(self):
        """Do a simple basic test"""
        sys.argv = self.test_basic_argv
        main()
        plotfiles = glob.glob(os.path.join(self._tempdir, "*.png"))
        if not plotfiles:
            self.fail("plot files were NOT created as they should have")

    def test_dpi(self):
        """Test setting dpi"""
        sys.argv = self.test_basic_argv + [
            "--dpi",
            "198.7",
        ]
        main()
        plotfiles = glob.glob(os.path.join(self._tempdir, "*.png"))
        if not plotfiles:
            self.fail("plot files were NOT created as they should have")

    def test_need_overwrite(self):
        """Ensure failure if output file exists but --overwrite not given"""
        sys.argv = self.test_basic_argv
        main()
        with self.assertRaisesRegex(
            FileExistsError, "File already exists but --overwrite not given"
        ):
            main()

    def test_outdir(self):
        """Test that --outdir option works"""
        outdir = os.path.join(self._tempdir, "abc123")
        sys.argv = [
            "mesh_plotter",
            "--input",
            self._infile,
            "--outdir",
            outdir,
        ]
        main()
        plotfiles = glob.glob(os.path.join(outdir, "*.png"))
        if not plotfiles:
            self.fail("plot files were NOT created as they should have")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
