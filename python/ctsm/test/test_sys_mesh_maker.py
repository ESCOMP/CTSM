#!/usr/bin/env python3
"""
System tests for mesh_maker

"""

import unittest
import os
import sys
import tempfile
import shutil

# pylint: disable=wrong-import-position
from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.mesh_maker import main

# pylint: disable=invalid-name


class SysTestMeshMaker(unittest.TestCase):
    """
    Basic class for testing mesh_maker.py.
    """

    def setUp(self):
        """Setup for all tests"""
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._testinputs_path = testinputs_path
        self._infile = os.path.join(
            testinputs_path,
            "surfdata_1x1_mexicocityMEX_hist_16pfts_Irrig_CMIP6_simyr2000_c221206_modified.nc",
        )
        self._tempdir = tempfile.mkdtemp()
        self.mesh_out = os.path.join(self._tempdir, "mesh_out.nc")

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_basic(self):
        """Do a simple basic test"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--output",
            self.mesh_out,
        ]
        main()

    def test_noplot_add_area_mask(self):
        """Do a simple basic test without plotting and also adding area and mask"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--no-plot",
            "--area",
            "AREA",
            "--mask",
            "PFTDATA_MASK",
            "--output",
            self.mesh_out,
        ]
        main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
