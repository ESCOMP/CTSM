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

    def test_region(self):
        """Do a basic test for a small regional grid"""
        self._infile = os.path.join(
            self._testinputs_path,
            "surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214_modified_with_crop.nc"
        )
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "LATIXY",
            "--lon",
            "LONGXY",
            "--no-plot",
            "--mask",
            "PFTDATA_MASK",
            "--output",
            self.mesh_out,
        ]
        main()

    def test_noplot_add_mask(self):
        """Do a simple basic test without plotting and also adding mask"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--no-plot",
            "--mask",
            "PFTDATA_MASK",
            "--output",
            self.mesh_out,
        ]
        main()

    def test_add_area(self):
        """Do a simple basic test without plotting and also adding area"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--area",
            "AREA",
            "--output",
            self.mesh_out,
        ]
        main()

    def test_noinput(self):
        """Test with an input file that does not exist"""
        sys.argv = [
            "mesh_maker",
            "--input",
            "zztop",
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(SystemExit, "Input file not found."):
            main()

    def test_nolongs(self):
        """Bad name for longitude"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "zztop",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(SystemExit, "Input file does not have variable named zztop"):
            main()

    def test_nolats(self):
        """Bad name for latitude"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "zztop",
            "--lon",
            "lsmlon",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(SystemExit, "Input file does not have variable named zztop"):
            main()

    def test_badareaname(self):
        """Bad name for area"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--area",
            "zztop",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(
            SystemExit, "Input file does not have area variable named zztop"
        ):
            main()

    def test_badmaskname(self):
        """Bad name for mask"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--mask",
            "zztop",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(
            SystemExit, "Input file does not have mask variable named zztop"
        ):
            main()

    def test_badareaunits(self):
        """Bad area units"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--area",
            "PCT_CROP",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(
            SystemExit, r"Area does NOT have the correct units of radians\^2 but has unitless"
        ):
            main()

    def test_badmaskvalues(self):
        """Bad mask values"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "lsmlat",
            "--lon",
            "lsmlon",
            "--mask",
            "LAKEDEPTH",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(SystemExit, "Mask variable is not within 0 to 1"):
            main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
