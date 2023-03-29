#!/usr/bin/env python3
"""
Advanced System tests for mesh_maker (requires the ctsm_pylib_wdask conda environment)

"""

import unittest
import os
import sys
import tempfile
import shutil
import xarray as xr

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
        infile = os.path.join(
            self._testinputs_path,
            "surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214_modified_with_crop.nc",
        )
        sys.argv = [
            "mesh_maker",
            "--input",
            infile,
            "--lat",
            "LATIXY",
            "--lon",
            "LONGXY",
            "--mask",
            "PFTDATA_MASK",
            "--output",
            self.mesh_out,
        ]
        main()

    def compare_mesh_files(self, mesh_out, expected):
        """Compare two mesh files that you expect to be equal"""
        self.assertEqual(
            mesh_out.dims["coordDim"], expected.dims["coordDim"], "coordDim not the same"
        )
        self.assertEqual(
            mesh_out.dims["origGridRank"],
            expected.dims["origGridRank"],
            "origGridRank not the same",
        )
        self.assertEqual(
            mesh_out.dims["nodeCount"], expected.dims["nodeCount"], "nodeCount not the same"
        )
        self.assertEqual(
            mesh_out.dims["elementCount"],
            expected.dims["elementCount"],
            "elementCount not the same",
        )
        self.assertEqual(
            mesh_out.dims["maxNodePElement"],
            expected.dims["maxNodePElement"],
            "maxNodePElement not the same",
        )
        equalorigGridDims = mesh_out.origGridDims == expected.origGridDims
        equalelementConn = mesh_out.elementConn == expected.elementConn
        equalnumElementConn = mesh_out.numElementConn == expected.numElementConn
        equalcenterCoords = mesh_out.centerCoords == expected.centerCoords
        equalelementMask = mesh_out.elementMask == expected.elementMask
        equalelementArea = mesh_out.elementArea == expected.elementArea
        self.assertTrue(equalorigGridDims.all, "origGridDims different")
        self.assertTrue(equalelementConn.all, "elementConn different")
        self.assertTrue(equalnumElementConn.all, "numElementConn different")
        self.assertTrue(equalcenterCoords.all, "centerCoords different")
        self.assertTrue(equalelementMask.all, "mask different")
        self.assertTrue(equalelementArea.all, "area different")
        print(mesh_out)
        print(expected)
        self.assertTrue(
            mesh_out.equals(expected), "Output mesh does not compare to the expected baseline file"
        )

    def test_domainfile_region_warea(self):
        """Do a basic test for a small regional grid with a domain file rather than a surfdata file including area"""
        infile = os.path.join(self._testinputs_path, "domain.lnd.5x5pt-amazon_navy.090715.nc")
        sys.argv = [
            "mesh_maker",
            "--input",
            infile,
            "--lat",
            "yc",
            "--lon",
            "xc",
            "--mask",
            "mask",
            "--area",
            "area",
            "--output",
            self.mesh_out,
        ]
        main()
        expected_mesh = os.path.join(
            self._testinputs_path, "ESMF_mesh_5x5pt_amazon_from_domain_c230308.nc"
        )
        mesh_out = xr.open_dataset(self.mesh_out)
        expected = xr.open_dataset(expected_mesh)
        self.compare_mesh_files(mesh_out, expected)

if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
