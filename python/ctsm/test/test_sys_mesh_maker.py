#!/usr/bin/env python3
"""
System tests for mesh_maker

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
from ctsm.mesh_maker import read_main

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
            testinputs_path, "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031_modified.nc"
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
            "LATIXY",
            "--lon",
            "LONGXY",
            "--output",
            self.mesh_out,
        ]
        main()

    def test_region(self):
        """Do a basic test for a small regional grid"""
        infile = os.path.join(
            self._testinputs_path,
            "surfdata_5x5_amazon_hist_78pfts_CMIP6_2000_c230517_modified_with_crop.nc",
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
            "LANDFRAC_PFT",
            "--output",
            self.mesh_out,
        ]
        main()

    def compare_mesh_files(self, mesh_out, expected):
        """Compare two mesh files that you expect to be equal"""
        self.assertEqual(
            mesh_out.dims["coordDim"], expected.dims["coordDim"], "coordDim not the same"
        )
        if "origGridRank" in mesh_out.variables and "origGridRank" in expected.variables:
            self.assertEqual(
                mesh_out.dims["origGridRank"],
                expected.dims["origGridRank"],
                "origGridRank not the same",
            )
            equalorigGridDims = mesh_out.origGridDims == expected.origGridDims
            self.assertTrue(equalorigGridDims.all, "origGridDims different")
            compare_files = True
        else:
            # don't compare files if origGridRank isn't on one of the files
            compare_files = False

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
        equalelementConn = mesh_out.elementConn == expected.elementConn
        equalnumElementConn = mesh_out.numElementConn == expected.numElementConn
        equalcenterCoords = mesh_out.centerCoords == expected.centerCoords
        equalelementMask = mesh_out.elementMask == expected.elementMask

        if "elementArea" in mesh_out.variables and "elementArea" in expected.variables:
            equalelementArea = mesh_out.elementArea == expected.elementArea
            self.assertTrue(equalelementArea.all, "area different")

        self.assertTrue(equalelementConn.all, "elementConn different")
        self.assertTrue(equalnumElementConn.all, "numElementConn different")
        self.assertTrue(equalcenterCoords.all, "centerCoords different")
        self.assertTrue(equalelementMask.all, "mask different")
        if compare_files:
            self.assertTrue(
                mesh_out.equals(expected),
                "Output mesh does not compare to the expected baseline file",
            )

    def test_domainfile_region_warea(self):
        """
        Do a basic test for a small regional grid with a domain file
        rather than a surfdata file including area
        """
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

    def test_domainfile_SAmerica_region_warea(self):
        """
        Do a basic test for a South America regional grid with a domain file
        rather than a surfdata file including area
        """
        infile = os.path.join(
            self._testinputs_path, "domain.lnd.fv0.9x1.25_gx1v7_f09_58x45_SouthAmerica_c230522.nc"
        )
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
            self._testinputs_path,
            "ESMF_mesh_fv0.9x1.25_gx1v7_f09_58x45_SouthAmerica_from_domain_c230522.nc",
        )
        mesh_out = xr.open_dataset(self.mesh_out)
        expected = xr.open_dataset(expected_mesh)
        self.compare_mesh_files(mesh_out, expected)

    def test_domainfile_f10_warea(self):
        """
        Do a test of converting the f10 domain file
        """
        infile = os.path.join(self._testinputs_path, "domain.lnd.fv10x15_gx3v7.180321.nc")
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
        with self.assertRaisesRegex(SystemExit, "Expected size for element connections is wrong"):
            main()

    def test_readfile(self):
        """
        Test that reading a file results in the same mesh as converting one
        """
        infile = os.path.join(
            self._testinputs_path, "ESMF_mesh_5x5pt_amazon_from_domain_c230308.nc"
        )
        sys.argv = [
            "mesh_maker",
            "--input",
            infile,
            "--lat",
            "yc",
            "--lon",
            "xc",
            "--output",
            self.mesh_out,
        ]
        read_main()
        mesh_out = xr.open_dataset(self.mesh_out)
        expected = xr.open_dataset(infile)
        self.compare_mesh_files(mesh_out, expected)

    def test_add_mask(self):
        """Do a simple basic test also adding mask"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "LATIXY",
            "--lon",
            "LONGXY",
            "--mask",
            "LANDFRAC_PFT",
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
            "LATIXY",
            "--lon",
            "LONGXY",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(SystemExit, "Input file not found."):
            main()

    def test_singlepoint_dies(self):
        """Test that a single point file dies because we don't need mesh files
        for single point cases"""
        infile = os.path.join(
            self._testinputs_path,
            "surfdata_1x1_mexicocityMEX_hist_16pfts_Irrig_CMIP6_simyr2000_c221206.nc",
        )
        sys.argv = [
            "mesh_maker",
            "--input",
            infile,
            "--lat",
            "LATIXY",
            "--lon",
            "LONGXY",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(
            SystemExit, r"No need to create a mesh file for a single point grid."
        ):
            main()

    def test_nolongs(self):
        """Bad name for longitude"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "LATIXY",
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
            "LONGXY",
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
            "LATIXY",
            "--lon",
            "LONGXY",
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
            "LATIXY",
            "--lon",
            "LONGXY",
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
            "LATIXY",
            "--lon",
            "LONGXY",
            "--area",
            "PCT_CROP",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(
            SystemExit, r"Area does NOT have the correct units of radians\^2 but has unitless"
        ):
            main()

    def test_missingreaunits(self):
        """Missing area units"""
        self._infile = os.path.join(
            self._testinputs_path,
            "surfdata_5x5_amazon_hist_78pfts_CMIP6_2000_c230517_modified_with_crop.nc",
        )
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "LATIXY",
            "--lon",
            "LONGXY",
            "--area",
            "PCT_CROP",
            "--output",
            self.mesh_out,
        ]
        with self.assertRaisesRegex(SystemExit, r"Units attribute is NOT on the area variable"):
            main()

    def test_badmaskvalues(self):
        """Bad mask values"""
        sys.argv = [
            "mesh_maker",
            "--input",
            self._infile,
            "--lat",
            "LATIXY",
            "--lon",
            "LONGXY",
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
