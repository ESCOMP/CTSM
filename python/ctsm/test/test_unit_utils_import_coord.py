#!/usr/bin/env python3

"""
Unit tests for utils.py functions related to importing coordinate variables
"""

import unittest
import os
import sys
import shutil

import tempfile
import xarray as xr
import numpy as np

# -- add python/ctsm  to path (needed if we want to run test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.path_utils import path_to_ctsm_root
from ctsm.ctsm_pylib_dependent_utils import import_coord_1d, import_coord_2d

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


# Allow as many public methods as needed...
# pylint: disable=too-many-public-methods
# Allow all the instance attributes that we need
# pylint: disable=too-many-instance-attributes
class TestUtilsImportCoord(unittest.TestCase):
    """
    Tests the importcoord* subroutines from utils.py
    """

    def setUp(self):
        """Setup for trying out the methods"""
        self._previous_dir = os.getcwd()
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._testinputs_path = testinputs_path
        self._tempdir = tempfile.mkdtemp()

        self._1d_lonlat_file = os.path.join(
            self._testinputs_path, "cropcals", "swh_rf_ggcmi_crop_calendar_phase3_v1.01.nc4"
        )
        self._2d_lonlat_file = os.path.join(
            self._testinputs_path,
            "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031_modified.nc",
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_importcoord1d(self):
        """
        Tests importing a 1-d lat/lon variable
        """
        ds = xr.open_dataset(self._1d_lonlat_file)
        lat, n_lat = import_coord_1d(ds, "lat")
        np.testing.assert_equal(n_lat, 360)
        np.testing.assert_array_equal(lat.values[:4], [89.75, 89.25, 88.75, 88.25])
        np.testing.assert_array_equal(lat.values[-4:], [-88.25, -88.75, -89.25, -89.75])

    def test_importcoord1d_attrs(self):
        """
        Tests attributes of an imported 1-d lat/lon variable
        """
        ds = xr.open_dataset(self._1d_lonlat_file)
        lat, _ = import_coord_1d(ds, "lat")
        # Unlike import_coord_2d, import_coord_1d doesn't rename the long name.
        expected_attributes = {
            "long_name": ds["lat"].attrs["long_name"],
            "units": "degrees_north",
        }
        self.assertDictEqual(lat.attrs, expected_attributes)

    def test_importcoord1d_too_many_dims(self):
        """
        Tests that 1d-importing function errors when given a 2d variable to import
        """
        ds = xr.open_dataset(self._2d_lonlat_file)
        with self.assertRaises(
            SystemExit,
            msg="Expected 1 dimension for LATIXY; found 2: ('lsmlat', 'lsmlon')",
        ):
            import_coord_1d(ds, "LATIXY")

    def test_importcoord2d(self):
        """
        Tests importing a 2-d lat/lon variable
        """
        ds = xr.open_dataset(self._2d_lonlat_file)
        lat, _ = import_coord_2d(ds, "lat", "LATIXY")
        expected_values = np.array([-13.9, -11.7, -9.5, -7.3, -5.1]).astype(np.float32)
        np.testing.assert_array_equal(lat.values, expected_values)

    def test_importcoord2d_attrs(self):
        """
        Tests attributes of an imported 2-d lat/lon variable
        """
        ds = xr.open_dataset(self._2d_lonlat_file)
        lat, _ = import_coord_2d(ds, "lat", "LATIXY")
        expected_attributes = {
            "long_name": "coordinate latitude",
            "units": "degrees_north",
        }
        self.assertDictEqual(lat.attrs, expected_attributes)

    def test_importcoord2d_rename_dim(self):
        """
        Tests renaming of an imported 2-d lat/lon variable
        """
        ds = xr.open_dataset(self._2d_lonlat_file)
        lat, _ = import_coord_2d(ds, "lat", "LATIXY")
        self.assertTupleEqual(lat.dims, ("lat",))

    def test_importcoord2d_no_dim_contains_coordName(self):
        """
        Tests that 2d-importing function errors when given a nonexistent dim name
        """
        ds = xr.open_dataset(self._2d_lonlat_file)
        ds = ds.rename({"lsmlat": "abc"})
        with self.assertRaises(
            SystemExit,
            msg="ERROR: Expected 1 dimension name containing lat; found 0: []",
        ):
            import_coord_2d(ds, "lat", "LATIXY")

    def test_importcoord2d_1_dim_containing(self):
        """
        Tests that 2d-importing function errors when given an ambiguous dim name
        """
        ds = xr.open_dataset(self._2d_lonlat_file)
        ds = ds.rename({"lsmlon": "lsmlat2"})
        with self.assertRaises(
            SystemExit,
            msg="Expected 1 dimension name containing lat; found 2: ['lsmlat', 'lsmlat2']",
        ):
            import_coord_2d(ds, "lat", "LATIXY")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
