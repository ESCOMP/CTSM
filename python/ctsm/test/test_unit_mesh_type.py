#!/usr/bin/env python3
"""
Unit tests for mesh_type

You can run this by:
    python -m unittest test_unit_mesh_type.py
"""

import os
import sys
import unittest
import numpy as np
import xarray as xr

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing

# from ctsm.site_and_regional.mesh_type import MeshType
from ctsm.site_and_regional.mesh_type import MeshType

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestMeshType(unittest.TestCase):
    """
    Basic class for testing mesh_type.py.
    """

    def setUp(self):
        """Setup for all tests"""
        lon0 = np.array([120.0])
        lat0 = np.array([45.0])
        x_dim = "lon"
        y_dim = "lat"
        self.lons = xr.DataArray(lon0, name="lon", dims=x_dim, coords={x_dim: lon0})
        self.lats = xr.DataArray(lat0, name="lat", dims=y_dim, coords={y_dim: lat0})

        self.mesh = MeshType(self.lats, self.lons)

    def test_read_file_fails_notXarrayDataset(self):
        """Test that read_file properly fails if not given an X-array Dataset"""
        with self.assertRaisesRegex(SystemExit, "Input file is not a X-Array DataSet type"):
            self.mesh.read_file(1.0)

    def test_read_file_fails_badvarlist(self):
        """Test that read_file properly fails if input dataset does not have required variables"""
        ds = xr.Dataset(
            {
                "PCT_NATVEG": xr.DataArray(
                    data=np.random.rand(1, 1),
                    dims=["lat", "lon"],
                    coords={"lat": self.lats, "lon": self.lons},
                    attrs={
                        "long_name": "total percent natural vegetation landunit",
                        "units": "unitless",
                    },
                ),
            },
            attrs={"Conventions": "test data only"},
        )
        with self.assertRaisesRegex(
            SystemExit, "Variable expected to be on an ESMF mesh file is NOT found: nodeCoords"
        ):
            self.mesh.read_file(ds)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
