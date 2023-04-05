#!/usr/bin/env python3
"""
Unit tests for mesh_type

You can run this by:
    python -m unittest test_unit_mesh_type.py
"""

import os
import sys
import numpy as np
import xarray as xr
import unittest

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.mesh_type import MeshType
from ctsm.site_and_regional.mesh_plot_type import MeshPlotType

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
        lons = xr.DataArray(lon0, name="lon", dims=x_dim, coords={x_dim: lon0})
        lats = xr.DataArray(lat0, name="lat", dims=y_dim, coords={y_dim: lat0})

        self.mesh = MeshPlotType(lats, lons)

    def test_read_file_fails_notXarrayDataset(self):
        """Test that read_file properly fails if not given an X-array Dataset"""
        with self.assertRaisesRegex(SystemExit, "Input file is not a X-Array DataSet type"):
            self.mesh.read_file( 1.0 )
