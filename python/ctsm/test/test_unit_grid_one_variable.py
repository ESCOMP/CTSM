#!/usr/bin/env python3

"""
Unit tests for grid_one_variable
"""

import unittest

import numpy as np
import xarray as xr

from ctsm import unit_testing
from ctsm.crop_calendars import grid_one_variable as g1v

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access

## Too many instant variables as part of the class (too many self.<varible> in the SetUp)
# pylint: disable=too-many-instance-attributes


class TestCreateFilledArray(unittest.TestCase):
    """Unit tests for create_filled_array"""

    def setUp(self):
        # Set up this_ds, which will provide us with sizes of dimensions in most cases
        lat_vals = [55.0, 56.0, 57.0]
        lat_da = xr.DataArray(
            data=lat_vals,
            dims=["lat"],
            coords={"lat": lat_vals},
        )
        lon_vals = [255.0, 256.0, 257.0]
        lon_da = xr.DataArray(
            data=lon_vals,
            dims=["lon"],
            coords={"lon": lon_vals},
        )
        self.this_ds = xr.Dataset(
            data_vars={
                "lat": lat_da,
                "lon": lon_da,
            }
        )

    def test_create_filled_array_fillNone(self):
        """
        Test create_filled_array() with fill_value None: Should be filled with NaN
        """

        fill_value = None
        thisvar_da_dummy = xr.DataArray()
        new_dims = ["lat", "lon"]

        result = g1v.create_filled_array(self.this_ds, fill_value, thisvar_da_dummy, new_dims)

        self.assertTrue(np.all(np.isnan(result)))

    def test_create_filled_array_fill1(self):
        """
        Test create_filled_array() with fill_value 1: Should be filled with 1
        """

        fill_value = 1.0
        thisvar_da_dummy = xr.DataArray()
        new_dims = ["lat", "lon"]

        result = g1v.create_filled_array(self.this_ds, fill_value, thisvar_da_dummy, new_dims)

        self.assertTrue(np.all(result == fill_value))

    def test_create_filled_array_fill0(self):
        """
        Test create_filled_array() with fill_value 0: Should be filled with 0
        """

        fill_value = 0.0
        thisvar_da_dummy = xr.DataArray()
        new_dims = ["lat", "lon"]

        result = g1v.create_filled_array(self.this_ds, fill_value, thisvar_da_dummy, new_dims)

        self.assertTrue(np.all(result == fill_value))
