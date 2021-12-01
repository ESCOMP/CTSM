#!/usr/bin/env python3

import unittest

import numpy as np

from ctsm import unit_testing
from ctsm.modify_fsurdat.modify_fsurdat import ModifyFsurdat

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable; in pylint can use: disable=invalid-name

class TestFsurdatModifier(unittest.TestCase):

    def test_getNotRectangle_lon1leLon2Lat1leLat2(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 <= lon_2 and lat_1 <= lat_2, expect not_rectangle to be
        False in a rectangle bounded by these lon/lat values
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # get cols, rows also
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=1, _max_lon=5, _min_lat=1, _max_lat=6)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = 2
        lon_2 = 4  # lon_1 < lon_2
        lat_1 = 2
        lat_2 = 3  # lat_1 < lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[lat_1-1:lat_2, lon_1-1:lon_2] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_lon1leLon2Lat1gtLat2(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 <= lon_2 and lat_1 > lat_2, expect not_rectangle to be
        False in two rectangles bounded by these lon/lat values, one to the
        north of lat_1 and one to the south of lat_2
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # get cols, rows also
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=1, _max_lon=5, _min_lat=1, _max_lat=6)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = 2
        lon_2 = 4  # lon_1 < lon_2
        lat_1 = 3
        lat_2 = 2  # lat_1 > lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[:lat_1, lon_1-1:lon_2] = 0
        compare[lat_2-1:, lon_1-1:lon_2] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_lon1gtLon2Lat1leLat2(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 > lon_2 and lat_1 <= lat_2, expect not_rectangle to be
        False in two rectangles bounded by these lon/lat values, one to the
        east of lat_1 and one to the west of lat_2
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # get cols, rows also
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=1, _max_lon=5, _min_lat=1, _max_lat=6)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        # TODO add out-of-bounds test
        #      does test still work with lon_1 = 355 and lon_2 = 5?
        lon_1 = 4
        lon_2 = 2  # lon_1 > lon_2
        lat_1 = 2
        lat_2 = 3  # lat_1 < lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[lat_1-1:lat_2, :lon_2] = 0
        compare[lat_1-1:lat_2, lon_1-1:] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_lon1gtLon2Lat1gtLat2(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 > lon_2 and lat_1 > lat_2, expect not_rectangle to be
        False in four rectangles bounded by these lon/lat values, in the
        top left, top right, bottom left, and bottom right of the domain
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # get cols, rows also
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=1, _max_lon=5, _min_lat=1, _max_lat=6)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = 4
        lon_2 = 2  # lon_1 > lon_2
        lat_1 = 3
        lat_2 = 2  # lat_1 > lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[:lat_1, :lon_2] = 0
        compare[:lat_1, lon_1-1:] = 0
        compare[lat_2-1:, :lon_2] = 0
        compare[lat_2-1:, lon_1-1:] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def _get_longxy_latixy(self, _min_lon, _max_lon, _min_lat, _max_lat):
        """
        Return longxy, latixy, cols, rows
        """
        cols = _max_lon - _min_lon + 1
        rows = _max_lat - _min_lat + 1

        long = np.arange(_min_lon, _max_lon + 1)
        self.assertEqual(min(long), _min_lon)
        self.assertEqual(max(long), _max_lon)
        longxy = long * np.ones((rows,cols))
        compare = np.repeat([long], rows, axis=0)
        # assert this to confirm intuitive understanding of these matrices
        self.assertIsNone(np.testing.assert_array_equal(longxy, compare))

        lati = np.arange(_min_lat, _max_lat + 1)
        self.assertEqual(min(lati), _min_lat)
        self.assertEqual(max(lati), _max_lat)
        latixy_transp = lati * np.ones((cols,rows))
        compare = np.repeat([lati], cols, axis=0)
        # assert this to confirm intuitive understanding of these matrices
        self.assertIsNone(np.testing.assert_array_equal(latixy_transp, compare))
        latixy = np.transpose(latixy_transp)

        return longxy, latixy, cols, rows

"""Unit tests for _get_not_rectangle
"""
if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
