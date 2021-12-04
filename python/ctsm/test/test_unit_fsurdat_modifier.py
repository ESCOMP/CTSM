#!/usr/bin/env python3

import unittest

import numpy as np

from ctsm import unit_testing
from ctsm.utils import lon_range_0_to_360
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
        Work with integer lon/lat values to keep the testing simple
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = 2  # expects min_lon < max_lon
        min_lat = 3  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=7, _min_lat=min_lat, _max_lat=8)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        lon_1 = 3
        lon_2 = 5  # lon_1 < lon_2
        lat_1 = 6
        lat_2 = 8  # lat_1 < lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[lat_1-min_lat:lat_2-min_lat+1, lon_1-min_lon:lon_2-min_lon+1] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_lon1leLon2Lat1gtLat2(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 <= lon_2 and lat_1 > lat_2, expect not_rectangle to be
        False in two rectangles bounded by these lon/lat values, one to the
        north of lat_1 and one to the south of lat_2
        Work with integer lon/lat values to keep the testing simple
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = -3  # expects min_lon < max_lon
        min_lat = -2  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=6, _min_lat=min_lat, _max_lat=5)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = 0
        lon_2 = 4  # lon_1 < lon_2
        lat_1 = 4
        lat_2 = 0  # lat_1 > lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[:lat_2-min_lat+1, lon_1-min_lon:lon_2-min_lon+1] = 0
        compare[lat_1-min_lat:, lon_1-min_lon:lon_2-min_lon+1] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_lon1gtLon2Lat1leLat2(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 > lon_2 and lat_1 <= lat_2, expect not_rectangle to be
        False in two rectangles bounded by these lon/lat values, one to the
        east of lat_1 and one to the west of lat_2
        Work with integer lon/lat values to keep the testing simple
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = 1  # expects min_lon < max_lon
        min_lat = 1  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=359, _min_lat=min_lat, _max_lat=90)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
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
        compare[lat_1-min_lat:lat_2-min_lat+1, :lon_2-min_lon+1] = 0
        compare[lat_1-min_lat:lat_2-min_lat+1, lon_1-min_lon:] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_lon1gtLon2Lat1gtLat2(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 > lon_2 and lat_1 > lat_2, expect not_rectangle to be
        False in four rectangles bounded by these lon/lat values, in the
        top left, top right, bottom left, and bottom right of the domain
        Work with integer lon/lat values to keep the testing simple
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = -8  # expects min_lon < max_lon
        min_lat = -9  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=5, _min_lat=min_lat, _max_lat=6)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = -1
        lon_2 = -6  # lon_1 > lon_2
        lat_1 = 0
        lat_2 = -3  # lat_1 > lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[:lat_2-min_lat+1, :lon_2-min_lon+1] = 0
        compare[:lat_2-min_lat+1, lon_1-min_lon:] = 0
        compare[lat_1-min_lat:, :lon_2-min_lon+1] = 0
        compare[lat_1-min_lat:, lon_1-min_lon:] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_lonsStraddle0deg(self):
        """
        Tests that not_rectangle is True and False in the grid cells expected
        according to the user-specified lon_1, lon_2, lat_1, lat_2
        When lon_1 > lon_2 and lat_1 > lat_2, expect not_rectangle to be
        False in four rectangles bounded by these lon/lat values, in the
        top left, top right, bottom left, and bottom right of the domain
        Work with integer lon/lat values to keep the testing simple
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = 0  # expects min_lon < max_lon
        min_lat = -5  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=359, _min_lat=min_lat, _max_lat=5)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = 355
        lon_2 = 5  # lon_1 > lon_2
        lat_1 = -4
        lat_2 = -6  # lat_1 > lat_2
        not_rectangle = ModifyFsurdat._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
            longxy=longxy, latixy=latixy)
        compare = np.ones((rows,cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[:lat_2-min_lat+1, :lon_2-min_lon+1] = 0
        compare[:lat_2-min_lat+1, lon_1-min_lon:] = 0
        compare[lat_1-min_lat:, :lon_2-min_lon+1] = 0
        compare[lat_1-min_lat:, lon_1-min_lon:] = 0
        self.assertIsNone(np.testing.assert_array_equal(not_rectangle, compare))

    def test_getNotRectangle_latsOutOfBound(self):
        """
        Tests that out-of-bound latitude values abort with message
        Out-of-bound longitudes already tested in test_unit_utils.py
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = 0  # expects min_lon < max_lon
        min_lat = -5  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=359, _min_lat=min_lat, _max_lat=5)

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = 355
        lon_2 = 5
        lat_1 = -91
        lat_2 = 91
        with self.assertRaisesRegex(SystemExit,
            "lat_1 and lat_2 need to be in the range -90 to 90"):
            _ = ModifyFsurdat._get_not_rectangle(
                lon_1=lon_1, lon_2=lon_2, lat_1=lat_1, lat_2=lat_2,
                longxy=longxy, latixy=latixy)

    def _get_longxy_latixy(self, _min_lon, _max_lon, _min_lat, _max_lat):
        """
        Return longxy, latixy, cols, rows
        """
        cols = _max_lon - _min_lon + 1
        rows = _max_lat - _min_lat + 1

        long = np.arange(_min_lon, _max_lon + 1)
        long = [lon_range_0_to_360(longitude) for longitude in long]
        longxy = long * np.ones((rows,cols))
        compare = np.repeat([long], rows, axis=0)  # alternative way to form
        # assert this to confirm intuitive understanding of these matrices
        self.assertIsNone(np.testing.assert_array_equal(longxy, compare))

        lati = np.arange(_min_lat, _max_lat + 1)
        self.assertEqual(min(lati), _min_lat)
        self.assertEqual(max(lati), _max_lat)
        latixy_transp = lati * np.ones((cols,rows))
        compare = np.repeat([lati], cols, axis=0)  # alternative way to form
        # assert this to confirm intuitive understanding of these matrices
        self.assertIsNone(np.testing.assert_array_equal(latixy_transp, compare))
        latixy = np.transpose(latixy_transp)

        return longxy, latixy, cols, rows

"""Unit tests for _get_not_rectangle
"""
if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
