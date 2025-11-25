#!/usr/bin/env python3

"""
Unit tests for _get_rectangle
"""

import unittest

import numpy as np
import xarray as xr

from ctsm import unit_testing
from ctsm.longitude import Longitude
from ctsm.modify_input_files.modify_fsurdat import ModifyFsurdat

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access

## Too many instant variables as part of the class (too many self.<varible> in the SetUp)
# pylint: disable=too-many-instance-attributes


class TestModifyFsurdat(unittest.TestCase):
    """Tests the setvar_lev functions and the
    _get_rectangle, check_varlist, and set_varlist methods
    """

    def setUp(self):
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        self.min_lon = 2  # expects min_lon < max_lon
        self.min_lat = 3  # expects min_lat < max_lat
        self.lon_type = 360
        longxy, latixy, self.cols, self.rows = self._get_longxy_latixy(
            _min_lon=self.min_lon,
            _max_lon=10,
            _min_lat=self.min_lat,
            _max_lat=12,
            lon_type=self.lon_type,
        )

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        self.lon_1 = Longitude(3, self.lon_type)
        self.lon_2 = Longitude(5, self.lon_type)  # lon_1 < lon_2
        self.lat_1 = 5
        self.lat_2 = 7  # lat_1 < lat_2

        # create xarray dataset containing lev0, lev1, and lev2 variables;
        # the fsurdat_modify tool reads variables like this from fsurdat file
        var_1d = np.arange(self.cols)
        var_lev0 = var_1d * np.ones((self.rows, self.cols))
        var_lev1 = var_1d * np.ones((self.cols, self.rows, self.cols))
        var_lev2 = var_1d * np.ones((self.rows, self.cols, self.rows, self.cols))
        var_lev3 = var_1d * np.ones((self.cols, self.rows, self.cols, self.rows, self.cols))
        my_data = xr.Dataset(
            data_vars={
                "LONGXY": (["x", "y"], longxy),  # use LONGXY as var_lev0
                "LATIXY": (["x", "y"], latixy),  # __init__ expects LONGXY, LATIXY
                "urbdens": (["numurbl"], var_1d),  # numurbl needs to be dimension
                "var_lev0": (["x", "y"], var_lev0),
                "var_lev1": (["w", "x", "y"], var_lev1),
                "var_lev2": (["v", "w", "x", "y"], var_lev2),
                "var_lev3": (["z", "v", "w", "x", "y"], var_lev3),
                "VAR_LEV0_UPPERCASE": (["x", "y"], var_lev0),
                "VAR_LEV1_UPPERCASE": (["w", "x", "y"], var_lev1),
                "VAR_LEV2_UPPERCASE": (["v", "w", "x", "y"], var_lev2),
                "VAR_LEV3_UPPERCASE": (["z", "v", "w", "x", "y"], var_lev3),
            }
        )

        # create ModifyFsurdat object
        self.modify_fsurdat = ModifyFsurdat(
            my_data=my_data,
            lon_1=self.lon_1,
            lon_2=self.lon_2,
            lat_1=self.lat_1,
            lat_2=self.lat_2,
            landmask_file=None,
            lat_dimname=None,
            lon_dimname=None,
        )

    def test_setvarLev(self):
        """
        Tests that setvar_lev0, setvar_lev1, and setvar_lev2 update values of
        variables within a rectangle defined by user-specified
        lon_1, lon_2, lat_1, lat_2
        """

        # initialize and then modify the comparison matrices
        comp_lev0 = self.modify_fsurdat.file.LONGXY
        comp_lev1 = self.modify_fsurdat.file.var_lev1
        comp_lev2 = self.modify_fsurdat.file.var_lev2
        val_for_rectangle = 1.5
        comp_lev0[
            self.lat_1 - self.min_lat : self.lat_2 - self.min_lat + 1,
            self.lon_1.get(self.lon_type)
            - self.min_lon : self.lon_2.get(self.lon_type)
            - self.min_lon
            + 1,
        ] = val_for_rectangle
        comp_lev1[
            ...,
            self.lat_1 - self.min_lat : self.lat_2 - self.min_lat + 1,
            self.lon_1.get(self.lon_type)
            - self.min_lon : self.lon_2.get(self.lon_type)
            - self.min_lon
            + 1,
        ] = val_for_rectangle
        comp_lev2[
            ...,
            self.lat_1 - self.min_lat : self.lat_2 - self.min_lat + 1,
            self.lon_1.get(self.lon_type)
            - self.min_lon : self.lon_2.get(self.lon_type)
            - self.min_lon
            + 1,
        ] = val_for_rectangle

        # test setvar
        self.modify_fsurdat.setvar_lev0("LONGXY", val_for_rectangle)
        np.testing.assert_array_equal(self.modify_fsurdat.file.LONGXY, comp_lev0)

        self.modify_fsurdat.setvar_lev1("var_lev1", val_for_rectangle, self.cols - 1)
        np.testing.assert_array_equal(self.modify_fsurdat.file.var_lev1, comp_lev1)

        self.modify_fsurdat.setvar_lev2("var_lev2", val_for_rectangle, self.cols - 1, self.rows - 1)
        np.testing.assert_array_equal(self.modify_fsurdat.file.var_lev2, comp_lev2)

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
            _min_lon=min_lon, _max_lon=7, _min_lat=min_lat, _max_lat=8, lon_type=self.lon_type
        )

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        lon_1 = Longitude(3, self.lon_type)
        lon_2 = Longitude(5, self.lon_type)  # lon_1 < lon_2
        lat_1 = 6
        lat_2 = 8  # lat_1 < lat_2
        rectangle = ModifyFsurdat._get_rectangle(
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            longxy=longxy,
            latixy=latixy,
        )
        not_rectangle = np.logical_not(rectangle)
        compare = np.ones((rows, cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[
            lat_1 - min_lat : lat_2 - min_lat + 1,
            lon_1.get(self.lon_type) - min_lon : lon_2.get(self.lon_type) - min_lon + 1,
        ] = 0
        np.testing.assert_array_equal(not_rectangle, compare)

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
            _min_lon=min_lon, _max_lon=6, _min_lat=min_lat, _max_lat=5, lon_type=self.lon_type
        )

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = Longitude(0, self.lon_type)
        lon_2 = Longitude(4, self.lon_type)  # lon_1 < lon_2
        lat_1 = 4
        lat_2 = 0  # lat_1 > lat_2
        rectangle = ModifyFsurdat._get_rectangle(
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            longxy=longxy,
            latixy=latixy,
        )
        not_rectangle = np.logical_not(rectangle)
        compare = np.ones((rows, cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[
            : lat_2 - min_lat + 1,
            lon_1.get(self.lon_type) - min_lon : lon_2.get(self.lon_type) - min_lon + 1,
        ] = 0
        compare[
            lat_1 - min_lat :,
            lon_1.get(self.lon_type) - min_lon : lon_2.get(self.lon_type) - min_lon + 1,
        ] = 0
        np.testing.assert_array_equal(not_rectangle, compare)

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
            _min_lon=min_lon, _max_lon=359, _min_lat=min_lat, _max_lat=90, lon_type=self.lon_type
        )

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = Longitude(4, self.lon_type)
        lon_2 = Longitude(2, self.lon_type)  # lon_1 > lon_2
        lat_1 = 2
        lat_2 = 3  # lat_1 < lat_2
        rectangle = ModifyFsurdat._get_rectangle(
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            longxy=longxy,
            latixy=latixy,
        )
        not_rectangle = np.logical_not(rectangle)
        compare = np.ones((rows, cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[lat_1 - min_lat : lat_2 - min_lat + 1, : lon_2.get(self.lon_type) - min_lon + 1] = 0
        compare[lat_1 - min_lat : lat_2 - min_lat + 1, lon_1.get(self.lon_type) - min_lon :] = 0
        np.testing.assert_array_equal(not_rectangle, compare)

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
            _min_lon=min_lon, _max_lon=5, _min_lat=min_lat, _max_lat=6, lon_type=180
        )

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_type = 180
        lon_1 = Longitude(-1, lon_type=lon_type)
        lon_2 = Longitude(-6, lon_type=lon_type)  # lon_1 > lon_2
        lat_1 = 0
        lat_2 = -3  # lat_1 > lat_2
        rectangle = ModifyFsurdat._get_rectangle(
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            longxy=longxy,
            latixy=latixy,
        )
        not_rectangle = np.logical_not(rectangle)
        compare = np.ones((rows, cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[: lat_2 - min_lat + 1, : lon_2.get(lon_type) - min_lon + 1] = 0
        compare[: lat_2 - min_lat + 1, lon_1.get(lon_type) - min_lon :] = 0
        compare[lat_1 - min_lat :, : lon_2.get(lon_type) - min_lon + 1] = 0
        compare[lat_1 - min_lat :, lon_1.get(lon_type) - min_lon :] = 0
        np.testing.assert_array_equal(not_rectangle, compare)

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
            _min_lon=min_lon, _max_lon=359, _min_lat=min_lat, _max_lat=5, lon_type=self.lon_type
        )

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = Longitude(355, self.lon_type)
        lon_2 = Longitude(5, self.lon_type)  # lon_1 > lon_2
        lat_1 = -4
        lat_2 = -6  # lat_1 > lat_2
        rectangle = ModifyFsurdat._get_rectangle(
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            longxy=longxy,
            latixy=latixy,
        )
        not_rectangle = np.logical_not(rectangle)
        compare = np.ones((rows, cols))
        # assert this to confirm intuitive understanding of these matrices
        self.assertEqual(np.size(not_rectangle), np.size(compare))

        # Hardwire where I expect not_rectangle to be False (0)
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple
        compare[: lat_2 - min_lat + 1, : lon_2.get(self.lon_type) - min_lon + 1] = 0
        compare[: lat_2 - min_lat + 1, lon_1.get(self.lon_type) - min_lon :] = 0
        compare[lat_1 - min_lat :, : lon_2.get(self.lon_type) - min_lon + 1] = 0
        compare[lat_1 - min_lat :, lon_1.get(self.lon_type) - min_lon :] = 0
        np.testing.assert_array_equal(not_rectangle, compare)

    def test_getNotRectangle_latsOutOfBounds(self):
        """
        Tests that out-of-bound latitude values abort with message
        Out-of-bound longitudes already tested in test_unit_utils.py
        """
        # get longxy, latixy that would normally come from an fsurdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = 0  # expects min_lon < max_lon
        min_lat = -5  # expects min_lat < max_lat
        longxy, latixy, _, _ = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=359, _min_lat=min_lat, _max_lat=5, lon_type=self.lon_type
        )

        # get not_rectangle from user-defined lon_1, lon_2, lat_1, lat_2
        # I have chosen the lon/lat ranges to match their corresponding index
        # values to keep this simple (see usage below)
        lon_1 = Longitude(355, self.lon_type)
        lon_2 = Longitude(5, self.lon_type)
        lat_1 = -91
        lat_2 = 91
        with self.assertRaisesRegex(
            SystemExit, "lat_1 and lat_2 need to be in the range -90 to 90"
        ):
            _ = ModifyFsurdat._get_rectangle(
                lon_1=lon_1,
                lon_2=lon_2,
                lat_1=lat_1,
                lat_2=lat_2,
                longxy=longxy,
                latixy=latixy,
            )

    def test_check_varlist_lists(self):
        """Test the check_varlist method for list for dimensions that works"""
        lev1list = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        lev2list = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        settings = {"var_lev1": lev1list, "var_lev2": lev2list}
        settings_new = self.modify_fsurdat.check_varlist(settings)
        self.assertEqual(
            settings_new, settings, "list of variable settings not identical as expected"
        )

    def test_check_varlist_lists_wrongsizes(self):
        """Test the check_varlist method for lists to gracefully fail when the sizes are wrong"""
        lev1list = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        settings = {"var_lev1": lev1list}
        with self.assertRaisesRegex(
            SystemExit,
            " Variable var_lev1 is 8 is of the wrong size."
            + " It should be = 9 in input settings dictionary",
        ):
            self.modify_fsurdat.check_varlist(settings)

    def test_get_numurb_dens(self):
        """Check that get num urban density types is correct"""
        self.assertEqual(
            self.modify_fsurdat.get_urb_dens(),
            9,
            "Default number of urban density types is correct",
        )

    def test_check_varlist_uppercase(self):
        """Test the check_varlist method for all the dimensions that
        works with allowuppercase option"""
        vallist = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        vallist2 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        expected = {
            "VAR_LEV0_UPPERCASE": 100.0,
            "VAR_LEV1_UPPERCASE": vallist,
            "VAR_LEV2_UPPERCASE": vallist2,
        }
        settings = {
            "var_lev0_uppercase": 100.0,
            "var_lev1_uppercase": vallist,
            "var_lev2_uppercase": vallist2,
        }
        settings_new = self.modify_fsurdat.check_varlist(settings, allow_uppercase_vars=True)
        self.assertEqual(
            expected,
            settings_new,
            "list of variable settings not converted to uppercase as expected",
        )

    def test_check_varlist_badvar(self):
        """Test the check_varlist method for a variable not on the file"""
        settings = {"badvar": 100.0}
        with self.assertRaisesRegex(
            SystemExit, "Variable badvar is NOT in the input settings dictionary"
        ):
            self.modify_fsurdat.check_varlist(settings)

    def test_check_varlist_badvar_uppercase(self):
        """Test the check_varlist method for a variable not on the file with allow uppercase"""
        settings = {"badvar": 100.0}
        with self.assertRaisesRegex(
            SystemExit, "Variable BADVAR is NOT in the input settings dictionary"
        ):
            self.modify_fsurdat.check_varlist(settings, allow_uppercase_vars=True)

    def test_set_varlist_toohighdim(self):
        """Test the set_varlist method for a variable of too high a dimension"""
        settings = {"var_lev3": 100.0}
        with self.assertRaisesRegex(
            SystemExit, "Variable var_lev3 is a higher dimension than currently allowed"
        ):
            self.modify_fsurdat.set_varlist(settings)

    def test_set_varlist_toohighdim_uppercase(self):
        """Test the set_varlist method for a variable of too high a dimension in uppercase"""
        settings = {"var_lev3_uppercase": 100.0}
        with self.assertRaisesRegex(
            SystemExit,
            "For higher dimensional vars, the variable needs to be expressed as a "
            + "list of values of the dimension size = 9 for variable=VAR_LEV3_UPPERCASE",
        ):
            self.modify_fsurdat.check_varlist(settings, allow_uppercase_vars=True)

    def _get_longxy_latixy(self, _min_lon, _max_lon, _min_lat, _max_lat, lon_type):
        """
        Return longxy, latixy, cols, rows
        """
        cols = _max_lon - _min_lon + 1
        rows = _max_lat - _min_lat + 1

        long = np.arange(_min_lon, _max_lon + 1)
        if lon_type == 180:
            long = [Longitude(longitude, lon_type).get(360) for longitude in long]
        longxy = long * np.ones((rows, cols))
        compare = np.repeat([long], rows, axis=0)  # alternative way to form
        # assert this to confirm intuitive understanding of these matrices
        np.testing.assert_array_equal(longxy, compare)

        lati = np.arange(_min_lat, _max_lat + 1)
        self.assertEqual(min(lati), _min_lat)
        self.assertEqual(max(lati), _max_lat)
        latixy_transp = lati * np.ones((cols, rows))
        compare = np.repeat([lati], cols, axis=0)  # alternative way to form
        # assert this to confirm intuitive understanding of these matrices
        np.testing.assert_array_equal(latixy_transp, compare)
        latixy = np.transpose(latixy_transp)

        return longxy, latixy, cols, rows


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
