"""
A class to handle longitude values, plus helper functions for that class
"""

import logging
from argparse import ArgumentTypeError
import numpy as np

logger = logging.getLogger(__name__)


def _check_lon_type_180(lon_in):
    """
    Checks value range of longitude with type 180
    """
    lon_min = np.min(lon_in)
    lon_max = np.max(lon_in)
    for lon in [lon_min, lon_max]:
        if not -180 <= lon <= 180:
            raise ValueError("(All values of) lon_in must be in the range [-180, 180]")


def _check_lon_type_360(lon_in):
    """
    Checks value range of longitude with type 360
    """
    lon_min = np.min(lon_in)
    lon_max = np.max(lon_in)
    for lon in [lon_min, lon_max]:
        if not 0 <= lon <= 360:
            raise ValueError("(All values of) lon_in must be in the range [0, 360]")


def _check_lon_value_given_type(lon_in, lon_type_in):
    """
    Checks value range of longitude
    """
    if lon_type_in == 180:
        _check_lon_type_180(lon_in)
    elif lon_type_in == 360:
        _check_lon_type_360(lon_in)
    else:
        raise RuntimeError(f"You must specify lon_type_in as 180 or 360, not {lon_type_in}")


def _convert_lon_type_180_to_360(lon_in):
    """
    Convert a longitude from type 180 to type 360
    """
    _check_lon_type_180(lon_in)
    lon_out = lon_in % 360
    logger.info(
        "Converting longitude from [-180, 180] to [0, 360]: %s to %s",
        str(lon_in),
        str(lon_out),
    )

    return lon_out


def detect_lon_type(lon_in):
    """
    Detect longitude type of a given numeric. If lon_in contains more than one number (as in a list
    or Numpy array), this function will assume all members are of the same type if (a) there is at
    least one unambiguous member and (b) all unambiguous members are of the same type.
    """
    lon_min = np.min(lon_in)
    lon_max = np.max(lon_in)
    if lon_min < -180:
        raise ValueError(f"(Minimum) longitude < -180: {lon_min}")
    if lon_max > 360:
        raise ValueError(f"(Maximum) longitude > 360: {lon_max}")
    min_type_180 = lon_min < 0
    max_type_360 = lon_max > 180
    if min_type_180 and max_type_360:
        raise RuntimeError("Longitude array contains values of both types 180 and 360")
    if not min_type_180 and not max_type_360:
        raise ArgumentTypeError("Longitude(s) ambiguous; could be type 180 or 360")
    if min_type_180:
        lon_type = 180
    else:
        lon_type = 360
    return lon_type


def _convert_lon_type_360_to_180(lon_in):
    """
    Convert a longitude from type 360 to type 180
    """
    _check_lon_type_360(lon_in)

    if lon_in <= 180:
        lon_out = lon_in
    else:
        lon_out = lon_in - 360

    logger.info(
        "Converting longitude from [-180, 180] to [0, 360]: %s to %s",
        str(lon_in),
        str(lon_out),
    )

    return lon_out


def check_other_is_lontype(other):
    """
    Used in comparison operators to throw an error if the "other" object being compared isn't also
    a Longitude object. This makes it so that comparing longitudes requires that both sides of the
    comparison must be Longitude objects and thus must have a specified longitude type (180 or 360).

    We could try to coerce non-Longitude `other` to Longitude, but that might result in
    situations where tests think everything works but code will fail if `other` is
    ambiguous.
    """
    if not isinstance(other, Longitude):
        raise TypeError(
            f"Comparison not supported between instances of 'Longitude' and '{type(other)}'"
        )


class Longitude:
    """
    A class to keep track of a longitude and its type
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, lon, lon_type):

        # Do not allow nesting Longitude objects
        if isinstance(lon, Longitude):
            raise TypeError("Trying to initialize a Longitude object with a Longitude object")

        # If lon_type comes in as a string, convert it to int
        if isinstance(lon_type, str):
            lon_type = int(lon_type)

        # Check that longitude is within bounds
        _check_lon_value_given_type(lon, lon_type)

        self._lon = lon
        self._lon_type = lon_type

    def _check_lons_same_type(self, other):
        """
        If you're comparing two Longitudes of different types in different hemispheres, then
        `lon1 > lon2` and `lon2 < lon1` will incorrectly give different answers! We could make it so
        that this doesn't fail as long as symmetricality isn't violated, but that might lead to
        unexpected failures in practice.
        """
        if self.lon_type() != other.lon_type():
            raise TypeError("Comparison not supported between Longitudes of different types")

    # __eq__ makes it so that == and != both work.
    def __eq__(self, other):
        check_other_is_lontype(other)
        return self._lon == other.get(self._lon_type)

    def __lt__(self, other):
        check_other_is_lontype(other)
        self._check_lons_same_type(other)
        return self._lon < other._lon

    def __gt__(self, other):
        check_other_is_lontype(other)
        self._check_lons_same_type(other)
        return self._lon > other._lon

    def __le__(self, other):
        check_other_is_lontype(other)
        self._check_lons_same_type(other)
        return self._lon <= other._lon

    def __ge__(self, other):
        check_other_is_lontype(other)
        self._check_lons_same_type(other)
        return self._lon >= other._lon

    def __str__(self):
        """
        We don't allow implicit string conversion because the user should always specify the
        Longitude type they want
        """
        raise NotImplementedError("Use Longitude.get_str() instead of implicit string conversion")

    def get(self, lon_type_out):
        """
        Get the longitude value, converting longitude type if needed
        """
        if lon_type_out == self._lon_type:
            return self._lon
        if lon_type_out == 180:
            return _convert_lon_type_360_to_180(self._lon)
        if lon_type_out == 360:
            return _convert_lon_type_180_to_360(self._lon)
        raise RuntimeError(f"Add handling for lon_type_out {lon_type_out}")

    def get_str(self, lon_type_out):
        """
        Get the longitude value as a string, converting longitude type if needed
        """
        lon_out = self.get(lon_type_out)
        # Use float() because the standard in CTSM filenames is to put .0 after whole-number values
        return str(float(lon_out))

    def lon_type(self):
        """
        Getter method for self._lon_type
        """
        return self._lon_type
