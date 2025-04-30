"""
A class to handle longitude values, plus helper functions for that class
"""

import logging
from argparse import ArgumentTypeError
from functools import total_ordering

logger = logging.getLogger(__name__)


def _check_lon_type_180(lon_in):
    """
    Checks value range of longitude with type 180
    """
    if not -180 <= lon_in <= 180:
        raise ValueError(f"lon_in needs to be in the range [-180, 180]: {lon_in}")


def _check_lon_type_360(lon_in):
    """
    Checks value range of longitude with type 360
    """
    if not 0 <= lon_in <= 360:
        raise ValueError(f"lon_in needs to be in the range [0, 360]: {lon_in}")


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


def _detect_lon_type(lon_in):
    if (lon_in < -180) or (lon_in > 360):
        raise ValueError(f"Longitude outside range [-180, 360]: {lon_in}")
    if lon_in < 0:
        lon_type = 180
    elif lon_in > 180:
        lon_type = 360
    else:
        msg = "When providing an ambiguous longitude, you must specify --lon-type 180 or 360"
        raise ArgumentTypeError(msg)
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


def convert_number_to_lon(other):
    """
    Try to convert a given number to Longitude type
    """
    try:
        other_type = _detect_lon_type(other)
    except ArgumentTypeError as err:
        msg = "Trying to compare Longitude class with an ambiguous or invalid value"
        raise TypeError(msg) from err
    other = Longitude(other, other_type)
    return other


@total_ordering  # Avoid having to specify all comparison operators
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

    # __eq__ makes it so that == and != both work.
    def __eq__(self, other):
        if not isinstance(other, Longitude):
            # We could try to coerce non-Longitude `other` to Longitude, but that might result in
            # situations where tests think everything works but code will fail if `other` is
            # ambiguous.
            raise TypeError(
                f"Comparison not supported between instances of 'Longitude' and '{type(other)}'"
            )
        return self._lon == other.get(self._lon_type)

    def __lt__(self, other):
        if not isinstance(other, Longitude):
            # We could try to coerce non-Longitude `other` to Longitude, but that might result in
            # situations where tests think everything works but code will fail if `other` is
            # ambiguous.
            raise TypeError(
                f"Comparison not supported between instances of 'Longitude' and '{type(other)}'"
            )
        if self.lon_type() != other.lon_type():
            # If you're comparing two Longitudes of different types in different hemispheres, then
            # `lon1 > lon2` and `lon2 < lon1` will give different answers! We could make it so that
            # this doesn't fail as long as symmetricality isn't violated, but that might lead to
            # unexpected failures in practice.
            raise TypeError("Comparison not supported between Longitudes of different types")
        return self._lon < other._lon

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
