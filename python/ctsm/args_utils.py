"""
General-purpose utilities for handling command-line
arguments and flags in ctsm python codes.
"""

from ctsm.config_utils import lon_range_0_to_360

# Types

# Types for command-lines error handling:

def plat_type(plat):
    """
    Function to define lat type for the parser
    and
    raise error if latitude is not between -90 and 90.

    Args:
        plat(str): latitude
    Raises:
        Error when plat (latitude) is not between -90 and 90.
    Returns:
        plat (float): latitude in float
    """
    plat_out = float(plat)
    if (plat_out < -90) or (plat_out > 90):
        raise argparse.ArgumentTypeError(
            "ERROR: Latitude should be between -90 and 90."
        )
    return plat_out

def plon_type(plon):
    """
    Function to define lon type for the parser and
    convert negative longitudes and
    raise error if lon is not between -180 and 360.

    Args:
        plon (str): longitude
    Raises:
        Error (ArgumentTypeError): when longitude is <-180 and >360.
    Returns:
        plon(float): converted longitude between 0 and 360
    """
    plon_out = float(plon)
    if -180 <= plon < 0:
        logger.debug("lon is: %f", plon)
        plon = plon % 360
        logger.debug("after modulo lon is: %f", plon)
    if plon < 0 or plon > 360:
        raise argparse.ArgumentTypeError("ERROR: Longitude of single point should be between 0 and "
                                         "360 or -180 and 180.")
    return plon




