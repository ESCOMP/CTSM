"""
General-purpose utilities for handling command-line
arguments and flags in ctsm python codes.
Types for command-lines error handling.
"""

import logging
import argparse

logger = logging.getLogger(__name__)


def plat_type(plat):
    """
    Function to define lat type for the parser
    and
    raise error if latitude is not between -90 and 90.

    Args:
        plat(str): latitude
    Raises:
        Error (ArgumentTypeError): when plat (latitude) is not between -90 and 90.
    Returns:
        plat_out (float): latitude in float
    """
    plat_out = float(plat)
    if plat_out < -90 or plat_out > 90:
        raise argparse.ArgumentTypeError("ERROR: Latitude should be between -90 and 90.")
    return plat_out
