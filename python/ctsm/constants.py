"""
Constants shared between CTSM's python modules

Python has no true constants. By convention the names are ALL_CAPS, and the typing.Final
qualifier lets a type checker (e.g., mypy) warn if one is reassigned. Nothing stops
reassignment at runtime, so treat these as read-only.

Constants relating to PFT numbers live in ctsm.pft_utils.
"""

from typing import Final

MONTHS_PER_YEAR: Final[int] = 12
