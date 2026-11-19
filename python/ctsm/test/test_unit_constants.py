#!/usr/bin/env python3

"""
Unit tests for ctsm.constants
"""

import unittest

from ctsm import unit_testing
from ctsm import constants


class TestConstants(unittest.TestCase):
    """Tests of the shared constants"""

    def test_months_per_year(self):
        """MONTHS_PER_YEAR is the number of calendar months in a year"""
        self.assertEqual(constants.MONTHS_PER_YEAR, 12)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
