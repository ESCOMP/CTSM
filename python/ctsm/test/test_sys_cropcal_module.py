#!/usr/bin/env python3

"""
System tests for cropcal_module.py
"""

import unittest
import os

from ctsm import unit_testing
from ctsm.crop_calendars.cropcal_module import import_max_gs_length

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access

## Too many instant variables as part of the class (too many self.<varible> in the SetUp)
# pylint: disable=too-many-instance-attributes


class TestImportMaxGsLength(unittest.TestCase):
    """Tests of import_max_gs_length()"""

    def setUp(self):
        self._paramfile_51 = os.path.join(
            os.path.dirname(__file__), "testinputs", "ctsm51_params.c211112.nc"
        )
        self._paramfile_60 = os.path.join(
            os.path.dirname(__file__), "testinputs", "ctsm60_params_cal115_c250813.nc"
        )

    def test_import_max_gs_length_ctsm51(self):
        """Basic test with ctsm51 paramfile"""
        import_max_gs_length(self._paramfile_51)

    def test_import_max_gs_length_ctsm60(self):
        """Basic test with ctsm60 paramfile"""
        import_max_gs_length(self._paramfile_60)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
