#!/usr/bin/env python3

"""Unit tests for add_tag_to_filename
"""

import unittest

from ctsm import unit_testing
from ctsm.utils import add_tag_to_filename

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestUtilsAddTag(unittest.TestCase):
    """Tests of utils: add_tag_to_filename"""

    def testSimple(self):
        """Simple test of surface dataset name"""

        fsurf_in = "surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_c221105.nc"
        fsurf_out = add_tag_to_filename(fsurf_in, "tag")
        expect_fsurf = "surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_tag_c221105.nc"
        self.assertEqual(expect_fsurf, fsurf_out, "Expect filenames to be as expected")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
