#!/usr/bin/env python3

"""Unit tests for lilac_make_runtime_inputs"""

import unittest

from ctsm import unit_testing
from ctsm.lilac_make_runtime_inputs import determine_bldnml_opts

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestMakeRuntimeInputs(unittest.TestCase):
    """Tests of lilac_make_runtime_inputs"""

    def test_buildnmlOpts_bgc(self):
        """Test determine_buildnml_opts with bgc_mode='bgc'"""
        bldnml_opts = determine_bldnml_opts(bgc_mode="bgc", crop="off", vichydro="off")
        self.assertRegex(bldnml_opts, r"^ *-bgc bgc *$")

    def test_buildnmlOpts_fates(self):
        """Test determine_buildnml_opts with bgc_mode='fates'"""
        bldnml_opts = determine_bldnml_opts(bgc_mode="fates", crop="off", vichydro="off")
        self.assertRegex(bldnml_opts, r"^ *-bgc fates +-no-megan *$")

    def test_buildnmlOpts_bgcCrop(self):
        """Test determine_buildnml_opts with bgc_mode='bgc' and crop on"""
        bldnml_opts = determine_bldnml_opts(bgc_mode="bgc", crop="on", vichydro="off")
        self.assertRegex(bldnml_opts, r"^ *-bgc bgc +-crop *$")

    def test_buildnmlOpts_spCrop_fails(self):
        """Test determine_buildnml_opts with bgc_mode='sp' and crop on: should fail"""
        with self.assertRaisesRegex(
            SystemExit, "setting crop to 'on' is only compatible with bgc_mode"
        ):
            _ = determine_bldnml_opts(bgc_mode="sp", crop="on", vichydro="off")

    def test_buildnmlOpts_spVic(self):
        """Test determine_buildnml_opts with bgc_mode='sp' and vic on"""
        bldnml_opts = determine_bldnml_opts(bgc_mode="sp", crop="off", vichydro="on")
        self.assertRegex(bldnml_opts, r"^ *-bgc sp +-vichydro *$")

    def test_buildnmlOpts_bgcVic(self):
        """Test determine_buildnml_opts with bgc_mode='bgc' and vic on: should fail"""
        with self.assertRaisesRegex(
            SystemExit,
            "setting vichydro to 'on' is only compatible with bgc_mode of 'sp'",
        ):
            _ = determine_bldnml_opts(bgc_mode="bgc", crop="off", vichydro="on")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
