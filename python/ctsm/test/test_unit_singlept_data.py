#!/usr/bin/env python3
"""
Unit tests for SinglePointCase

You can run this by:
    python -m unittest test_unit_singlept_data.py
"""

import unittest
import argparse
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.single_point_case import SinglePointCase

# pylint: disable=invalid-name


class TestSinglePointCase(unittest.TestCase):
    """
    Basic class for testing SinglePointCase class in single_point_case.py.
    """

    plat = 20.1
    plon = 50.5
    site_name = None
    create_domain = True
    create_surfdata = True
    create_landuse = True
    create_datm = True
    create_user_mods = True
    dom_pft = [8]
    evenly_split_cropland = False
    pct_pft = None
    num_pft = 16
    cth = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    cbh = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    include_nonveg = False
    uni_snow = True
    cap_saturation = True
    out_dir = os.getcwd()
    overwrite = False

    def test_create_tag_noname(self):
        """
        Test create_tag when site_name is NOT given.
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )

        single_point.create_tag()
        self.assertEqual(single_point.tag, "50.5_20.1")

    def test_create_tag_name(self):
        """
        Test create_tag when site_name is given.
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.site_name = "foo"
        single_point.create_tag()
        self.assertEqual(single_point.tag, "foo")

    def test_check_dom_pft_too_big(self):
        """
        Test check_dom_pft
        When one of the given dom_pft(s) are bigger than 78
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [16, 36, 79]
        with self.assertRaisesRegex(argparse.ArgumentTypeError, "values for --dompft should*"):
            single_point.check_dom_pft()

    def test_check_dom_pft_too_small(self):
        """
        Test check_dom_pft
        When one of the given dom_pft(s) are bigger than 1
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [16, 36, -1]
        with self.assertRaisesRegex(argparse.ArgumentTypeError, "values for --dompft should*"):
            single_point.check_dom_pft()

    def test_check_dom_pft_numpft(self):
        """
        Test check_dom_pft
        When dom_pft > 15 but no crop (aka num_pft =<15)
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [15, 53]
        single_point.num_pft = 16
        with self.assertRaisesRegex(argparse.ArgumentTypeError, "Please use --crop*"):
            single_point.check_dom_pft()

    def test_check_dom_pft_mixed_range(self):
        """
        Test check_dom_pft
        Test if all dom_pft(s) are in the same range of either 1-15 or 16-78
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [1, 5, 15]
        single_point.num_pft = 78
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "You are subsetting using mixed land*"
        ):
            single_point.check_dom_pft()

    def test_check_nonveg_nodompft(self):
        """
        Test check_nonveg
        If include_nonveg =False and no dompft it should complain.
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = None
        single_point.include_nonveg = False
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError,
            "To include non-veg land units, you need to specify*",
        ):
            single_point.check_nonveg()

    def test_check_pct_pft_notsamenumbers(self):
        """
        Test check_pct_pft
        Check if pct_pft is the same length as dom_pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [1, 5]
        single_point.pct_pft = [0.5]
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "Please provide the same number of inputs*"
        ):
            single_point.check_pct_pft()

    def test_check_pct_pft_sum_not1(self):
        """
        Test check_pct_pft
        Check if pct_pft adds up to 1 or 100.
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [1, 5]
        single_point.pct_pft = [0.1, 0.5]
        with self.assertRaisesRegex(
            argparse.ArgumentTypeError, "Sum of --pctpft values should be equal to 1*"
        ):
            single_point.check_pct_pft()

    def test_check_pct_pft_fraction_topct(self):
        """
        Test check_pct_pft
        Check if pct_pft is corretly converted to percent.
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [1, 5, 8]
        single_point.pct_pft = [0.5, 0.4, 0.1]
        single_point.check_pct_pft()
        self.assertEqual(single_point.pct_pft, [50, 40, 10])


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
