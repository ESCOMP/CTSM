#!/usr/bin/env python3

"""Unit tests for add_tag_to_filename"""

import unittest

from unittest.mock import patch
from datetime import date
from ctsm import unit_testing

from ctsm import utils

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestUtilsAddTag(unittest.TestCase):
    """Tests of utils: add_tag_to_filename"""

    @staticmethod
    def _fake_today():
        """Set the fake date to Halloween"""
        return date(year=2022, month=10, day=31)

    def testSimple(self):
        """Simple test of surface dataset name"""

        fsurf_in = "surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_c221105.nc"
        with patch("ctsm.utils.date") as mock_date:
            mock_date.today.side_effect = self._fake_today

            fsurf_out = utils.add_tag_to_filename(fsurf_in, "tag")
            fsurf_out2 = utils.add_tag_to_filename(fsurf_in, "tag", replace_res=True)

        expect_fsurf = "surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_tag_c221031.nc"
        self.assertEqual(expect_fsurf, fsurf_out, "Expect filenames to be as expected")
        expect_fsurf2 = "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc"
        self.assertEqual(expect_fsurf2, fsurf_out2, "Expect filenames to be as expected")

    def testSimpleLanduse(self):
        """Simple test of landuse dataset name"""

        landuse_in = "landuse.timeseries_0.9x1.25_hist_78pfts_CMIP6_simyr1850-2015_c190214.nc"
        with patch("ctsm.utils.date") as mock_date:
            mock_date.today.side_effect = self._fake_today

            landuse_out = utils.add_tag_to_filename(landuse_in, "tag")
            landuse_out2 = utils.add_tag_to_filename(landuse_in, "tag", replace_res=True)

        expect_landuse = (
            "landuse.timeseries_0.9x1.25_hist_78pfts_CMIP6_simyr1850-2015_tag_c221031.nc"
        )
        self.assertEqual(expect_landuse, landuse_out, "Expect filenames to be as expected")
        expect_landuse2 = "landuse.timeseries_tag_hist_78pfts_CMIP6_simyr1850-2015_c221031.nc"
        self.assertEqual(expect_landuse2, landuse_out2, "Expect filenames to be as expected")

    def testSimpleDatmDomain(self):
        """Simple test of datm domain dataset name"""

        file_in = "domain.lnd.360x720_gswp3.0v1.c170606.nc"
        with patch("ctsm.utils.date") as mock_date:
            mock_date.today.side_effect = self._fake_today

            file_out = utils.add_tag_to_filename(file_in, "tag")

        expect_filename = "domain.lnd.360x720_gswp3.0v1_tag_c221031.nc"
        self.assertEqual(expect_filename, file_out, "Expect filenames to be as expected")

    def testSimpleDomain(self):
        """Simple test of domain dataset name"""

        file_in = "domain.lnd.fv0.9x1.25_gx1v7.151020.nc"
        with patch("ctsm.utils.date") as mock_date:
            mock_date.today.side_effect = self._fake_today

            file_out = utils.add_tag_to_filename(file_in, "tag")

        expect_filename = "domain.lnd.fv0.9x1.25_gx1v7_tag_c221031.nc"
        self.assertEqual(expect_filename, file_out, "Expect filenames to be as expected")

    def testSurfReplaceListDomain(self):
        """Simple test of list of surface dataset name with replace_res option"""

        files_in = [
            "surfdata_48x96_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc",
            "surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc",
            "surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc",
            "surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304.nc",
            "surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304.nc",
            "surfdata_4x5_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc",
            "surfdata_10x15_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc",
            "surfdata_10x15_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc",
            "surfdata_0.125nldas2_hist_16pfts_Irrig_CMIP6_simyr2005_c190412.nc",
            "surfdata_64x128_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc",
            "surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr2000_c190214.nc",
            "surfdata_1.9x2.5_hist_78pfts_CMIP6_simyr2000_c190304.nc",
            "surfdata_0.125x0.125_hist_78pfts_CMIP6_simyr2005_c190624.nc",
            "surfdata_10x15_hist_78pfts_CMIP6_simyr2000_c190214.nc",
            "surfdata_4x5_hist_78pfts_CMIP6_simyr2000_c190214.nc",
            "surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr1850_c190304.nc",
            "surfdata_10x15_hist_16pfts_Irrig_CMIP6_simyr1850_c190214.nc",
            "surfdata_4x5_hist_16pfts_Irrig_CMIP6_simyr1850_c190214.nc",
            "surfdata_48x96_hist_78pfts_CMIP6_simyr1850_c190214.nc",
            "surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr1850_c190214.nc",
            "surfdata_1.9x2.5_hist_78pfts_CMIP6_simyr1850_c190304.nc",
            "surfdata_10x15_hist_78pfts_CMIP6_simyr1850_c190214.nc",
            "surfdata_4x5_hist_78pfts_CMIP6_simyr1850_c190214.nc",
            "surfdata_ne0np4.ARCTICGRIS.ne30x8_hist_78pfts_CMIP6_simyr2000_c200426.nc",
            "surfdata_C96_hist_78pfts_CMIP6_simyr1850_c200317.nc",
            "surfdata_C96_hist_78pfts_CMIP6_simyr1850_c20221108.nc",
            "surfdata_0.9x1.25_hist_16pfts_nourb_CMIP6_simyrPtVg_c181114.nc",
        ]
        expect_filenames = [
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2005_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr2005_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_16pfts_Irrig_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr2000_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_78pfts_CMIP6_simyr1850_c221031.nc",
            "surfdata_tag_hist_16pfts_nourb_CMIP6_simyrPtVg_c221031.nc",
        ]
        self.assertEqual(
            len(files_in), len(expect_filenames), "length of arrays does not match as expected"
        )
        for i, file_in in enumerate(files_in):

            with patch("ctsm.utils.date") as mock_date:
                mock_date.today.side_effect = self._fake_today

                file_out = utils.add_tag_to_filename(file_in, "tag", replace_res=True)

            self.assertEqual(expect_filenames[i], file_out, "Expect filenames to be as expected")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
