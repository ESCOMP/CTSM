#!/usr/bin/env python3

"""
Unit tests for generate_gdds.py
"""

import unittest
import os
import argparse

import numpy as np

from ctsm import unit_testing
from ctsm.crop_calendars import generate_gdds as gg

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access

## Too many instant variables as part of the class (too many self.<varible> in the SetUp)
# pylint: disable=too-many-instance-attributes


class TestGenerateGddsArgs(unittest.TestCase):
    """Tests the generate_gdds.py argument parsing"""

    def setUp(self):
        self._input_dir = os.path.join("dummy", "path", "to", "inputdir")
        self._sdates_file = os.path.join("dummy", "path", "to", "sdates")
        self._hdates_file = os.path.join("dummy", "path", "to", "hdates")
        self._paramfile = os.path.join("dummy", "path", "to", "paramfile")

    def test_generate_gdds_args_reqd_shortnames(self):
        """Basic test with all required inputs, short arg names"""
        args = [
            "-i",
            self._input_dir,
            "-1",
            "1986",
            "-n",
            "1987",
            "-sd",
            self._sdates_file,
            "-hd",
            self._hdates_file,
            "--paramfile",
            self._paramfile,
        ]
        gg._parse_args(args)

        # Again, with capital -N option
        args = [
            "-i",
            self._input_dir,
            "-1",
            "1986",
            "-N",
            "1987",
            "-sd",
            self._sdates_file,
            "-hd",
            self._hdates_file,
            "--paramfile",
            self._paramfile,
        ]
        gg._parse_args(args)

    def test_generate_gdds_args_reqd_longnames(self):
        """Basic test with all required inputs, long arg names"""
        args = [
            "--input-dir",
            self._input_dir,
            "--first-season",
            "1986",
            "--last-season",
            "1987",
            "--sdates-file",
            self._sdates_file,
            "--hdates-file",
            self._hdates_file,
            "--paramfile",
            self._paramfile,
        ]
        gg._parse_args(args)

    def test_generate_gdds_args_mxmat_from_hdatefile(self):
        """Test with option to get max season length from hdates file"""
        args = [
            "--input-dir",
            self._input_dir,
            "--first-season",
            "1986",
            "--last-season",
            "1987",
            "--sdates-file",
            self._sdates_file,
            "--hdates-file",
            self._hdates_file,
            "--max-season-length-from-hdates-file",
        ]
        gg._parse_args(args)

    def test_generate_gdds_args_error_with_paramfile_and_nomxmat(self):
        """Should error if both --paramfile and --max-season-length-from-hdates-file are given"""
        args = [
            "--input-dir",
            self._input_dir,
            "--first-season",
            "1986",
            "--last-season",
            "1987",
            "--sdates-file",
            self._sdates_file,
            "--hdates-file",
            self._hdates_file,
            "--max-season-length-from-hdates-file",
            "--paramfile",
            self._paramfile,
        ]
        with self.assertRaises(SystemExit):
            gg._parse_args(args)

    def test_generate_gdds_args_error_with_nomxmat_and_cushion(self):
        """Should error if both --max-season-length-cushion and --max-season-length-from-hdates-file
        are given"""
        args = [
            "--input-dir",
            self._input_dir,
            "--first-season",
            "1986",
            "--last-season",
            "1987",
            "--sdates-file",
            self._sdates_file,
            "--hdates-file",
            self._hdates_file,
            "--max-season-length-from-hdates-file",
            "--max-season-length-cushion",
            "14",
        ]
        with self.assertRaises(argparse.ArgumentError):
            gg._parse_args(args)

    def test_generate_gdds_args_ok_with_nomxmat_and_cushion0(self):
        """As test_generate_gdds_args_error_with_nomxmat_and_cushion, but cushion 0 is ok"""
        args = [
            "--input-dir",
            self._input_dir,
            "--first-season",
            "1986",
            "--last-season",
            "1987",
            "--sdates-file",
            self._sdates_file,
            "--hdates-file",
            self._hdates_file,
            "--max-season-length-from-hdates-file",
            "--max-season-length-cushion",
            "0",
        ]
        gg._parse_args(args)


class TestGetMaxGsLengths(unittest.TestCase):
    """Tests get_max_growing_season_lengths()"""

    def setUp(self):
        self._paramfile_51 = os.path.join(
            os.path.dirname(__file__), "testinputs", "ctsm51_params.c211112.nc"
        )
        self._paramfile_60 = os.path.join(
            os.path.dirname(__file__), "testinputs", "ctsm60_params_cal115_c250813.nc"
        )

        # Default arguments
        self.no_mxmats = False
        self.paramfile = self._paramfile_60
        self.cushion = 0

    def test_generate_gdds_get_mxmats_ctsm51(self):
        """Test importing from a ctsm51 paramfile (no fail)"""
        paramfile = self._paramfile_51
        gg._get_max_growing_season_lengths(self.no_mxmats, paramfile, self.cushion)

    def test_generate_gdds_get_mxmats_ctsm60(self):
        """Test importing from a ctsm60 paramfile (no fail)"""
        paramfile = self._paramfile_60
        gg._get_max_growing_season_lengths(self.no_mxmats, paramfile, self.cushion)

    def test_generate_gdds_get_mxmats_none(self):
        """Test not importing from a paramfile (should return None)"""
        max_season_length_from_hdates_file = True
        paramfile = None
        result = gg._get_max_growing_season_lengths(
            max_season_length_from_hdates_file, paramfile, self.cushion
        )
        self.assertIsNone(result)

    def test_generate_gdds_get_mxmats_values(self):
        """Check values with no cushion"""
        mxmats = gg._get_max_growing_season_lengths(self.no_mxmats, self.paramfile, self.cushion)

        # Check values
        self.assertTrue(np.isinf(mxmats["needleleaf_evergreen_temperate_tree"]))
        self.assertEqual(mxmats["temperate_corn"], 165)
        self.assertEqual(mxmats["miscanthus"], 210)

    def test_generate_gdds_get_mxmats_cushion14(self):
        """As test_generate_gdds_get_mxmats_values, but cushion 14"""
        cushion = 14
        mxmats = gg._get_max_growing_season_lengths(self.no_mxmats, self.paramfile, cushion)

        # Check values
        self.assertTrue(np.isinf(mxmats["needleleaf_evergreen_temperate_tree"]))
        self.assertEqual(mxmats["temperate_corn"], 165 - cushion)
        self.assertEqual(mxmats["miscanthus"], 210 - cushion)

    def test_generate_gdds_get_mxmats_cushionneg14(self):
        """As test_generate_gdds_get_mxmats_values, but cushion -14"""
        cushion = -14
        mxmats = gg._get_max_growing_season_lengths(self.no_mxmats, self.paramfile, cushion)

        # Check values
        self.assertTrue(np.isinf(mxmats["needleleaf_evergreen_temperate_tree"]))
        self.assertEqual(mxmats["temperate_corn"], 165 - cushion)
        self.assertEqual(mxmats["miscanthus"], 210 - cushion)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
