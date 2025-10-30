#!/usr/bin/env python3

"""
Unit tests for generate_gdds.py
"""

import unittest
import os

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

    def test_generate_gdds_get_mxmats_ctsm51(self):
        """Test importing from a ctsm51 paramfile (no fail)"""
        paramfile = self._paramfile_51
        gg._get_max_growing_season_lengths(self.no_mxmats, paramfile)

    def test_generate_gdds_get_mxmats_ctsm60(self):
        """Test importing from a ctsm60 paramfile (no fail)"""
        paramfile = self._paramfile_60
        gg._get_max_growing_season_lengths(self.no_mxmats, paramfile)

    def test_generate_gdds_get_mxmats_none(self):
        """Test not importing from a paramfile (should return None)"""
        max_season_length_from_hdates_file = True
        paramfile = None
        result = gg._get_max_growing_season_lengths(max_season_length_from_hdates_file, paramfile)
        self.assertIsNone(result)

    def test_generate_gdds_get_mxmats_values(self):
        """Check values"""
        mxmats = gg._get_max_growing_season_lengths(self.no_mxmats, self.paramfile)

        # Check values
        self.assertTrue(np.isinf(mxmats["needleleaf_evergreen_temperate_tree"]))
        self.assertEqual(mxmats["temperate_corn"], 165)
        self.assertEqual(mxmats["miscanthus"], 210)
        


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
