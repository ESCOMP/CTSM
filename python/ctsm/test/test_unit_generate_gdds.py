#!/usr/bin/env python3

"""
Unit tests for generate_gdds.py
"""

import unittest
import os

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
        ]
        gg._parse_args(args)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
