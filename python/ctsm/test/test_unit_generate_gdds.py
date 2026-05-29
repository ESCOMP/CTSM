#!/usr/bin/env python3

"""
Unit tests for generate_gdds.py
"""

import unittest
import os
import argparse
import tempfile
import shutil
import re

import numpy as np
import xarray as xr
from cftime import DatetimeNoLeap

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
        """
        Should error if both --max-season-length-cushion and --max-season-length-from-hdates-file
        are given
        """
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


class TestGetTimeSliceLists(unittest.TestCase):
    """Tests for _get_time_slice_lists()"""

    def test_generate_gdds_get_time_slice_lists(self):
        """Test that _get_time_slice_lists works with two different years"""
        h1_slices, h2_slices = gg._get_time_slice_lists(1986, 1987)

        # Check h1 slices (annual timesteps - single day)
        expected_h1 = [
            slice("1987-01-01", "1987-01-01"),
            slice("1988-01-01", "1988-01-01"),
            slice("1989-01-01", "1989-01-01"),
        ]
        self.assertEqual(h1_slices, expected_h1)

        # Check h2 slices (daily timesteps - full year)
        # For daily, starts at first_season (1986), not first_season + 1
        expected_h2 = [
            slice("1986-01-02", "1987-01-01"),
            slice("1987-01-02", "1988-01-01"),
            slice("1988-01-02", "1989-01-01"),
        ]
        self.assertEqual(h2_slices, expected_h2)

    def test_generate_gdds_get_time_slice_lists_1yr(self):
        """Test that _get_time_slice_lists works with the same year"""
        h1_slices, h2_slices = gg._get_time_slice_lists(1987, 1987)

        # Check h1 slices
        expected_h1 = [
            slice("1988-01-01", "1988-01-01"),
            slice("1989-01-01", "1989-01-01"),
        ]
        self.assertEqual(h1_slices, expected_h1)

        # Check h2 slices
        # For daily, starts at first_season (1987), not first_season + 1
        expected_h2 = [
            slice("1987-01-02", "1988-01-01"),
            slice("1988-01-02", "1989-01-01"),
        ]
        self.assertEqual(h2_slices, expected_h2)

    def test_generate_gdds_get_time_slice_list_valueerror(self):
        """Test that _get_time_slice_list raises ValueError if last < first"""
        with self.assertRaisesRegex(ValueError, "first_season.* > last_season"):
            gg._get_time_slice_lists(1987, 1986)

    def test_generate_gdds_get_time_slice_lists_typeerror_first(self):
        """Test that _get_time_slice_lists raises TypeError if not given integer first season"""
        with self.assertRaisesRegex(
            TypeError, r"_get_time_slice_list\(\) arguments must be integers"
        ):
            gg._get_time_slice_lists(1986.3, 1987)

    def test_generate_gdds_get_time_slice_lists_typeerror_last(self):
        """Test that _get_time_slice_lists raises TypeError if not given integer last season"""
        with self.assertRaisesRegex(
            TypeError, r"_get_time_slice_list\(\) arguments must be integers"
        ):
            gg._get_time_slice_lists(1986, None)

    def test_generate_gdds_get_time_slice_lists_lengths_match(self):
        """Test that h1 and h2 slice lists have the same length"""
        h1_slices, h2_slices = gg._get_time_slice_lists(2000, 2005)
        self.assertEqual(len(h1_slices), len(h2_slices))
        # Should be last_season - first_season + 2
        self.assertEqual(len(h1_slices), 2005 - 2000 + 2)

    def test_generate_gdds_get_time_slice_lists_h1_start_equals_stop(self):
        """Test that h1 slices are single-day (start == stop)"""
        h1_slices, _ = gg._get_time_slice_lists(2000, 2002)
        for s in h1_slices:  # pylint: disable=not-an-iterable
            self.assertEqual(s.start, s.stop)

    def test_generate_gdds_get_time_slice_lists_h2_year_long(self):
        """Test that h2 slices span one year"""
        _, h2_slices = gg._get_time_slice_lists(2000, 2002)
        for s in h2_slices:  # pylint: disable=not-an-iterable
            # Start should be Jan 2 of year Y
            self.assertIn("-01-02", s.start)
            # Stop should be Jan 1 of year Y+1
            self.assertIn("-01-01", s.stop)
            # Extract years and verify they're consecutive
            start_year = int(s.start[:4])
            stop_year = int(s.stop[:4])
            self.assertEqual(stop_year, start_year + 1)


class TestGetHistoryYrRange(unittest.TestCase):
    """Tests for _get_history_yr_range()"""

    def test_get_history_yr_range_annual(self):
        """Test _get_history_yr_range with annual frequency"""
        result = gg._get_history_yr_range(1986, 1987, "annual")
        # For annual: first_season + 1 through last_season + 2
        expected = range(1987, 1990)
        self.assertEqual(result, expected)

    def test_get_history_yr_range_daily(self):
        """Test _get_history_yr_range with daily frequency"""
        result = gg._get_history_yr_range(1986, 1987, "daily")
        # For daily: first_season through last_season + 1
        expected = range(1986, 1989)
        self.assertEqual(result, expected)

    def test_get_history_yr_range_annual_single_year(self):
        """Test _get_history_yr_range with annual frequency and single year"""
        result = gg._get_history_yr_range(2000, 2000, "annual")
        # Should give 2001, 2002
        expected = range(2001, 2003)
        self.assertEqual(result, expected)

    def test_get_history_yr_range_daily_single_year(self):
        """Test _get_history_yr_range with daily frequency and single year"""
        result = gg._get_history_yr_range(2000, 2000, "daily")
        # Should give 2000, 2001
        expected = range(2000, 2002)
        self.assertEqual(result, expected)

    def test_get_history_yr_range_unknown_freq(self):
        """Test _get_history_yr_range with unknown frequency"""
        with self.assertRaises(NotImplementedError):
            gg._get_history_yr_range(2000, 2001, "monthly")

    def test_get_history_yr_range_lengths_match(self):
        """Test that annual and daily ranges have the same length"""
        annual_range = gg._get_history_yr_range(2000, 2005, "annual")
        daily_range = gg._get_history_yr_range(2000, 2005, "daily")
        self.assertEqual(len(annual_range), len(daily_range))
        # Should be last_season - first_season + 2
        self.assertEqual(len(annual_range), 2005 - 2000 + 2)


class TestGetFileLists(unittest.TestCase):
    """Tests of _get_file_lists()"""

    def setUp(self):
        """
        Set up and change to temporary directory
        """
        self.prev_dir = os.getcwd()
        self.temp_dir = tempfile.mkdtemp()
        os.chdir(self.temp_dir)

    def tearDown(self):
        """
        Delete temporary directory
        """
        os.chdir(self.prev_dir)
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _create_test_file(self, filename):
        """Helper to create an empty test file with time coordinate"""
        filepath = os.path.join(self.temp_dir, filename)

        # Extract date from filename using regex (format: *.h#i.YYYY-MM-DD-*.nc)
        match = re.search(r"(\d{4})-(\d{2})-(\d{2})", filename)
        if match:
            year, month, day = match.groups()
            time_val = DatetimeNoLeap(int(year), int(month), int(day), has_year_zero=True)
        else:
            raise ValueError(f"Could not extract date from filename: {filename}")

        # Create a simple dataset with time coordinate
        time = xr.DataArray([time_val], dims=["time"], name="time")
        ds = xr.Dataset({"time": time})
        ds.to_netcdf(filepath)

        return filepath

    def test_get_file_lists_single_year(self):
        """Test _get_file_lists with a single year of data"""
        # Create h1 and h2 files for 2000 and 2001
        # Also need h2 file for 1999 since daily starts at first_season
        h2_file_1999 = self._create_test_file("test.clm2.h2i.1999-01-02-00000.nc")
        h1_file_2000 = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        h2_file_2000 = self._create_test_file("test.clm2.h2i.2000-01-02-00000.nc")
        h1_file_2001 = self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc")

        # Get time slice lists for first_season=1999, last_season=1999
        # This will give us slices for 2000 and 2001 (h1), and 1999, 2000, 2001 (h2)
        time_slice_lists_list = gg._get_time_slice_lists(1999, 1999)

        h1_file_lists, h2_file_lists = gg._get_file_lists(
            self.temp_dir, time_slice_lists_list, logger=None
        )

        # Should have two lists (one for each year: 2000, 2001)
        self.assertEqual(len(h1_file_lists), 2)
        self.assertEqual(len(h2_file_lists), 2)

        # Check contents of file lists
        # pylint: disable=unsubscriptable-object
        self.assertEqual(len(h1_file_lists[0]), 1)
        self.assertEqual(h1_file_lists[0], [h1_file_2000])
        self.assertEqual(len(h2_file_lists[0]), 1)
        self.assertEqual(h2_file_lists[0], [h2_file_1999])
        self.assertEqual(len(h1_file_lists[1]), 1)
        self.assertEqual(h1_file_lists[1], [h1_file_2001])
        self.assertEqual(len(h2_file_lists[1]), 1)
        self.assertEqual(h2_file_lists[1], [h2_file_2000])

    def test_get_file_lists_multiple_years(self):
        """Test _get_file_lists with multiple years of data"""
        # Create h1 and h2 files for 2000-2002
        # Also need h2 file for 1999 since daily starts at first_season
        h1_files = []
        h2_files = [self._create_test_file("test.clm2.h2i.1999-01-02-00000.nc")]
        for year in [2000, 2001, 2002]:
            h1_files.append(self._create_test_file(f"test.clm2.h1i.{year}-01-01-00000.nc"))
            h2_files.append(self._create_test_file(f"test.clm2.h2i.{year}-01-02-00000.nc"))

        # Get time slice lists for first_season=1999, last_season=2000
        # This will give us slices for 2000, 2001, 2002 (h1) and 1999, 2000, 2001 (h2)
        time_slice_lists_list = gg._get_time_slice_lists(1999, 2000)

        h1_file_lists, h2_file_lists = gg._get_file_lists(
            self.temp_dir, time_slice_lists_list, logger=None
        )

        # Should have three lists (one for each year: 2000, 2001, 2002)
        self.assertEqual(len(h1_file_lists), 3)
        self.assertEqual(len(h2_file_lists), 3)

        # Check contents of file lists
        # pylint: disable=unsubscriptable-object
        for i in range(3):
            self.assertEqual(len(h1_file_lists[i]), 1)
            self.assertEqual(len(h2_file_lists[i]), 1)
            self.assertEqual(h1_file_lists[i], [h1_files[i]])
            self.assertEqual(h2_file_lists[i], [h2_files[i]])

    def test_get_file_lists_multiple_files_per_slice(self):
        """Test _get_file_lists when multiple files fall within a time slice"""
        # Create h1 files for 2000 and 2001 (annual)
        h1_file_2000 = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        h1_file_2001 = self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc")

        # Create multiple h2 files for 1999 and 2000 (daily throughout the year)
        h2_files_1999 = []
        for month in ["01", "06", "12"]:
            h2_files_1999.append(self._create_test_file(f"test.clm2.h2i.1999-{month}-15-00000.nc"))

        h2_files_2000 = []
        for month in ["01", "06", "12"]:
            h2_files_2000.append(self._create_test_file(f"test.clm2.h2i.2000-{month}-15-00000.nc"))

        # Get time slice lists for first_season=1999, last_season=1999
        time_slice_lists_list = gg._get_time_slice_lists(1999, 1999)

        h1_file_lists, h2_file_lists = gg._get_file_lists(
            self.temp_dir, time_slice_lists_list, logger=None
        )

        # Should have two lists (for 2000 and 2001)
        self.assertEqual(len(h1_file_lists), 2)
        self.assertEqual(len(h2_file_lists), 2)

        # Check contents of file lists
        # pylint: disable=unsubscriptable-object
        self.assertEqual(len(h1_file_lists[0]), 1)
        self.assertEqual(h1_file_lists[0], [h1_file_2000])
        self.assertEqual(len(h2_file_lists[0]), 3)
        self.assertEqual(h2_file_lists[0], sorted(h2_files_1999))

        # Check second year (2001)
        self.assertEqual(len(h1_file_lists[1]), 1)
        self.assertEqual(h1_file_lists[1], [h1_file_2001])
        self.assertEqual(len(h2_file_lists[1]), 3)
        self.assertEqual(h2_file_lists[1], sorted(h2_files_2000))

    def test_get_file_lists_no_h1_files(self):
        """Test _get_file_lists when h1 files are missing"""

        time_slice_lists_list = gg._get_time_slice_lists(1999, 1999)

        # Should raise FileNotFoundError when h1 files are not found
        with self.assertRaisesRegex(FileNotFoundError, "No files found matching patterns.*h1.*"):
            gg._get_file_lists(self.temp_dir, time_slice_lists_list, logger=None)

    def test_get_file_lists_no_h2_files(self):
        """Test _get_file_lists when h2 files are missing"""
        growing_season = 2000

        # Create h1 files with data from growing_season year and the next one
        self._create_test_file(f"test.clm2.h1i.{growing_season + 1}-01-01-00000.nc")
        self._create_test_file(f"test.clm2.h1i.{growing_season + 2}-01-01-00000.nc")
        # But don't create h2 files

        time_slice_lists_list = gg._get_time_slice_lists(growing_season, growing_season)

        # Should raise FileNotFoundError when h2 files are not found
        with self.assertRaisesRegex(FileNotFoundError, "No files found matching patterns.*h2.*"):
            gg._get_file_lists(self.temp_dir, time_slice_lists_list, logger=None)

    def test_get_file_lists_h1_missing_a_time_slice(self):
        """
        Test _get_file_lists when h1 files exist but not for one of the needed time slices

        We will be simulating a need for processing the 2000 growing season only. This will require
        data from 2000 and also, because seasons can extend into the next calendar year, 2001.

        Because of how CESM timestamps annual output files, it will be looking for h1 timestamps
        with the years 2001-01-01 (data from 2000) and 2002-01-01 (data from 2001).

        In this test, we will only create one file for h1. We expect an error to be thrown about
        that before the check of the available h2 time steps happens.
        """
        growing_season = 2000

        # Create h1 files with data from growing_season year BUT NOT the next one
        self._create_test_file(f"test.clm2.h1i.{growing_season + 1}-01-01-00000.nc")

        # Get required time slice lists for h1 and h2 files
        time_slice_lists_list = gg._get_time_slice_lists(growing_season, growing_season)

        # Make sure the required time slice list is what we expect
        assert time_slice_lists_list[0] == [
            slice(f"{growing_season+1}-01-01", f"{growing_season+1}-01-01", None),
            slice(f"{growing_season+2}-01-01", f"{growing_season+2}-01-01", None),
        ]

        # Should raise FileNotFoundError when h1 files have no timesteps in slice
        with self.assertRaisesRegex(FileNotFoundError, "No h1 timesteps found in"):
            gg._get_file_lists(self.temp_dir, time_slice_lists_list, logger=None)

    def test_get_file_lists_h2_missing_a_time_slice(self):
        """
        Test _get_file_lists when h2 files exist but not for one of the needed time slices

        We will be simulating a need for processing the 2000 growing season only. This will require
        data from 2000 and also, because seasons can extend into the next calendar year, 2001.

        Because of how CESM timestamps daily output files, it will be looking for h2 timestamps
        starting 2000-01-02 (data from 2000-01-01) through 2002-01-01 (data from 2001-12-31).

        In this test, we will create all necessary files for h1 but only one of the two required for
        h2.

        (As test_get_file_lists_h1_missing_a_time_slice but for h2 instead.)
        """
        growing_season = 2000

        # Create h1 files with data from growing_season year and the next one
        self._create_test_file(f"test.clm2.h1i.{growing_season + 1}-01-01-00000.nc")
        self._create_test_file(f"test.clm2.h1i.{growing_season + 2}-01-01-00000.nc")

        # Create h2 file with data from growing_season year BUT NOT the next one
        self._create_test_file(f"test.clm2.h2i.{growing_season+1}-01-02-00000.nc")

        # Get required time slice lists for h1 and h2 files
        time_slice_lists_list = gg._get_time_slice_lists(growing_season, growing_season)

        # Make sure the required time slice lists are what we expect
        assert time_slice_lists_list[0] == [
            slice(f"{growing_season+1}-01-01", f"{growing_season+1}-01-01", None),
            slice(f"{growing_season+2}-01-01", f"{growing_season+2}-01-01", None),
        ]
        assert time_slice_lists_list[1] == [
            slice(f"{growing_season}-01-02", f"{growing_season+1}-01-01", None),
            slice(f"{growing_season+1}-01-02", f"{growing_season+2}-01-01", None),
        ]

        # Should raise FileNotFoundError when h2 files have no timesteps in slice
        with self.assertRaisesRegex(FileNotFoundError, "No h2 timesteps found in"):
            gg._get_file_lists(self.temp_dir, time_slice_lists_list, logger=None)

    def test_get_file_lists_partial_overlap(self):
        """Test _get_file_lists when some time slices have files and others don't"""
        # Create h1 and h2 files for 2000 only
        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("test.clm2.h2i.2000-01-02-00000.nc")

        # Request time slices for 2000, 2001, 2002 (only 2000 has files)
        time_slice_lists_list = gg._get_time_slice_lists(1999, 2000)

        # Should raise FileNotFoundError when second time slice has no files
        with self.assertRaisesRegex(FileNotFoundError, "No h1 timesteps found in"):
            gg._get_file_lists(self.temp_dir, time_slice_lists_list, logger=None)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
