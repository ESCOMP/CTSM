#!/usr/bin/env python3

"""
Unit tests for generate_gdds.py and generate_gdds_functions.py
"""

import unittest
import os
import argparse
import tempfile
import shutil
import logging
import re

import numpy as np
import xarray as xr
from cftime import DatetimeNoLeap, DatetimeAllLeap

from ctsm import unit_testing
from ctsm.crop_calendars import generate_gdds as gg
from ctsm.crop_calendars import generate_gdds_functions as gf

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

    def test_generate_gdds_get_time_slice_lists_h1_single_day(self):
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


class TestCheckGridMatch(unittest.TestCase):
    """Tests check_grid_match()"""

    def test_check_grid_match_true_npnp(self):
        """Test check_grid_match() with two matching numpy arrays"""
        np0 = np.array([0, 1, 2, np.pi])
        match, max_abs_diff = gf.check_grid_match(np0, np0)
        self.assertTrue(match)
        self.assertEqual(max_abs_diff, 0.0)

    def test_check_grid_match_true_dada(self):
        """Test check_grid_match() with two matching DataArrays"""
        np0 = np.array([0, 1, 2, np.pi])
        da0 = xr.DataArray(data=np0)
        match, max_abs_diff = gf.check_grid_match(da0, da0)
        self.assertTrue(match)
        self.assertEqual(max_abs_diff, 0.0)

    def test_check_grid_match_false_npnp(self):
        """Test check_grid_match() with two non-matching numpy arrays"""
        np0 = np.array([0, 1, 2, np.pi])
        np1 = np0.copy()
        diff = 2 * gf.GRID_TOL_DEG
        np1[0] = np0[0] + diff
        match, max_abs_diff = gf.check_grid_match(np0, np1)
        self.assertFalse(match)
        self.assertEqual(max_abs_diff, diff)

    def test_check_grid_match_false_dada(self):
        """Test check_grid_match() with two non-matching DataArrays"""
        np0 = np.array([0, 1, 2, np.pi])
        np1 = np0.copy()
        diff = 2 * gf.GRID_TOL_DEG
        np1[0] = np0[0] + diff
        da0 = xr.DataArray(data=np0)
        da1 = xr.DataArray(data=np1)
        match, max_abs_diff = gf.check_grid_match(da0, da1)
        self.assertFalse(match)
        self.assertEqual(max_abs_diff, diff)

    def test_check_grid_match_falseneg_npnp(self):
        """As test_check_grid_match_false_npnp, but with diff in negative direction"""
        np0 = np.array([0, 1, 2, np.pi])
        np1 = np0.copy()
        diff = -2 * gf.GRID_TOL_DEG
        np1[0] = np0[0] + diff
        match, max_abs_diff = gf.check_grid_match(np0, np1)
        self.assertFalse(match)
        self.assertEqual(max_abs_diff, abs(diff))

    def test_check_grid_match_matchnans_true_npnp(self):
        """Test check_grid_match() with two numpy arrays that have nans and match"""
        np0 = np.array([np.nan, 1, 2, np.pi])
        with self.assertWarnsRegex(RuntimeWarning, r"NaN\(s\) in grid"):
            match, max_abs_diff = gf.check_grid_match(np0, np0)
        self.assertTrue(match)
        self.assertEqual(max_abs_diff, 0.0)

    def test_check_grid_match_matchnans_true_dada(self):
        """Test check_grid_match() with two DataArrays that have nans and match"""
        np0 = np.array([np.nan, 1, 2, np.pi])
        da0 = xr.DataArray(data=np0)
        with self.assertWarnsRegex(RuntimeWarning, r"NaN\(s\) in grid"):
            match, max_abs_diff = gf.check_grid_match(da0, da0)
        self.assertTrue(match)
        self.assertEqual(max_abs_diff, 0.0)

    def test_check_grid_match_matchnans_false_npnp(self):
        """Test check_grid_match() with two numpy arrays with nans that DON'T match"""
        np0 = np.array([np.nan, 1, 2, np.pi])
        np1 = np.array([np.nan, 1, np.nan, np.pi])
        with self.assertWarnsRegex(RuntimeWarning, r"NaN\(s\) in grid don't match"):
            match, max_abs_diff = gf.check_grid_match(np0, np1)
        self.assertFalse(match)
        self.assertIsNone(max_abs_diff)

    def test_check_grid_match_matchnans_falseshape_npnp(self):
        """Test check_grid_match() with two numpy arrays that have different shapes"""
        np0 = np.array([0, 1, 2, np.pi])
        np1 = np.array([0, 1, 2, np.pi, 4])
        match, max_abs_diff = gf.check_grid_match(np0, np1)
        self.assertFalse(match)
        self.assertIsNone(max_abs_diff)

    def test_check_grid_match_matchnans_falseshape_dada(self):
        """Test check_grid_match() with two DataArrays that have different shapes"""
        np0 = np.array([0, 1, 2, np.pi])
        np1 = np.array([0, 1, 2, np.pi, 4])
        da0 = xr.DataArray(data=np0)
        da1 = xr.DataArray(data=np1)
        match, max_abs_diff = gf.check_grid_match(da0, da1)
        self.assertFalse(match)
        self.assertIsNone(max_abs_diff)


class TestFindInstHistFiles(unittest.TestCase):
    """Tests of find_inst_hist_files()"""

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
        """Helper to create an empty test file"""
        filepath = os.path.join(self.temp_dir, filename)
        with open(filepath, "a", encoding="utf-8"):
            pass
        return filepath

    def test_find_inst_hist_files(self):
        """Test finding only h2 files when h1i files present too"""
        # Create test files
        file1 = self._create_test_file("test.clm2.h2i.2000-01-02-00000.nc")
        file2 = self._create_test_file("test.clm2.h2i.2001-01-02-00000.nc")
        # Create h1 file that should not be found
        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")

        result = gf.find_inst_hist_files(self.temp_dir, h=2)

        # Should find only h2i files
        self.assertEqual(len(result), 2)
        self.assertIn(file1, result)
        self.assertIn(file2, result)

    def test_find_inst_hist_files_prefer_nc_over_base(self):
        """Test that .nc files are preferred over .nc.base files"""
        # Create both .nc and .nc.base files
        file_nc = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        file_nc_base = self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc.base")

        result = gf.find_inst_hist_files(self.temp_dir, h=1)

        # Should find .nc files first (pattern order preference)
        self.assertIn(file_nc, result)
        self.assertNotIn(file_nc_base, result)
        # Should have only 1 file (the .nc file, not the .nc.base)
        self.assertEqual(len(result), 1)

    def test_find_inst_hist_files_base_only(self):
        """Test finding files when only .nc.base files exist"""
        # Create only .nc.base files
        file1 = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc.base")
        file2 = self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc.base")

        result = gf.find_inst_hist_files(self.temp_dir, h=1)

        # Should find .nc.base files when no .nc files exist
        self.assertEqual(len(result), 2)
        self.assertIn(file1, result)
        self.assertIn(file2, result)

    def test_find_inst_hist_files_multiple_months_same_year(self):
        """Test finding multiple files from the same year"""
        # Create multiple files from 2000
        file1 = self._create_test_file("test.clm2.h2i.2000-01-01-00000.nc")
        file2 = self._create_test_file("test.clm2.h2i.2000-01-15-00000.nc")
        file3 = self._create_test_file("test.clm2.h2i.2000-01-31-00000.nc")
        # Create file from different year
        file4 = self._create_test_file("test.clm2.h2i.2001-01-01-00000.nc")

        result = gf.find_inst_hist_files(self.temp_dir, h=2)

        # Should find all files
        self.assertEqual(len(result), 4)
        self.assertIn(file1, result)
        self.assertIn(file2, result)
        self.assertIn(file3, result)
        self.assertIn(file4, result)

    def test_find_inst_hist_files_no_files_found(self):
        """Test error when no matching files are found"""
        # Create a non-matching file
        self._create_test_file("test.clm2.h0.2000-01-01-00000.nc")

        # Should raise a FileNotFoundError error
        with self.assertRaisesRegex(FileNotFoundError, "No files found matching patterns"):
            gf.find_inst_hist_files(self.temp_dir, h=1)

    def test_find_inst_hist_files_different_case_names(self):
        """Test that RuntimeError is raised when files from different case names are found"""
        # Create files with different case names
        self._create_test_file("case1.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("case2.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("longcasename.clm2.h1i.2000-01-01-00000.nc")

        # Should raise RuntimeError due to multiple case names
        with self.assertRaisesRegex(RuntimeError, "Found files from multiple case names"):
            gf.find_inst_hist_files(self.temp_dir, h=1)

    def test_find_inst_hist_files_different_case_names_with_logger(self):
        """
        Test that RuntimeError is raised when files from different case names are found, with logger
        """
        # Create a logger
        logger = logging.getLogger("test_logger_case_names")
        logger.setLevel(logging.DEBUG)

        # Create files with different case names
        self._create_test_file("case1.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("case2.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("longcasename.clm2.h1i.2000-01-01-00000.nc")

        # Should raise RuntimeError due to multiple case names, even with logger
        with self.assertRaisesRegex(RuntimeError, "Found files from multiple case names"):
            gf.find_inst_hist_files(self.temp_dir, h=1, logger=logger)

    def test_find_inst_hist_files_no_files_found_with_logger(self):
        """Test error when no matching files are found, with logger"""
        # Create a logger
        logger = logging.getLogger("test_logger_no_files")
        logger.setLevel(logging.DEBUG)

        # Create a non-matching file
        self._create_test_file("test.clm2.h0.2000-01-01-00000.nc")

        # Should raise a FileNotFoundError even with logger
        with self.assertRaisesRegex(FileNotFoundError, "No files found matching patterns"):
            gf.find_inst_hist_files(self.temp_dir, h=1, logger=logger)

    def test_find_inst_hist_files_h_str_with_logger(self):
        """Test that TypeError is raised when h is a string, with logger"""
        # Create a logger
        logger = logging.getLogger("test_logger_h_str")
        logger.setLevel(logging.DEBUG)

        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")

        with self.assertRaisesRegex(TypeError, "must be an integer, not"):
            gf.find_inst_hist_files(self.temp_dir, h="1", logger=logger)

    def test_find_inst_hist_files_h_float_with_logger(self):
        """Test that TypeError is raised when h is a float, with logger"""
        # Create a logger
        logger = logging.getLogger("test_logger_h_float")
        logger.setLevel(logging.DEBUG)

        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")

        with self.assertRaisesRegex(TypeError, "must be an integer, not"):
            gf.find_inst_hist_files(self.temp_dir, h=1.0, logger=logger)


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


class TestGenInstDailyYear(unittest.TestCase):
    """Tests of gen_inst_daily_year()"""

    def test_gen_inst_daily_year_noleap(self):
        """Test gen_inst_daily_year with DatetimeNoLeap calendar"""
        result = gf.gen_inst_daily_year(2000, DatetimeNoLeap)

        # Should have 365 timesteps (Jan 2 through Jan 1 of next year, inclusive)
        self.assertEqual(len(result), 365)

        # Check first and last dates
        self.assertEqual(result.values[0], DatetimeNoLeap(2000, 1, 2, has_year_zero=True))
        self.assertEqual(result.values[-1], DatetimeNoLeap(2001, 1, 1, has_year_zero=True))

        # Verify no Feb 29
        dates_str = [str(d) for d in result.values]
        self.assertNotIn("2000-02-29", dates_str)

        # Check that all timestamps are at exactly midnight
        for t in result.values:
            assert t.hour == 0
            assert t.minute == 0
            assert t.second == 0
            assert t.microsecond == 0

    def test_gen_inst_daily_year_leap(self):
        """Test gen_inst_daily_year with a leap year"""
        cal_type = DatetimeAllLeap
        year = 2004
        result = gf.gen_inst_daily_year(year, cal_type)

        # Should have 366 timesteps (Jan 2 through Jan 1 of next year, inclusive, with leap day)
        self.assertEqual(len(result), 366)

        # Check first and last dates
        self.assertEqual(result.values[0], cal_type(year, 1, 2, has_year_zero=True))
        self.assertEqual(result.values[-1], cal_type(year + 1, 1, 1, has_year_zero=True))

        # Verify Feb 29 is there
        dates_str = [str(d) for d in result.values]
        self.assertIn(f"{year}-02-29 00:00:00", dates_str)

        # Check that all timestamps are at exactly midnight
        for t in result.values:
            assert t.hour == 0
            assert t.minute == 0
            assert t.second == 0
            assert t.microsecond == 0

    def test_gen_inst_daily_year_consecutive_days(self):
        """Test that gen_inst_daily_year produces consecutive days"""
        result = gf.gen_inst_daily_year(2000, DatetimeNoLeap)

        # Check that each day is exactly one day after the previous
        for i in range(1, len(result)):
            diff = result.values[i] - result.values[i - 1]
            self.assertEqual(diff.days, 1)

    def test_gen_inst_daily_year_different_year(self):
        """Test gen_inst_daily_year with a different year"""
        result = gf.gen_inst_daily_year(1987, DatetimeNoLeap)

        # Should still have 365 timesteps
        self.assertEqual(len(result), 365)

        # Check first and last dates
        self.assertEqual(result.values[0], DatetimeNoLeap(1987, 1, 2, has_year_zero=True))
        self.assertEqual(result.values[-1], DatetimeNoLeap(1988, 1, 1, has_year_zero=True))


class TestCheckTimeDa(unittest.TestCase):
    """Tests of _check_time_da()"""

    def test_check_time_da_annual_correct(self):
        """Test _check_time_da with correct annual data"""
        time_val = DatetimeNoLeap(2000, 1, 1, has_year_zero=True)
        time_da = xr.DataArray([time_val], dims=["time"], coords={"time": [time_val]})

        # Should not raise an error
        gf._check_time_da("annual", 2000, time_da, logger=None)

    def test_check_time_da_annual_incorrect(self):
        """Test _check_time_da with incorrect annual data"""
        # Wrong date (Jan 2 instead of Jan 1)
        time_val = DatetimeNoLeap(2000, 1, 2, has_year_zero=True)
        time_da = xr.DataArray([time_val], dims=["time"], coords={"time": [time_val]})

        # Should raise AssertionError (via error())
        with self.assertRaises(AssertionError):
            gf._check_time_da("annual", 2000, time_da, logger=None)

    def test_check_time_da_daily_correct(self):
        """Test _check_time_da with correct daily data"""
        time_da = gf.gen_inst_daily_year(2000, DatetimeNoLeap)

        # Should not raise an error
        gf._check_time_da("daily", 2000, time_da, logger=None)

    def test_check_time_da_daily_correct_leap(self):
        """As test_check_time_da_daily_correct, but for a leap year"""
        time_da = gf.gen_inst_daily_year(2000, DatetimeAllLeap)

        # Should not raise an error
        gf._check_time_da("daily", 2000, time_da, logger=None)

    def test_check_time_da_daily_missing_day(self):
        """Test _check_time_da with daily data missing a day"""
        time_da = gf.gen_inst_daily_year(2000, DatetimeNoLeap)
        # Remove one day from the middle (day 100)
        time_da_missing = xr.concat(
            [time_da.isel(time=slice(0, 100)), time_da.isel(time=slice(101, None))],
            dim="time",
        )

        # Should raise AssertionError (via error())
        with self.assertRaises(AssertionError):
            gf._check_time_da("daily", 2000, time_da_missing, logger=None)

    def test_check_time_da_daily_extra_day(self):
        """Test _check_time_da with daily data having an extra day"""
        time_da = gf.gen_inst_daily_year(2000, DatetimeNoLeap)
        # Add an extra day
        extra_day = DatetimeNoLeap(2001, 1, 3, has_year_zero=True)
        time_da = xr.concat(
            [time_da, xr.DataArray([extra_day], dims=["time"], coords={"time": [extra_day]})],
            dim="time",
        )

        # Should raise AssertionError (via error())
        with self.assertRaises(AssertionError):
            gf._check_time_da("daily", 2000, time_da, logger=None)

    def test_check_time_da_unknown_freq(self):
        """Test _check_time_da with unknown frequency"""
        time_val = DatetimeNoLeap(2000, 1, 1, has_year_zero=True)
        time_da = xr.DataArray([time_val], dims=["time"], coords={"time": [time_val]})

        # Should raise NotImplementedError for unknown frequency
        with self.assertRaises(NotImplementedError):
            gf._check_time_da("monthly", 2000, time_da, logger=None)


class TestCheckFileLists(unittest.TestCase):
    """Tests of check_file_lists()"""

    def setUp(self):
        """Set up and change to temporary directory"""
        self.prev_dir = os.getcwd()
        self.temp_dir = tempfile.mkdtemp()
        os.chdir(self.temp_dir)

    def tearDown(self):
        """Delete temporary directory"""
        os.chdir(self.prev_dir)
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _create_test_file_with_time(self, filename, time_values):
        """Helper to create a test file with specific time values"""
        filepath = os.path.join(self.temp_dir, filename)
        time = xr.DataArray(time_values, dims=["time"], name="time")
        ds = xr.Dataset({"time": time})
        ds.to_netcdf(filepath)
        return filepath

    def test_check_file_lists_annual_correct(self):
        """Test check_file_lists with correct annual data"""
        # Create files with correct annual timesteps
        file1 = self._create_test_file_with_time(
            "test.h1i.2000.nc", [DatetimeNoLeap(2000, 1, 1, has_year_zero=True)]
        )
        file2 = self._create_test_file_with_time(
            "test.h1i.2001.nc", [DatetimeNoLeap(2001, 1, 1, has_year_zero=True)]
        )

        history_yr_range = range(2000, 2002)
        h_file_lists = [[file1], [file2]]
        time_slice_list = [
            slice("2000-01-01", "2000-01-01"),
            slice("2001-01-01", "2001-01-01"),
        ]

        # Should not raise an error
        gf.check_file_lists(history_yr_range, h_file_lists, time_slice_list, "annual", logger=None)

    def test_check_file_lists_annual_correct_extrafile(self):
        """Test check_file_lists with correct annual data but an extra file"""
        # Create files with correct annual timesteps
        file1 = self._create_test_file_with_time(
            "test.h1i.2000.nc", [DatetimeNoLeap(2000, 1, 1, has_year_zero=True)]
        )
        file2 = self._create_test_file_with_time(
            "test.h1i.2001.nc", [DatetimeNoLeap(2001, 1, 1, has_year_zero=True)]
        )
        file3 = self._create_test_file_with_time(
            "test.h1i.2002.nc", [DatetimeNoLeap(2002, 1, 1, has_year_zero=True)]
        )

        history_yr_range = range(2000, 2002)
        h_file_lists = [[file1], [file2], [file3]]
        time_slice_list = [
            slice("2000-01-01", "2000-01-01"),
            slice("2001-01-01", "2001-01-01"),
        ]

        # Should not raise an error
        gf.check_file_lists(history_yr_range, h_file_lists, time_slice_list, "annual", logger=None)

    def test_check_file_lists_annual_incorrect(self):
        """Test check_file_lists with incorrect annual data (wrong day)"""
        # Create file with wrong date (Jan 2 instead of Jan 1)
        file1 = self._create_test_file_with_time(
            "test.h1i.2000.nc",
            [DatetimeNoLeap(2000, 1, 2, has_year_zero=True)],  # Wrong day
        )

        history_yr_range = range(2000, 2001)
        h_file_lists = [[file1]]
        # Use a slice that will include Jan 2
        time_slice_list = [slice("2000-01-01", "2000-01-02")]

        # Should raise AssertionError (has Jan 2 instead of Jan 1)
        with self.assertRaises(AssertionError):
            gf.check_file_lists(
                history_yr_range, h_file_lists, time_slice_list, "annual", logger=None
            )

    def test_check_file_lists_daily_correct(self):
        """Test check_file_lists with correct daily data"""
        # Create file with correct daily timesteps
        time_da = gf.gen_inst_daily_year(2000, DatetimeNoLeap)
        file1 = self._create_test_file_with_time("test.h2i.2000.nc", time_da.values)

        history_yr_range = range(2000, 2001)
        h_file_lists = [[file1]]
        time_slice_list = [slice("2000-01-02", "2001-01-01")]

        # Should not raise an error
        gf.check_file_lists(history_yr_range, h_file_lists, time_slice_list, "daily", logger=None)

    def test_check_file_lists_daily_correct_leap(self):
        """As test_check_file_lists_daily_correct but with a leap day"""
        # Create file with correct daily timesteps
        time_da = gf.gen_inst_daily_year(2000, DatetimeAllLeap)
        file1 = self._create_test_file_with_time("test.h2i.2000.nc", time_da.values)

        history_yr_range = range(2000, 2001)
        h_file_lists = [[file1]]
        time_slice_list = [slice("2000-01-02", "2001-01-01")]

        # Should not raise an error
        gf.check_file_lists(history_yr_range, h_file_lists, time_slice_list, "daily", logger=None)

    def test_check_file_lists_daily_missing_day(self):
        """Test check_file_lists with daily data missing a day"""
        # Create file with incomplete daily timesteps
        time_da = gf.gen_inst_daily_year(2000, DatetimeNoLeap)
        # Remove one day from the middle
        time_da_incomplete = xr.concat(
            [time_da.isel(time=slice(0, 100)), time_da.isel(time=slice(101, None))],
            dim="time",
        )
        file1 = self._create_test_file_with_time("test.h2i.2000.nc", time_da_incomplete.values)

        history_yr_range = range(2000, 2001)
        h_file_lists = [[file1]]
        time_slice_list = [slice("2000-01-02", "2001-01-01")]

        # Should raise AssertionError
        with self.assertRaises(AssertionError):
            gf.check_file_lists(
                history_yr_range, h_file_lists, time_slice_list, "daily", logger=None
            )

    def test_check_file_lists_multiple_files_per_year(self):
        """Test check_file_lists with multiple files per year"""
        # Create two files that together have all daily timesteps for 2000
        time_da_full = gf.gen_inst_daily_year(2000, DatetimeNoLeap)
        time_da_first_half = time_da_full.isel(time=slice(0, 182))
        time_da_second_half = time_da_full.isel(time=slice(182, None))

        file1 = self._create_test_file_with_time("test.h2i.2000a.nc", time_da_first_half.values)
        file2 = self._create_test_file_with_time("test.h2i.2000b.nc", time_da_second_half.values)

        history_yr_range = range(2000, 2001)
        h_file_lists = [[file1, file2]]
        time_slice_list = [slice("2000-01-02", "2001-01-01")]

        # Should not raise an error
        gf.check_file_lists(history_yr_range, h_file_lists, time_slice_list, "daily", logger=None)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
