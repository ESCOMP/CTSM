#!/usr/bin/env python3

"""
Unit tests for generate_gdds_functions.py
"""

import unittest
import os
import tempfile
import shutil
import logging

import numpy as np
import xarray as xr
from cftime import DatetimeNoLeap, DatetimeAllLeap

from ctsm import unit_testing
from ctsm.crop_calendars import generate_gdds_functions as gf

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access

## Too many instant variables as part of the class (too many self.<varible> in the SetUp)
# pylint: disable=too-many-instance-attributes


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
