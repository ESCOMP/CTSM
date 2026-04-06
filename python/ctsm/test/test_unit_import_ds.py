#!/usr/bin/env python3

"""
Unit tests for import_ds.py
"""

import unittest
import os
import tempfile
import shutil

import xarray as xr
from cftime import DatetimeNoLeap

from ctsm import unit_testing
from ctsm.crop_calendars import import_ds

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


def _make_timestep(str_in):
    """
    Because of float imprecision, microseconds should be specified like:
        1,8 instead of 1.000008
    and:
        1,800000 instead of 1.8
    """
    h = minute = s = us = 0

    str_in_split = str_in.split(" ")
    y, month, d = str_in_split[0].split("-")
    if len(str_in_split) > 1:
        h, minute, s = str_in_split[1].split(":")
        if "," in s:
            s, us = s.split(",")
    inputs = [int(x) for x in [y, month, d, h, minute, s, us]]
    return DatetimeNoLeap(*inputs, has_year_zero=True)


class TestMakeTimestep(unittest.TestCase):
    """Test this test module's _make_timestep() function"""

    def test_make_timestep_ymd(self):
        """Test with YYYY-MM-DD"""
        self.assertEqual(
            _make_timestep("1987-07-24"),
            DatetimeNoLeap(1987, 7, 24, 0, 0, 0, 0, has_year_zero=True),
        )

    def test_make_timestep_hms(self):
        """Test with YYYY-MM-DD hh:mm:ss"""
        self.assertEqual(
            _make_timestep("1987-07-24 09:25:07"),
            DatetimeNoLeap(1987, 7, 24, 9, 25, 7, 0, has_year_zero=True),
        )

    def test_make_timestep_microsec_leadzeros(self):
        """Test with microseconds with leading zeros"""
        self.assertEqual(
            _make_timestep("1987-07-24 09:25:07,8"),
            DatetimeNoLeap(1987, 7, 24, 9, 25, 7, 8, has_year_zero=True),
        )

    def test_make_timestep_microsec_noleadzeros(self):
        """Test with microseconds without leading zeros"""
        self.assertEqual(
            _make_timestep("1987-07-24 09:25:07,800000"),
            DatetimeNoLeap(1987, 7, 24, 9, 25, 7, 800000, has_year_zero=True),
        )


class TestGetFilesInTimeSlice(unittest.TestCase):
    """Tests of get_files_in_time_slice()"""

    def setUp(self):
        """
        Set up and change to temporary directory
        """
        self.prev_dir = os.getcwd()
        self.temp_dir = tempfile.mkdtemp()
        os.chdir(self.temp_dir)

    def tearDown(self):
        """
        Delete temporary directory and any files within
        """
        os.chdir(self.prev_dir)
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _create_annual_test_files(self, years):
        """
        Helper method to create test files with one timestep per file (annual)

        Args:
            years: List of years to create files for

        Returns:
            List of filenames created
        """
        filelist = []
        for year in years:
            filename = os.path.join(self.temp_dir, f"test_{year}.nc")
            filelist.append(filename)

            # Create a simple dataset with one time step
            time = xr.DataArray(
                [_make_timestep(f"{year}-01-01")],
                dims=["time"],
                name="time",
            )
            ds = xr.Dataset({"time": time})
            ds.to_netcdf(filename)

        return filelist

    def _create_monthly_test_files(self, year_month_list):
        """
        Helper method to create test files with multiple timesteps per file (monthly)

        Args:
            year_month_list: List of tuples (year, list_of_months) where each file
                           contains multiple monthly timesteps

        Returns:
            List of filenames created
        """
        filelist = []
        for year, months in year_month_list:
            filename = os.path.join(self.temp_dir, f"test_{year}.nc")
            filelist.append(filename)

            # Create a dataset with multiple monthly time steps
            timesteps = [_make_timestep(f"{year}-{month:02d}-15 12:00:00") for month in months]
            time = xr.DataArray(timesteps, dims=["time"], name="time")
            ds = xr.Dataset({"time": time})
            ds.to_netcdf(filename)

        return filelist

    def _create_daily_test_file(self, year, month, days, *, hour=0, minute=0, second=0):
        """
        Helper method to create a test file with daily timesteps

        Args:
            year: Year for the file
            month: Month for the file
            days: List of days to include
            hour: Hour of day (default 0)
            minute: Minute of hour (default 0)
            second: Second of minute (default 0)

        Returns:
            Filename created
        """
        filename = os.path.join(self.temp_dir, f"test_{year}_{month:02d}.nc")

        # Create a dataset with daily time steps
        time_str = f"{hour:02d}:{minute:02d}:{second:02d}" if hour or minute or second else ""
        if time_str:
            timesteps = [_make_timestep(f"{year}-{month:02d}-{day:02d} {time_str}") for day in days]
        else:
            timesteps = [_make_timestep(f"{year}-{month:02d}-{day:02d}") for day in days]
        time = xr.DataArray(timesteps, dims=["time"], name="time")
        ds = xr.Dataset({"time": time})
        ds.to_netcdf(filename)

        return filename

    def test_get_files_in_time_slice_middle(self):
        """Test get_files_in_time_slice with a slice in the middle of the range"""
        years = [2000, 2001, 2002, 2003, 2004]
        filelist = self._create_annual_test_files(years)

        time_slice = slice(_make_timestep("2001-01-01"), _make_timestep("2003-01-01"))
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_2001.nc", "test_2002.nc", "test_2003.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_middle_slice_just_strings(self):
        """
        As test_get_files_in_time_slice_middle, but with the slice containing strings instead of
        actual cftime timestamps
        """
        years = [2000, 2001, 2002, 2003, 2004]
        filelist = self._create_annual_test_files(years)

        time_slice = slice("2001-01-01", "2003-01-01")
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_2001.nc", "test_2002.nc", "test_2003.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_from_beginning(self):
        """Test get_files_in_time_slice with unbounded start (from beginning)"""
        years = [2000, 2001, 2002, 2003, 2004]
        filelist = self._create_annual_test_files(years)

        time_slice = slice(None, _make_timestep("2001-01-01"))
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_2000.nc", "test_2001.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_to_end(self):
        """Test get_files_in_time_slice with unbounded end (to the end)"""
        years = [2000, 2001, 2002, 2003, 2004]
        filelist = self._create_annual_test_files(years)

        time_slice = slice(_make_timestep("2003-01-01"), None)
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_2003.nc", "test_2004.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_all_files(self):
        """Test get_files_in_time_slice with unbounded slice (all files)"""
        years = [2000, 2001, 2002, 2003, 2004]
        filelist = self._create_annual_test_files(years)

        time_slice = slice(None, None)
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        # For this test, compare full paths since expected is also full paths
        self.assertEqual(result, filelist)

    def test_get_files_in_time_slice_no_match(self):
        """Test get_files_in_time_slice with no matching files (should raise FileNotFoundError)"""
        years = [2000, 2001, 2002, 2003, 2004]
        filelist = self._create_annual_test_files(years)

        time_slice = slice(_make_timestep("2010-01-01"), _make_timestep("2011-01-01"))
        with self.assertRaises(FileNotFoundError):
            import_ds.get_files_in_time_slice(filelist, time_slice)

    def test_get_files_in_time_slice_monthly_multiple_per_file(self):
        """Test get_files_in_time_slice with monthly data, multiple timesteps per file"""
        # Create files with monthly data: each file has 12 months
        year_month_list = [
            (2000, list(range(1, 13))),  # Jan-Dec 2000
            (2001, list(range(1, 13))),  # Jan-Dec 2001
            (2002, list(range(1, 13))),  # Jan-Dec 2002
        ]
        filelist = self._create_monthly_test_files(year_month_list)

        # Select from mid-2000 to mid-2001
        time_slice = slice(_make_timestep("2000-06-01"), _make_timestep("2001-08-01"))
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_2000.nc", "test_2001.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_monthly_partial_overlap(self):
        """Test get_files_in_time_slice with monthly data selecting partial year"""
        year_month_list = [
            (2000, list(range(1, 13))),
            (2001, list(range(1, 13))),
            (2002, list(range(1, 13))),
        ]
        filelist = self._create_monthly_test_files(year_month_list)

        # Select only within 2001
        time_slice = slice(_make_timestep("2001-03-01"), _make_timestep("2001-09-01"))
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_2001.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_with_hours_minutes_seconds(self):
        """Test get_files_in_time_slice with timesteps including hours, minutes, seconds"""
        # Create daily files with specific times
        filelist = []
        filelist.append(
            self._create_daily_test_file(2000, 6, list(range(1, 31)), hour=6, minute=30, second=15)
        )
        filelist.append(
            self._create_daily_test_file(2000, 7, list(range(1, 32)), hour=6, minute=30, second=15)
        )
        filelist.append(
            self._create_daily_test_file(2000, 8, list(range(1, 32)), hour=6, minute=30, second=15)
        )

        # Select from mid-June to mid-July with specific time
        time_slice = slice(
            _make_timestep("2000-06-15 06:30:15"),
            _make_timestep("2000-07-20 06:30:15"),
        )
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_2000_06.nc", "test_2000_07.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_mixed_frequencies(self):
        """Test get_files_in_time_slice with files containing different numbers of timesteps"""
        filelist = []
        # File 1: Single timestep (annual)
        filename1 = os.path.join(self.temp_dir, "test_1999.nc")
        time1 = xr.DataArray([_make_timestep("1999-07-01 00:00:00")], dims=["time"], name="time")
        xr.Dataset({"time": time1}).to_netcdf(filename1)
        filelist.append(filename1)

        # File 2: Monthly timesteps
        filename2 = os.path.join(self.temp_dir, "test_2000.nc")
        timesteps2 = [_make_timestep(f"2000-{m:02d}-15 12:00:00") for m in range(1, 13)]
        time2 = xr.DataArray(timesteps2, dims=["time"], name="time")
        xr.Dataset({"time": time2}).to_netcdf(filename2)
        filelist.append(filename2)

        # File 3: Daily timesteps for one month
        filelist.append(
            self._create_daily_test_file(2001, 1, list(range(1, 32)), hour=3, minute=0, second=0)
        )

        # Select from late 1999 to mid-2000
        time_slice = slice(_make_timestep("1999-06-01"), _make_timestep("2000-08-01"))
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        expected = ["test_1999.nc", "test_2000.nc"]
        self.assertEqual(result_basenames, expected)

    def test_get_files_in_time_slice_exact_boundary_match(self):
        """Test get_files_in_time_slice with exact timestamp boundary matches at file edges"""
        filelist = []
        filelist.append(
            self._create_daily_test_file(
                2000, 12, list(range(1, 32)), hour=23, minute=59, second=59
            )
        )
        filelist.append(
            self._create_daily_test_file(2001, 1, list(range(1, 32)), hour=23, minute=59, second=59)
        )
        filelist.append(
            self._create_daily_test_file(2001, 2, list(range(1, 29)), hour=23, minute=59, second=59)
        )

        # Select from last timestep of December file to first timestep of February file
        # This tests exact boundary matching at file edges
        time_slice = slice(
            _make_timestep("2000-12-31 23:59:59"),  # Last timestep in first file
            _make_timestep("2001-02-01 23:59:59"),  # First timestep in third file
        )
        result = import_ds.get_files_in_time_slice(filelist, time_slice)
        result_basenames = [os.path.basename(f) for f in result]
        # Should include all three files since boundaries match exactly
        expected = ["test_2000_12.nc", "test_2001_01.nc", "test_2001_02.nc"]
        self.assertEqual(result_basenames, expected)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
