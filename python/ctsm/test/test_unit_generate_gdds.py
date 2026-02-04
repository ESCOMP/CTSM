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
from cftime import DatetimeNoLeap

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


class TestGetTimeSliceList(unittest.TestCase):
    """Tests for _get_time_slice_list()"""

    def test_generate_gdds_get_time_slice_list(self):
        """Test that _get_time_slice_list works with two different years"""
        season_list = [1986, 1987]
        result = gg._get_time_slice_list(season_list[0], season_list[-1])
        expected = [
            slice("1987-01-01", "1987-12-31"),
            slice("1988-01-01", "1988-12-31"),
            slice("1989-01-01", "1989-12-31"),
        ]
        assert result == expected

    def test_generate_gdds_get_time_slice_list_1yr(self):
        """Test that _get_time_slice_list works with the same year"""
        result = gg._get_time_slice_list(1987, 1987)
        expected = [
            slice("1988-01-01", "1988-12-31"),
            slice("1989-01-01", "1989-12-31"),
        ]
        assert result == expected

    def test_generate_gdds_get_time_slice_list_valueerror(self):
        """Test that _get_time_slice_list raises ValueError if last < first"""
        with self.assertRaisesRegex(ValueError, "first_season.* > last_season"):
            gg._get_time_slice_list(1987, 1986)

    def test_generate_gdds_get_time_slice_list_typeerror_first(self):
        """Test that _get_time_slice_list raises TypeError if not given integer first season"""
        with self.assertRaisesRegex(
            TypeError, r"_get_time_slice_list\(\) arguments must be integers"
        ):
            gg._get_time_slice_list(1986.3, 1987)

    def test_generate_gdds_get_time_slice_list_typeerror_last(self):
        """Test that _get_time_slice_list raises TypeError if not given integer last season"""
        with self.assertRaisesRegex(
            TypeError, r"_get_time_slice_list\(\) arguments must be integers"
        ):
            gg._get_time_slice_list(1986, None)


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

    def test_find_inst_hist_files_h1_no_year(self):
        """Test finding h1 files without specifying year"""
        # Create test files
        file1 = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        file2 = self._create_test_file("test.clm2.h1i.2000-02-01-00000.nc")
        file3 = self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc")

        result = gf.find_inst_hist_files(self.temp_dir, h=1, this_year=None)

        # Should find all h1i files
        self.assertEqual(len(result), 3)
        self.assertIn(file1, result)
        self.assertIn(file2, result)
        self.assertIn(file3, result)

    def test_find_inst_hist_files_h2_no_year(self):
        """Test finding h2 files without specifying year"""
        # Create test files
        file1 = self._create_test_file("test.clm2.h2i.2000-01-01-00000.nc")
        file2 = self._create_test_file("test.clm2.h2i.2001-01-01-00000.nc")
        # Create h1 file that should not be found
        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")

        result = gf.find_inst_hist_files(self.temp_dir, h=2, this_year=None)

        # Should find only h2i files
        self.assertEqual(len(result), 2)
        self.assertIn(file1, result)
        self.assertIn(file2, result)

    def test_find_inst_hist_files_with_year(self):
        """Test finding files for a specific year"""
        # Create test files
        file_2000 = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        file_2001 = self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc")
        file_2002 = self._create_test_file("test.clm2.h1i.2002-01-01-00000.nc")

        result = gf.find_inst_hist_files(self.temp_dir, h=1, this_year=2001)

        # Should find only 2001 file
        self.assertEqual(len(result), 1)
        self.assertIn(file_2001, result)
        self.assertNotIn(file_2000, result)
        self.assertNotIn(file_2002, result)

    def test_find_inst_hist_files_base_extension(self):
        """Test finding files with .nc.base extension"""
        # Create test files with .nc.base extension
        file1 = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc.base")
        file2 = self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc.base")

        result = gf.find_inst_hist_files(self.temp_dir, h=1, this_year=None)

        # Should find .nc.base files
        self.assertEqual(len(result), 2)
        self.assertIn(file1, result)
        self.assertIn(file2, result)

    def test_find_inst_hist_files_prefer_nc_over_base(self):
        """Test that .nc files are preferred over .nc.base files"""
        # Create both .nc and .nc.base files
        file_nc = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        file_nc_base = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc.base")

        result = gf.find_inst_hist_files(self.temp_dir, h=1, this_year=None)

        # Should find .nc files first (pattern order preference)
        self.assertIn(file_nc, result)
        self.assertNotIn(file_nc_base, result)

    def test_find_inst_hist_files_multiple_months_same_year(self):
        """Test finding multiple files from the same year"""
        # Create multiple files from 2000
        file1 = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        file2 = self._create_test_file("test.clm2.h1i.2000-01-15-00000.nc")
        file3 = self._create_test_file("test.clm2.h1i.2000-01-31-00000.nc")
        # Create file from different year
        self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc")

        result = gf.find_inst_hist_files(self.temp_dir, h=1, this_year=2000)

        # Should find all January 2000 files
        self.assertEqual(len(result), 3)
        self.assertIn(file1, result)
        self.assertIn(file2, result)
        self.assertIn(file3, result)

    def test_find_inst_hist_files_no_files_found(self):
        """Test error when no matching files are found"""
        # Create a non-matching file
        self._create_test_file("test.clm2.h0.2000-01-01-00000.nc")

        # Should raise a FileNotFoundError error
        with self.assertRaises(FileNotFoundError):
            gf.find_inst_hist_files(self.temp_dir, h=1, this_year=None)

    def test_find_inst_hist_files_different_case_names(self):
        """Test that RuntimeError is raised when files from different case names are found"""
        # Create files with different case names
        self._create_test_file("case1.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("case2.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("longcasename.clm2.h1i.2000-01-01-00000.nc")

        # Should raise RuntimeError due to multiple case names
        with self.assertRaises(RuntimeError):
            gf.find_inst_hist_files(self.temp_dir, h=1, this_year=2000)

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
        with self.assertRaises(RuntimeError):
            gf.find_inst_hist_files(self.temp_dir, h=1, this_year=2000, logger=logger)

    def test_find_inst_hist_files_no_files_found_with_logger(self):
        """Test error when no matching files are found, with logger"""
        # Create a logger
        logger = logging.getLogger("test_logger_no_files")
        logger.setLevel(logging.DEBUG)

        # Create a non-matching file
        self._create_test_file("test.clm2.h0.2000-01-01-00000.nc")

        # Should raise a FileNotFoundError even with logger
        with self.assertRaises(FileNotFoundError):
            gf.find_inst_hist_files(self.temp_dir, h=1, this_year=None, logger=logger)

    def test_find_inst_hist_files_h_str_with_logger(self):
        """Test that TypeError is raised when h is a string, with logger"""
        # Create a logger
        logger = logging.getLogger("test_logger_h_str")
        logger.setLevel(logging.DEBUG)

        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")

        with self.assertRaises(TypeError):
            gf.find_inst_hist_files(self.temp_dir, h="1", this_year=2000, logger=logger)

    def test_find_inst_hist_files_h_float_with_logger(self):
        """Test that TypeError is raised when h is a float, with logger"""
        # Create a logger
        logger = logging.getLogger("test_logger_h_float")
        logger.setLevel(logging.DEBUG)

        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")

        with self.assertRaises(TypeError):
            gf.find_inst_hist_files(self.temp_dir, h=1.0, this_year=2000, logger=logger)


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
        # Create h1 and h2 files for 2000
        h1_file = self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        h2_file = self._create_test_file("test.clm2.h2i.2000-01-01-00000.nc")

        time_slice_list = [slice("2000-01-01", "2000-12-31")]

        h1_file_lists, h2_file_lists = gg._get_file_lists(
            self.temp_dir, time_slice_list, logger=None
        )

        # Should have one list for each time slice
        self.assertEqual(len(h1_file_lists), 1)
        self.assertEqual(len(h2_file_lists), 1)

        # Check contents of file lists
        # pylint: disable=unsubscriptable-object
        self.assertEqual(len(h1_file_lists[0]), 1)
        self.assertEqual(len(h2_file_lists[0]), 1)
        self.assertEqual(h1_file_lists[0], [h1_file])
        self.assertEqual(h2_file_lists[0], [h2_file])

    def test_get_file_lists_multiple_years(self):
        """Test _get_file_lists with multiple years of data"""
        # Create h1 and h2 files for 2000-2002
        h1_files = []
        h2_files = []
        for year in [2000, 2001, 2002]:
            h1_files.append(self._create_test_file(f"test.clm2.h1i.{year}-01-01-00000.nc"))
            h2_files.append(self._create_test_file(f"test.clm2.h2i.{year}-01-01-00000.nc"))

        time_slice_list = [
            slice("2000-01-01", "2000-12-31"),
            slice("2001-01-01", "2001-12-31"),
            slice("2002-01-01", "2002-12-31"),
        ]

        h1_file_lists, h2_file_lists = gg._get_file_lists(
            self.temp_dir, time_slice_list, logger=None
        )

        # Should have one list for each time slice
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
        # Create multiple h1 and h2 files for 2000
        h1_files = []
        h2_files = []
        for month in ["01", "06", "12"]:
            h1_files.append(self._create_test_file(f"test.clm2.h1i.2000-{month}-01-00000.nc"))
            h2_files.append(self._create_test_file(f"test.clm2.h2i.2000-{month}-01-00000.nc"))

        time_slice_list = [slice("2000-01-01", "2000-12-31")]

        h1_file_lists, h2_file_lists = gg._get_file_lists(
            self.temp_dir, time_slice_list, logger=None
        )

        # Should have one list for the time slice
        self.assertEqual(len(h1_file_lists), 1)
        self.assertEqual(len(h2_file_lists), 1)

        # Check contents of file lists (should be sorted)
        # pylint: disable=unsubscriptable-object
        self.assertEqual(len(h1_file_lists[0]), 3)
        self.assertEqual(len(h2_file_lists[0]), 3)
        self.assertEqual(h1_file_lists[0], sorted(h1_files))
        self.assertEqual(h2_file_lists[0], sorted(h2_files))

    def test_get_file_lists_no_h1_files(self):
        """Test _get_file_lists when h1 files are missing"""
        # Create only h2 files
        self._create_test_file("test.clm2.h2i.2000-01-01-00000.nc")

        time_slice_list = [slice("2000-01-01", "2000-12-31")]

        # Should raise FileNotFoundError when h1 files are not found
        with self.assertRaises(FileNotFoundError):
            gg._get_file_lists(self.temp_dir, time_slice_list, logger=None)

    def test_get_file_lists_no_h2_files(self):
        """Test _get_file_lists when h2 files are missing"""
        # Create only h1 files
        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")

        time_slice_list = [slice("2000-01-01", "2000-12-31")]

        # Should raise FileNotFoundError when h2 files are not found
        with self.assertRaises(FileNotFoundError):
            gg._get_file_lists(self.temp_dir, time_slice_list, logger=None)

    def test_get_file_lists_h1_outside_time_slice(self):
        """Test _get_file_lists when h1 files exist but have no timesteps in the slice"""
        # Create h1 files for 2000 and h2 files for 2001
        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("test.clm2.h2i.2001-01-01-00000.nc")

        # Request time slice for 2001 (h1 files exist but are outside the slice)
        time_slice_list = [slice("2001-01-01", "2001-12-31")]

        # Should raise FileNotFoundError when h1 files have no timesteps in slice
        with self.assertRaisesRegex(FileNotFoundError, "h1"):
            gg._get_file_lists(self.temp_dir, time_slice_list, logger=None)

    def test_get_file_lists_h2_outside_time_slice(self):
        """Test _get_file_lists when h2 files exist but have no timesteps in the slice"""
        # Create h1 files for 2001 and h2 files for 2000
        self._create_test_file("test.clm2.h1i.2001-01-01-00000.nc")
        self._create_test_file("test.clm2.h2i.2000-01-01-00000.nc")

        # Request time slice for 2001 (h2 files exist but are outside the slice)
        time_slice_list = [slice("2001-01-01", "2001-12-31")]

        # Should raise FileNotFoundError when h2 files have no timesteps in slice
        with self.assertRaisesRegex(FileNotFoundError, "h2"):
            gg._get_file_lists(self.temp_dir, time_slice_list, logger=None)

    def test_get_file_lists_partial_overlap(self):
        """Test _get_file_lists when some time slices have files and others don't"""
        # Create h1 and h2 files for 2000 only
        self._create_test_file("test.clm2.h1i.2000-01-01-00000.nc")
        self._create_test_file("test.clm2.h2i.2000-01-01-00000.nc")

        # Request time slices for 2000 and 2001
        time_slice_list = [
            slice("2000-01-01", "2000-12-31"),
            slice("2001-01-01", "2001-12-31"),
        ]

        # Should raise FileNotFoundError when second time slice has no files
        with self.assertRaises(FileNotFoundError):
            gg._get_file_lists(self.temp_dir, time_slice_list, logger=None)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
