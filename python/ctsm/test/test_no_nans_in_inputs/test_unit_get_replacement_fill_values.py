#!/usr/bin/env python3
"""
Unit tests for get_replacement_fill_values.py script.
"""

from unittest.mock import MagicMock

import pytest

from ctsm.no_nans_in_inputs import namelist_utils



class TestHowNetcdfIsReferencedInFile:
    """Tests of how_netcdf_is_referenced_in_file()."""

    def return_input(self, x):
        """Take one input argument and return it"""
        return x

    @pytest.fixture(autouse=True)
    def mock_convert_to_absolute_path(self, monkeypatch):
        """Mock convert_to_absolute_path() to just return what it was given"""
        mock = MagicMock(side_effect=lambda x, *args, **kwargs: x)
        monkeypatch.setattr(namelist_utils, "convert_to_absolute_path", mock)
        return mock

    @pytest.fixture(autouse=True)
    def mock_replace_env_vars_in_netcdf_paths(self, monkeypatch):
        """Mock _replace_env_vars_in_netcdf_paths() to just return what it was given"""
        mock = MagicMock(side_effect=lambda x, *args, **kwargs: x)
        monkeypatch.setattr(namelist_utils, "_replace_env_vars_in_netcdf_paths", mock)
        return mock

    def test_how_netcdf_is_referenced_in_file_1found_once(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """Test how_netcdf_is_referenced_in_file() for one netCDF file present once in text file"""
        nc_file = "file.nc"
        monkeypatch.setattr(
            namelist_utils,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: [nc_file],
        )
        set_of_how_this_netcdf_appears = namelist_utils.how_netcdf_is_referenced_in_file(
            "dummy", nc_file
        )
        assert mock_convert_to_absolute_path.call_count == 2
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 1
        assert set_of_how_this_netcdf_appears == {nc_file}

    def test_how_netcdf_is_referenced_in_file_1found_twice_same(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """
        Test how_netcdf_is_referenced_in_file() for one netCDF file present twice in text file in
        the exact same way
        """
        nc_file = "file.nc"
        monkeypatch.setattr(
            namelist_utils,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: [nc_file, nc_file],
        )
        set_of_how_this_netcdf_appears = namelist_utils.how_netcdf_is_referenced_in_file(
            "dummy", nc_file
        )
        assert mock_convert_to_absolute_path.call_count == 3
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 2
        assert set_of_how_this_netcdf_appears == {nc_file}

    def test_how_netcdf_is_referenced_in_file_1found_twice_diff(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """
        Test how_netcdf_is_referenced_in_file() for one netCDF file present twice in text file in
        different ways
        """
        nc_file = "file.nc"
        nc_file2 = "abc123" + nc_file
        monkeypatch.setattr(
            namelist_utils,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: [nc_file, nc_file2],
        )
        set_of_how_this_netcdf_appears = namelist_utils.how_netcdf_is_referenced_in_file(
            "dummy", nc_file
        )
        assert mock_convert_to_absolute_path.call_count == 3
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 2
        assert set_of_how_this_netcdf_appears == {nc_file}

    def test_how_netcdf_is_referenced_in_file_2found(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """Test how_netcdf_is_referenced_in_file() for two netCDF files present in text file"""
        nc_file = "file.nc"
        nc_files = [nc_file, "file2.nc"]
        monkeypatch.setattr(
            namelist_utils,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: nc_files,
        )
        set_of_how_this_netcdf_appears = namelist_utils.how_netcdf_is_referenced_in_file(
            "dummy", nc_file
        )
        assert mock_convert_to_absolute_path.call_count == 3
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 2
        assert set_of_how_this_netcdf_appears == {nc_file}
