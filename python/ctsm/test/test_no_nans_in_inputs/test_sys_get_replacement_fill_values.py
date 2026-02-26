#!/usr/bin/env python3
"""
System tests for get_replacement_fill_values.py script.

Tests functions that require file I/O.
"""
# pylint: disable=too-few-public-methods

import os
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR
from ctsm.no_nans_in_inputs.get_replacement_fill_values import (
    file_has_nan_fill,
    main as main_func,
)

from ctsm.no_nans_in_inputs.json_io import create_empty_progress_dict_onefile


class TestFileHasNanFill:
    """Test the file_has_nan_fill function."""

    @pytest.mark.parametrize(
        "fill_value, expected",
        [
            (np.float32(np.nan), True),
            (np.float64(np.nan), True),
            (np.float32(-999.0), False),
            (np.int32(-999), False),
            (None, False),
        ],
    )
    def test_fill_value_detection(self, tmp_path, fill_value, expected):
        """Test that file_has_nan_fill correctly detects NaN vs non-NaN vs absent fill values."""
        test_file = tmp_path / "test.nc"
        dtype = np.float32 if fill_value is None else type(fill_value)
        ds = xr.Dataset(
            {
                "temp": xr.DataArray(
                    np.array([1, 2], dtype=dtype),
                    dims=["time"],
                ),
            }
        )
        # Use encoding to set (or suppress) the _FillValue
        encoding = {"temp": {ATTR: fill_value}}
        ds.to_netcdf(str(test_file), encoding=encoding)
        any_nan_fill, vars_with_nan_fills = file_has_nan_fill(str(test_file))
        assert any_nan_fill == expected
        assert len(vars_with_nan_fills) == int(expected)


class TestMain:
    """Test the main function of get_replacement_fill_values.py."""

    @patch("sys.argv", ["get_replacement_fill_values.py"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.namelist_utils._extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch(
        "ctsm.no_nans_in_inputs.namelist_utils._extract_file_path_set_from_usernl",
        return_value=set(),
    )
    @patch("os.path.exists", return_value=True)
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.file_has_nan_fill",
        return_value=(True, ["temp"]),
    )
    @patch("ctsm.no_nans_in_inputs.user_inputs.collect_new_fill_values")
    @patch(
        "ctsm.no_nans_in_inputs.namelist_utils.find_user_nl_files",
        return_value=["user_nl_clm"],
    )
    def test_main_happy_path(
        self,
        mock_find_user_nl_files,
        mock_collect,
        mock_file_has_nan,
        mock_exists,
        mock_extract_from_usernl,
        mock_extract_from_xml,
        mock_check_write,
        tmp_path,
        mock_xml_file_path,
    ):  # pylint: disable=unused-argument
        """Test main function with a single matching file."""
        # Setup mock dataset
        mock_ds = MagicMock(spec=xr.Dataset)
        # Make the mock dataset iterable
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()

        assert result == 0
        assert mock_extract_from_xml.call_count == 2
        assert mock_extract_from_usernl.call_count == 1
        path_to_test_file_rel = "lnd/clm2/test.nc"
        path_to_test_file_abs = os.path.join(str(tmp_path), path_to_test_file_rel)
        mock_collect.assert_called_once()
        mock_check_write.assert_called_once()
        mock_find_user_nl_files.assert_called_once()

        # Check the arguments passed to collect_new_fill_values
        args, kwargs = mock_collect.call_args
        progress = args[0]
        assert len(progress) == 1
        expected_dict = create_empty_progress_dict_onefile()
        expected_dict["found_in_files"] = {mock_xml_file_path: {path_to_test_file_rel}}
        expected_dict["vars_with_nan_fills"] = ["temp"]
        assert progress[path_to_test_file_abs] == expected_dict
        assert "delete_if_none_filled" in kwargs
        assert not kwargs["delete_if_none_filled"]
        assert "dry_run" in kwargs
        assert not kwargs["dry_run"]

    @patch("sys.argv", ["get_replacement_fill_values.py", "--dry-run"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.namelist_utils._extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch("os.path.exists", return_value=True)
    @patch(
        "ctsm.no_nans_in_inputs.netcdf_utils.get_vars_with_nan_fills",
        return_value=["temp"],
    )
    @patch("ctsm.no_nans_in_inputs.user_inputs.collect_new_fill_values")
    @patch("ctsm.no_nans_in_inputs.namelist_utils.find_user_nl_files", return_value=[])
    def test_main_dry_run(
        self,
        mock_find_user_nl_files,
        mock_collect,
        mock_get_vars_with_nan,
        mock_exists,
        mock_extract,
        mock_check_write,
        capsys,
        tmp_path,
        mock_xml_file_path,
    ):  # pylint: disable=unused-argument
        """Test main function with a single matching file under --dry-run"""
        # Setup mock dataset
        mock_ds = MagicMock(spec=xr.Dataset)
        # Make the mock dataset iterable
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()

        assert result == 0
        assert mock_extract.call_count == 2
        path_to_test_file_rel = "lnd/clm2/test.nc"
        path_to_test_file_abs = os.path.join(str(tmp_path), path_to_test_file_rel)
        mock_get_vars_with_nan.assert_called_with(path_to_test_file_abs)
        mock_collect.assert_called_once()
        mock_find_user_nl_files.assert_called_once()

        # Shouldn't check write access in dry run
        mock_check_write.assert_not_called()

        # Check the arguments passed to collect_new_fill_values
        args, kwargs = mock_collect.call_args
        progress = args[0]
        expected_dict = create_empty_progress_dict_onefile()
        expected_dict["found_in_files"] = {mock_xml_file_path: {path_to_test_file_rel}}
        expected_dict["vars_with_nan_fills"] = ["temp"]
        assert len(progress) == 1
        assert progress[path_to_test_file_abs] == expected_dict
        assert "delete_if_none_filled" in kwargs
        assert not kwargs["delete_if_none_filled"]
        assert "dry_run" in kwargs
        assert kwargs["dry_run"]

        # Check stdout
        stdout = capsys.readouterr().out
        assert "Checking write access" not in stdout

    @patch("sys.argv", ["get_replacement_fill_values.py"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.namelist_utils._extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch(
        "ctsm.no_nans_in_inputs.netcdf_utils.get_vars_with_nan_fills",
        return_value=["temp"],
    )
    @patch("ctsm.no_nans_in_inputs.user_inputs.collect_new_fill_values")
    @patch("ctsm.no_nans_in_inputs.namelist_utils.find_user_nl_files", return_value=[])
    def test_main_missing_file(
        self,
        mock_find_user_nl_files,
        mock_collect,
        mock_get_vars_with_nan,
        mock_extract,
        mock_check_write,
        capsys,
        tmp_path,
        mock_xml_file_path,
    ):  # pylint: disable=unused-argument
        """Test main function with a file that doesn't exist"""
        # Setup mock dataset
        mock_ds = MagicMock(spec=xr.Dataset)
        # Make the mock dataset iterable
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()

        assert result == 0
        mock_extract.assert_called_once()
        path_to_test_file_rel = "lnd/clm2/test.nc"
        path_to_test_file_abs = os.path.join(str(tmp_path), path_to_test_file_rel)
        mock_get_vars_with_nan.assert_not_called()

        # Check the arguments passed to collect_new_fill_values
        mock_collect.assert_called_once()
        args, _ = mock_collect.call_args
        progress = args[0]
        assert progress == {}

        # Check stdout
        stdout = capsys.readouterr().out
        assert "1\tFiles not found" in stdout
        assert f"Not found: '{path_to_test_file_abs}'" in stdout

    @patch("sys.argv", ["get_replacement_fill_values.py", "--delete-if-none-filled"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.namelist_utils._extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch("os.path.exists", return_value=True)
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.file_has_nan_fill",
        return_value=(True, ["temp"]),
    )
    @patch("ctsm.no_nans_in_inputs.user_inputs.collect_new_fill_values")
    def test_main_with_delete_flag(
        self,
        mock_collect,
        mock_file_has_nan,
        mock_exists,
        mock_extract,
        mock_check_write,
    ):  # pylint: disable=unused-argument
        """Test main function with --delete-if-none-filled flag."""
        mock_ds = MagicMock(spec=xr.Dataset)
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()
        assert result == 0

        mock_collect.assert_called_once()
        _, kwargs = mock_collect.call_args
        assert "delete_if_none_filled" in kwargs
        assert kwargs["delete_if_none_filled"]

    @patch("sys.argv", ["get_replacement_fill_values.py"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=False,
    )
    def test_main_no_write_access_exits(
        self, mock_check_write, capsys
    ):  # pylint: disable=unused-argument
        """Test that main exits if there's no write access."""
        with pytest.raises(SystemExit) as e:
            main_func()
        assert e.value.code == 1
        captured = capsys.readouterr()
        assert "No write access" in captured.err
