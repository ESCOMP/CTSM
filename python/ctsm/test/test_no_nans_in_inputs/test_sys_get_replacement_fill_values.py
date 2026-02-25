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

from ctsm.no_nans_in_inputs.constants import (
    ATTR,
    USER_REQ_SKIP_FILE,
    USER_REQ_QUIT,
)
from ctsm.no_nans_in_inputs.get_replacement_fill_values import (
    collect_new_fill_values,
    file_has_nan_fill,
    get_vars_with_nan_fills,
    main as main_func,
    show_ncdump_for_variable,
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


class TestGetVarsWithNanFills:
    """Test the get_vars_with_nan_fills function."""

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
        """Test that get_vars_with_nan_fills correctly detects NaN/non-NaN/absent fill values."""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        dtype = np.float32 if fill_value is None else type(fill_value)
        ds = xr.Dataset(
            {
                var_name: xr.DataArray(
                    np.array([1, 2], dtype=dtype),
                    dims=["time"],
                ),
            }
        )
        # Use encoding to set (or suppress) the _FillValue
        encoding = {var_name: {ATTR: fill_value}}
        ds.to_netcdf(str(test_file), encoding=encoding)
        vars_with_nan_fills = get_vars_with_nan_fills(str(test_file))
        assert (var_name in vars_with_nan_fills) == expected


class TestShowNcdumpForVariable:
    """Test the show_ncdump_for_variable function."""

    TEST_VAR = "temp"

    @pytest.fixture
    def test_netcdf_file(self, tmp_path):
        """Create a temporary NetCDF file for testing."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                self.TEST_VAR: xr.DataArray(
                    np.array([1.0, 2.0, 3.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan), "long_name": "temperature"},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds.close()
        return str(test_file)

    def test_prints_matching_lines(self, test_netcdf_file, capsys):
        """Test that matching lines from ncdump are printed."""
        show_ncdump_for_variable(test_netcdf_file, self.TEST_VAR)
        captured = capsys.readouterr()
        assert self.TEST_VAR in captured.out
        assert "Lines matching" in captured.out

    def test_prints_no_match_message(self, test_netcdf_file, capsys):
        """Test that a message is printed when no lines match."""
        show_ncdump_for_variable(test_netcdf_file, "nonexistent_var")
        captured = capsys.readouterr()
        assert "No lines found matching" in captured.out

    def test_none_file_path(self, capsys):
        """Test that None file path prints a message and returns."""
        show_ncdump_for_variable(None, self.TEST_VAR)
        captured = capsys.readouterr()
        assert "No file path available" in captured.out

    def test_nonexistent_file(self, capsys):
        """Test that a nonexistent file prints an error."""
        show_ncdump_for_variable("/nonexistent/file.nc", self.TEST_VAR)
        captured = capsys.readouterr()
        assert "Error running ncdump" in captured.out


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
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
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
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.get_vars_with_nan_fills",
        return_value=["temp"],
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
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
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.get_vars_with_nan_fills",
        return_value=["temp"],
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
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
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
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


class TestCollectNewFillValues:
    """Test the collect_new_fill_values function using real files."""

    def _create_test_netcdf(self, path, var_dict):
        """Helper function to create a NetCDF file for testing."""
        ds = xr.Dataset()
        for var_name, properties in var_dict.items():
            ds[var_name] = xr.DataArray(
                properties["data"],
                dims=["time"],
                attrs=properties["attrs"],
            )
        ds.to_netcdf(path)
        ds.close()

    @patch("ctsm.no_nans_in_inputs.json_io.save_progress")
    @patch("builtins.input", side_effect=["-999"])
    def test_user_provides_number(self, mock_input, mock_save, tmp_path, capsys):
        """Test happy path where user provides a numeric fill value."""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        self._create_test_netcdf(
            test_file,
            {
                var_name: {
                    "data": np.array([1.0, 2.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan), "long_name": "temperature"},
                }
            },
        )
        expected_fill_value = -999.0
        expected_result = {str(test_file): create_empty_progress_dict_onefile()}
        assert not expected_result[str(test_file)]["new_fill_values"]

        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        result = collect_new_fill_values(progress)

        # Check final result
        expected_result[str(test_file)]["new_fill_values"] = {var_name: expected_fill_value}
        expected_result[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert result == expected_result

        # Check what was saved
        mock_save.assert_called_once()
        saved_data = mock_save.call_args[0][0]
        saved_fill_value = saved_data[str(test_file)]["new_fill_values"][var_name]
        assert saved_fill_value == expected_fill_value
        assert isinstance(saved_fill_value, float)

    @patch("ctsm.no_nans_in_inputs.json_io.save_progress")
    @patch("builtins.input", side_effect=[USER_REQ_SKIP_FILE])
    def test_user_skips_file(self, mock_input, mock_save, tmp_path):
        """Test that requesting skipfile skips the current file."""
        test_file = tmp_path / "test.nc"
        self._create_test_netcdf(
            test_file,
            {
                "temp": {
                    "data": np.array([1.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan)},
                },
                "precip": {
                    "data": np.array([2.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan)},
                },
            },
        )
        progress = {str(test_file): create_empty_progress_dict_onefile()}
        result = collect_new_fill_values(progress)

        # 'temp' is processed, user skips, neither it nor 'precip' are in results
        expected_result = {str(test_file): create_empty_progress_dict_onefile()}
        assert not expected_result[str(test_file)]["new_fill_values"]
        assert result == expected_result
        # Progress was never saved
        mock_save.assert_not_called()

    @patch("ctsm.no_nans_in_inputs.json_io.save_progress")
    @patch("builtins.input", side_effect=["-100", USER_REQ_SKIP_FILE])
    def test_user_enters_then_skips_file(self, mock_input, mock_save, tmp_path):
        """Test that requesting skipfile skips the current file but saves entered results."""
        test_file = tmp_path / "test.nc"
        data = {
            "temp": {
                "data": np.array([1.0], dtype=np.float32),
                "attrs": {ATTR: np.float32(np.nan)},
            },
            "precip": {
                "data": np.array([2.0], dtype=np.float32),
                "attrs": {ATTR: np.float32(np.nan)},
            },
        }
        self._create_test_netcdf(
            test_file,
            data,
        )
        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = data.keys()
        result = collect_new_fill_values(progress)

        # 'temp' is processed, user enters so it saves, then precip is skipped
        expected_result = {str(test_file): create_empty_progress_dict_onefile()}
        expected_result[str(test_file)]["new_fill_values"] = {"temp": -100.0}
        expected_result[str(test_file)]["vars_with_nan_fills"] = data.keys()
        assert result == expected_result
        # Progress was saved only after temp
        mock_save.assert_called_once()

    @patch("ctsm.no_nans_in_inputs.json_io.save_progress")
    @patch("builtins.input", side_effect=[USER_REQ_QUIT])
    def test_user_quits(self, mock_input, mock_save, tmp_path, mock_progress_file):
        """Test that requesting quit exits the collection loop."""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        self._create_test_netcdf(
            test_file,
            {var_name: {"data": [1.0], "attrs": {ATTR: np.nan}}},
        )
        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert not progress[str(test_file)]["new_fill_values"]

        assert not os.path.exists(str(mock_progress_file))

        with pytest.raises(SystemExit) as e:
            collect_new_fill_values(progress)
        assert e.value.code == 0

        # No value was collected, so save_progress is not called
        mock_save.assert_not_called()

    @patch("ctsm.no_nans_in_inputs.json_io.save_progress")
    def test_dry_run(self, mock_save, tmp_path, capsys):
        """Test --dry-run"""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        self._create_test_netcdf(
            test_file,
            {
                var_name: {
                    "data": np.array([1.0, 2.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan), "long_name": "temperature"},
                }
            },
        )

        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        result = collect_new_fill_values(progress, dry_run=True)

        # Check that nothing was saved
        mock_save.assert_not_called()

        # Check that result is what we expect
        expected = {str(test_file): create_empty_progress_dict_onefile()}
        expected[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert result == expected

        # Check that result does print some info but not everything
        stdout = capsys.readouterr().out
        assert "Variable: temp" in stdout
        assert "new fill value(s)" not in stdout
