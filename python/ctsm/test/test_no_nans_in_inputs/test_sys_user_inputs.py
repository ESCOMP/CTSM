"""System tests of user_inputs (anything touching filesystem)"""

import os
from unittest.mock import patch

import pytest
import numpy as np
import xarray as xr

from ctsm.no_nans_in_inputs.constants import (
    ATTR,
    USER_REQ_QUIT,
    USER_REQ_SKIP_FILE,
)
from ctsm.no_nans_in_inputs.user_inputs import collect_new_fill_values
from ctsm.no_nans_in_inputs.json_io import NoNanFillValueProgress


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

    @patch("builtins.input", side_effect=["-999"])
    def test_user_provides_number(
        self, mock_input, tmp_path, capsys
    ):  # pylint: disable=unused-argument
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
        expected_result = NoNanFillValueProgress(progress_file=None)
        assert not expected_result[str(test_file)]["new_fill_values"]

        progress_file = tmp_path / "progress.json"
        progress = NoNanFillValueProgress(progress_file=progress_file)
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        result = collect_new_fill_values(progress)

        # Check final result
        expected_result[str(test_file)]["new_fill_values"] = {var_name: expected_fill_value}
        expected_result[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert result == expected_result

        # Check what was saved
        saved_data = NoNanFillValueProgress(progress_file=progress_file, load_without_asking=True)
        d = saved_data[str(test_file)]["new_fill_values"]
        saved_fill_value = d[var_name]
        assert saved_fill_value == expected_fill_value
        assert isinstance(saved_fill_value, float)

    @patch("ctsm.no_nans_in_inputs.json_io.NoNanFillValueProgress.save")
    @patch("builtins.input", side_effect=[USER_REQ_SKIP_FILE])
    def test_user_skips_file(
        self, mock_input, mock_save, tmp_path
    ):  # pylint: disable=unused-argument
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
        progress_file = tmp_path / "progress.json"
        progress = NoNanFillValueProgress(progress_file=progress_file)
        result = collect_new_fill_values(progress)

        # 'temp' is processed, user skips, neither it nor 'precip' are in results
        expected_result = NoNanFillValueProgress(progress_file=None)
        assert result == expected_result
        # Progress was never saved
        mock_save.assert_not_called()

    @patch("ctsm.no_nans_in_inputs.json_io.NoNanFillValueProgress.save")
    @patch("builtins.input", side_effect=["-100", USER_REQ_SKIP_FILE])
    def test_user_enters_then_skips_file(
        self, mock_input, mock_save, tmp_path
    ):  # pylint: disable=unused-argument
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
        progress = NoNanFillValueProgress(progress_file=None)
        progress[str(test_file)]["vars_with_nan_fills"] = data.keys()
        result = collect_new_fill_values(progress)

        # 'temp' is processed, user enters so it saves, then precip is skipped
        expected_result = NoNanFillValueProgress(progress_file=None)
        expected_result[str(test_file)]["new_fill_values"] = {"temp": -100.0}
        expected_result[str(test_file)]["vars_with_nan_fills"] = data.keys()
        assert result == expected_result
        # Progress was saved only after temp
        mock_save.assert_called_once()

    @patch("ctsm.no_nans_in_inputs.json_io.NoNanFillValueProgress.save")
    @patch("builtins.input", side_effect=[USER_REQ_QUIT])
    def test_user_quits(
        self, mock_input, mock_save, tmp_path, mock_progress_file
    ):  # pylint: disable=unused-argument
        """Test that requesting quit exits the collection loop."""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        self._create_test_netcdf(
            test_file,
            {var_name: {"data": [1.0], "attrs": {ATTR: np.nan}}},
        )
        progress = NoNanFillValueProgress(progress_file=None)
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert not progress[str(test_file)]["new_fill_values"]

        assert not os.path.exists(str(mock_progress_file))

        with pytest.raises(SystemExit) as e:
            collect_new_fill_values(progress)
        assert e.value.code == 0

        # No value was collected, so save_progress is not called
        mock_save.assert_not_called()

    @patch("ctsm.no_nans_in_inputs.json_io.NoNanFillValueProgress.save")
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

        progress = NoNanFillValueProgress(progress_file=None)
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        result = collect_new_fill_values(progress, dry_run=True)

        # Check that nothing was saved
        mock_save.assert_not_called()

        # Check that result is what we expect
        expected = NoNanFillValueProgress(progress_file=None)
        expected[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert result == expected

        # Check that result does print some info but not everything
        stdout = capsys.readouterr().out
        assert "Variable: temp" in stdout
        assert "new fill value(s)" not in stdout
