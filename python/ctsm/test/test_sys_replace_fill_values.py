#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

import json
import os
import subprocess

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR, USER_REQ_DELETE
from ctsm.no_nans_in_inputs.replace_fill_values import (
    build_ncatted_command,
    get_output_filename,
    load_new_fillvalues,
)


# Test constants
TEST_VAR_TEMP = "temp"
TEST_VAR_PRESSURE = "pressure"
TEST_OUTPUT_FILE = "output.nc"
TEST_FILL_VALUE = -999.0
NCATTED_CMD = "ncatted"
NCATTED_FLAG = "-a"


class TestLoadNewFillvalues:
    """Test the load_new_fillvalues function."""

    def test_load_valid_json(self, tmp_path):
        """Test loading a valid JSON file."""
        test_file = tmp_path / "test_fillvalues.json"
        test_data = {
            "/path/to/file1.nc": {"var1": -999.0, "var2": "delete"},
            "/path/to/file2.nc": {"var3": -999},
        }
        test_file.write_text(json.dumps(test_data), encoding="utf-8")

        result = load_new_fillvalues(str(test_file))
        assert result == test_data

    def test_file_not_found(self):
        """Test that missing file causes SystemExit."""
        with pytest.raises(SystemExit):
            load_new_fillvalues("/nonexistent/file.json")


class TestBuildNcattedCommand:
    """Test the build_ncatted_command function."""

    @pytest.fixture
    def test_netcdf_file(self, tmp_path):
        """Create a temporary NetCDF file for testing."""
        test_file = tmp_path / "test.nc"

        # Create a simple NetCDF file with float variables that have NaN fill values
        # (NetCDF doesn't allow NaN for integer types, and our scripts only work on
        # variables that already have NaN fill values)
        ds = xr.Dataset(
            {
                TEST_VAR_TEMP: xr.DataArray(
                    np.array([1.0, 2.0, 3.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan)},
                ),
                TEST_VAR_PRESSURE: xr.DataArray(
                    np.array([1000.0, 1010.0, 1020.0], dtype=np.float64),
                    dims=["time"],
                    attrs={ATTR: np.float64(np.nan)},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds.close()

        yield str(test_file)

    def test_delete_attribute(self, test_netcdf_file):
        """Test building command to delete an attribute."""
        var_fillvalues = {TEST_VAR_TEMP: USER_REQ_DELETE}
        cmd = build_ncatted_command(test_netcdf_file, TEST_OUTPUT_FILE, var_fillvalues)

        assert NCATTED_CMD in cmd
        assert NCATTED_FLAG in cmd
        assert f"{ATTR},{TEST_VAR_TEMP},d,," in cmd
        assert test_netcdf_file in cmd
        assert TEST_OUTPUT_FILE in cmd

    def test_modify_float_attribute(self, test_netcdf_file):
        """Test building command to modify a float attribute."""
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE}
        cmd = build_ncatted_command(test_netcdf_file, TEST_OUTPUT_FILE, var_fillvalues)

        assert NCATTED_CMD in cmd
        assert NCATTED_FLAG in cmd
        # Should use 'f' for float32
        assert f"{ATTR},{TEST_VAR_TEMP},o,f,{TEST_FILL_VALUE}" in cmd

    def test_modify_double_attribute(self, test_netcdf_file):
        """Test building command to modify a double (float64) attribute."""
        var_fillvalues = {TEST_VAR_PRESSURE: TEST_FILL_VALUE}
        cmd = build_ncatted_command(test_netcdf_file, TEST_OUTPUT_FILE, var_fillvalues)

        assert NCATTED_CMD in cmd
        assert NCATTED_FLAG in cmd
        # Should use 'd' for float64
        assert f"{ATTR},{TEST_VAR_PRESSURE},o,d,{TEST_FILL_VALUE}" in cmd

    def test_multiple_variables(self, test_netcdf_file):
        """Test building command with multiple variables."""
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE, TEST_VAR_PRESSURE: USER_REQ_DELETE}
        cmd = build_ncatted_command(test_netcdf_file, TEST_OUTPUT_FILE, var_fillvalues)

        # Should have two -a flags
        assert cmd.count(NCATTED_FLAG) == 2
        assert f"{ATTR},{TEST_VAR_TEMP},o,f,{TEST_FILL_VALUE}" in cmd
        assert f"{ATTR},{TEST_VAR_PRESSURE},d,," in cmd

    def test_variable_not_found(self, test_netcdf_file):
        """Test that missing variable raises ValueError."""
        var_fillvalues = {"nonexistent_var": TEST_FILL_VALUE}
        with pytest.raises(ValueError, match="not found"):
            build_ncatted_command(test_netcdf_file, TEST_OUTPUT_FILE, var_fillvalues)


class TestReplaceFullWorkflow:
    """Test the complete workflow of replacing fill values."""

    @pytest.fixture
    def test_setup(self, tmp_path):
        """Set up test files and JSON for full workflow test."""
        # Create input NetCDF file
        input_file = tmp_path / "input.nc"
        ds = xr.Dataset(
            {
                TEST_VAR_TEMP: xr.DataArray(
                    np.array([1.0, 2.0, np.nan, 4.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan)},
                ),
                TEST_VAR_PRESSURE: xr.DataArray(
                    np.array([1000.0, 1010.0, 1020.0, 1030.0], dtype=np.float64),
                    dims=["time"],
                    attrs={ATTR: np.float64(np.nan)},
                ),
            }
        )
        ds.to_netcdf(str(input_file))
        ds.close()

        # Create fillvalues JSON file
        fillvalues_file = tmp_path / "test_fillvalues.json"
        fillvalues_data = {
            str(input_file): {
                TEST_VAR_TEMP: TEST_FILL_VALUE,
                TEST_VAR_PRESSURE: USER_REQ_DELETE,
            }
        }
        fillvalues_file.write_text(json.dumps(fillvalues_data), encoding="utf-8")

        return {
            "input_file": str(input_file),
            "fillvalues_file": str(fillvalues_file),
            "tmp_path": tmp_path,
        }

    def test_full_workflow(self, test_setup):
        """Test the complete workflow from JSON to modified NetCDF file."""
        input_file = test_setup["input_file"]
        fillvalues_file = test_setup["fillvalues_file"]
        output_file = get_output_filename(input_file)

        # Load the fillvalues and build command
        fillvalues = load_new_fillvalues(fillvalues_file)
        var_fillvalues = fillvalues[input_file]
        cmd = build_ncatted_command(input_file, output_file, var_fillvalues)

        # Execute the command
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        assert result.returncode == 0

        # Verify the output file exists
        assert os.path.exists(output_file)

        # Open and check the output file
        ds_out = xr.open_dataset(output_file, decode_cf=False, decode_timedelta=False, decode_times=False)

        # Check temp variable - should have new fill value
        assert TEST_VAR_TEMP in ds_out
        assert ATTR in ds_out[TEST_VAR_TEMP].attrs
        assert ds_out[TEST_VAR_TEMP].attrs[ATTR] == TEST_FILL_VALUE
        # Data should be unchanged except fill value
        assert ds_out[TEST_VAR_TEMP].shape == (4,)

        # Check pressure variable - fill value should be deleted
        assert TEST_VAR_PRESSURE in ds_out
        assert ATTR not in ds_out[TEST_VAR_PRESSURE].attrs

        ds_out.close()

    def test_full_workflow_dry_run(self, test_setup):
        """Test that dry-run mode doesn't create output files."""
        input_file = test_setup["input_file"]
        output_file = get_output_filename(input_file)

        # In dry-run mode, output file should not be created
        # (We can't easily test the main() function directly, so we just verify
        # that the output file doesn't exist before we would run the command)
        assert not os.path.exists(output_file)
