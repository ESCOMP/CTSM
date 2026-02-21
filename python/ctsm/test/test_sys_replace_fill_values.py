#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

import json
import os
import tempfile

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR, USER_REQ_DELETE
from ctsm.no_nans_in_inputs.replace_fill_values import (
    build_ncatted_command,
    load_new_fillvalues,
)


class TestLoadNewFillvalues:
    """Test the load_new_fillvalues function."""

    def test_load_valid_json(self):
        """Test loading a valid JSON file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            test_data = {
                "/path/to/file1.nc": {"var1": -999.0, "var2": "delete"},
                "/path/to/file2.nc": {"var3": -999},
            }
            json.dump(test_data, f)
            temp_file = f.name

        try:
            result = load_new_fillvalues(temp_file)
            assert result == test_data
        finally:
            os.unlink(temp_file)

    def test_file_not_found(self):
        """Test that missing file causes SystemExit."""
        with pytest.raises(SystemExit):
            load_new_fillvalues("/nonexistent/file.json")


class TestBuildNcattedCommand:
    """Test the build_ncatted_command function."""

    @pytest.fixture
    def test_netcdf_file(self):
        """Create a temporary NetCDF file for testing."""
        temp_dir = tempfile.mkdtemp()
        test_file = os.path.join(temp_dir, "test.nc")

        # Create a simple NetCDF file with float variables that have NaN fill values
        # (NetCDF doesn't allow NaN for integer types, and our scripts only work on
        # variables that already have NaN fill values)
        ds = xr.Dataset(
            {
                "temp": xr.DataArray(
                    np.array([1.0, 2.0, 3.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan)},
                ),
                "pressure": xr.DataArray(
                    np.array([1000.0, 1010.0, 1020.0], dtype=np.float64),
                    dims=["time"],
                    attrs={ATTR: np.float64(np.nan)},
                ),
            }
        )
        ds.to_netcdf(test_file)
        ds.close()

        yield test_file

        # Cleanup
        if os.path.exists(test_file):
            os.unlink(test_file)
        os.rmdir(temp_dir)

    def test_delete_attribute(self, test_netcdf_file):
        """Test building command to delete an attribute."""
        var_fillvalues = {"temp": USER_REQ_DELETE}
        cmd = build_ncatted_command(test_netcdf_file, "output.nc", var_fillvalues)

        assert "ncatted" in cmd
        assert "-a" in cmd
        assert f"{ATTR},temp,d,," in cmd
        assert test_netcdf_file in cmd
        assert "output.nc" in cmd

    def test_modify_float_attribute(self, test_netcdf_file):
        """Test building command to modify a float attribute."""
        var_fillvalues = {"temp": -999.0}
        cmd = build_ncatted_command(test_netcdf_file, "output.nc", var_fillvalues)

        assert "ncatted" in cmd
        assert "-a" in cmd
        # Should use 'f' for float32
        assert f"{ATTR},temp,o,f,-999.0" in cmd

    def test_modify_double_attribute(self, test_netcdf_file):
        """Test building command to modify a double (float64) attribute."""
        var_fillvalues = {"pressure": -999.0}
        cmd = build_ncatted_command(test_netcdf_file, "output.nc", var_fillvalues)

        assert "ncatted" in cmd
        assert "-a" in cmd
        # Should use 'd' for float64
        assert f"{ATTR},pressure,o,d,-999.0" in cmd

    def test_multiple_variables(self, test_netcdf_file):
        """Test building command with multiple variables."""
        var_fillvalues = {"temp": -999.0, "pressure": USER_REQ_DELETE}
        cmd = build_ncatted_command(test_netcdf_file, "output.nc", var_fillvalues)

        # Should have two -a flags
        assert cmd.count("-a") == 2
        assert f"{ATTR},temp,o,f,-999.0" in cmd
        assert f"{ATTR},pressure,d,," in cmd

    def test_variable_not_found(self, test_netcdf_file):
        """Test that missing variable raises ValueError."""
        var_fillvalues = {"nonexistent_var": -999.0}
        with pytest.raises(ValueError, match="not found"):
            build_ncatted_command(test_netcdf_file, "output.nc", var_fillvalues)
