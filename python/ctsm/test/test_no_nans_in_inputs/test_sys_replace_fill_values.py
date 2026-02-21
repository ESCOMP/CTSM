#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

import json
import os
import subprocess
import xml.etree.ElementTree as ET

from netCDF4 import Dataset  # pylint: disable=no-name-in-module
import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR, USER_REQ_DELETE
from ctsm.no_nans_in_inputs.replace_fill_values import (
    build_ncatted_command,
    execute_command,
    get_output_filename,
    load_new_fillvalues,
    main,
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

    def test_same_input_output(self, test_netcdf_file):
        """Test that same input and output files raises ValueError."""
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE}
        with pytest.raises(ValueError, match="Input and output files are the same"):
            build_ncatted_command(test_netcdf_file, test_netcdf_file, var_fillvalues)

    def test_output_already_exists(self, test_netcdf_file, tmp_path):
        """Test that existing output file is detected."""
        # Create an output file that already exists
        output_file = tmp_path / "output.nc"
        output_file.write_text("dummy content")

        # The function should still build the command (skip logic is in main())
        # but we can verify the output file exists
        assert os.path.exists(str(output_file))

        # Building command should still work - skip logic is in main()
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE}
        cmd = build_ncatted_command(test_netcdf_file, str(output_file), var_fillvalues)
        assert NCATTED_CMD in cmd


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
        ds_out = xr.open_dataset(
            output_file, decode_cf=False, decode_timedelta=False, decode_times=False
        )

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

    def test_full_workflow_dry_run(self, test_setup, monkeypatch):
        """Test that dry-run mode doesn't create output files."""
        input_file = test_setup["input_file"]
        fillvalues_file = test_setup["fillvalues_file"]
        output_file = get_output_filename(input_file)

        # Mock sys.argv to simulate running with --dry-run
        monkeypatch.setattr(
            "sys.argv",
            ["replace_fill_values.py", "--fillvalues-file", fillvalues_file, "--dry-run"],
        )

        # Run main in dry-run mode
        result = main()
        assert result == 0

        # Output file should not have been created
        assert not os.path.exists(output_file)

    def test_skip_existing_without_overwrite(self, test_setup, capsys, monkeypatch):
        """Test that existing output files are skipped without --overwrite flag."""
        input_file = test_setup["input_file"]
        fillvalues_file = test_setup["fillvalues_file"]
        output_file = get_output_filename(input_file)

        # Create a dummy output file to simulate it already existing
        with open(output_file, "w", encoding="utf-8") as f:
            f.write("dummy existing file")

        original_mtime = os.path.getmtime(output_file)

        # Mock sys.argv to simulate running without --overwrite
        monkeypatch.setattr(
            "sys.argv", ["replace_fill_values.py", "--fillvalues-file", fillvalues_file]
        )

        # Run main (should skip the file)
        result = main()
        assert result == 0

        # Output file should be unchanged
        assert os.path.getmtime(output_file) == original_mtime

        # Check that skip message was printed
        captured = capsys.readouterr()
        assert "Skipping (output exists)" in captured.out
        assert "Use --overwrite" in captured.out

    def test_process_with_overwrite(self, test_setup, monkeypatch, create_mock_xml_file):
        """Test that existing output files are overwritten with --overwrite flag."""
        input_file = test_setup["input_file"]
        fillvalues_file = test_setup["fillvalues_file"]
        output_file = get_output_filename(input_file)

        # Create mock XML file with our test input file
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <test_paramfile>{input_file}</test_paramfile>
</namelist_defaults>
"""
        )

        # Create a dummy output file
        with open(output_file, "w", encoding="utf-8") as f:
            f.write("dummy existing file")

        assert os.path.exists(output_file)

        # Mock sys.argv to simulate running with --overwrite (XML file is auto-mocked)
        monkeypatch.setattr(
            "sys.argv",
            ["replace_fill_values.py", "--fillvalues-file", fillvalues_file, "--overwrite"],
        )

        # Run main (should replace the file)
        result = main()
        assert result == 0

        # Output file should now be a valid NetCDF file
        ds_out = xr.open_dataset(
            output_file, decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert TEST_VAR_TEMP in ds_out
        assert ATTR in ds_out[TEST_VAR_TEMP].attrs
        assert ds_out[TEST_VAR_TEMP].attrs[ATTR] == TEST_FILL_VALUE
        ds_out.close()

        # Verify XML was updated
        tree = ET.parse(xml_path)
        root = tree.getroot()
        test_paramfile = root.find("test_paramfile")
        assert test_paramfile is not None
        assert test_paramfile.text == output_file

    def test_skip_symlink_output(self, test_setup, capsys, monkeypatch):
        """Test that symlinked output files are never overwritten, even with --overwrite."""
        input_file = test_setup["input_file"]
        fillvalues_file = test_setup["fillvalues_file"]
        output_file = get_output_filename(input_file)
        tmp_path = test_setup["tmp_path"]

        # Create a target file and symlink to it
        target_file = tmp_path / "target.nc"
        target_file.write_text("target file content")
        os.symlink(str(target_file), output_file)

        assert os.path.islink(output_file)
        original_content = target_file.read_text()

        # Mock sys.argv with --overwrite (should still skip symlinks)
        monkeypatch.setattr(
            "sys.argv",
            ["replace_fill_values.py", "--fillvalues-file", fillvalues_file, "--overwrite"],
        )

        # Run main (should skip the symlink)
        result = main()
        assert result == 0

        # Target file should be unchanged
        assert target_file.read_text() == original_content

        # Check that warning was printed
        captured = capsys.readouterr()
        assert "WARNING: Output file is a symlink - SKIPPING" in captured.out
        assert "Symlinks will never be overwritten" in captured.out


class TestExecuteCommand:
    """Test the execute_command function."""

    @pytest.mark.parametrize(
        "netcdf_format",
        [
            "NETCDF4",
            "NETCDF4_CLASSIC",
            "NETCDF3_CLASSIC",
        ],
    )
    def test_execute_preserves_format(self, tmp_path, netcdf_format):
        """Test execute_command preserves NetCDF format from input to output."""
        # Create input file with specified format
        input_file = str(tmp_path / f"test_input_{netcdf_format}.nc")
        ds = xr.Dataset(
            {
                TEST_VAR_TEMP: xr.DataArray(
                    np.array([1.0, 2.0, 3.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan)},
                ),
            }
        )
        ds.to_netcdf(input_file, format=netcdf_format)
        ds.close()

        output_file = str(tmp_path / f"test_output_{netcdf_format}.nc")
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE}
        cmd = build_ncatted_command(input_file, output_file, var_fillvalues)

        # Execute the command without XML update
        files_processed = execute_command(None, input_file, output_file, cmd)

        # Should have processed 1 file
        assert files_processed == 1

        # Output file should exist
        assert os.path.exists(output_file)

        # Verify the output file is valid NetCDF with correct fill value
        ds_out = xr.open_dataset(
            output_file, decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert TEST_VAR_TEMP in ds_out
        assert ATTR in ds_out[TEST_VAR_TEMP].attrs
        assert ds_out[TEST_VAR_TEMP].attrs[ATTR] == TEST_FILL_VALUE

        # Verify NetCDF format matches input
        ds_in = xr.open_dataset(
            input_file, decode_cf=False, decode_timedelta=False, decode_times=False
        )

        # Check the actual file format using netCDF4 library
        nc_in = Dataset(input_file, "r")
        nc_out = Dataset(output_file, "r")

        # Get the actual format from the files
        input_format = nc_in.data_model
        output_format = nc_out.data_model

        nc_in.close()
        nc_out.close()
        ds_out.close()
        ds_in.close()

        # Verify formats match
        assert output_format == input_format
        assert output_format == netcdf_format

    def test_execute_without_xml_update(self, tmp_path, create_mock_xml_file):
        """Test that XML file is not touched when xml_file is None."""
        # Create input file
        input_file = str(tmp_path / "test_input.nc")
        ds = xr.Dataset(
            {
                TEST_VAR_TEMP: xr.DataArray(
                    np.array([1.0, 2.0, 3.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan)},
                ),
            }
        )
        ds.to_netcdf(input_file)
        ds.close()

        output_file = str(tmp_path / "test_output.nc")
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE}
        cmd = build_ncatted_command(input_file, output_file, var_fillvalues)

        # Create mock XML and read original content
        xml_path = create_mock_xml_file()
        with open(xml_path, "r", encoding="utf-8") as f:
            original_xml = f.read()

        # Execute without xml_file parameter
        files_processed = execute_command(None, input_file, output_file, cmd)

        assert files_processed == 1

        # XML file should be unchanged
        with open(xml_path, "r", encoding="utf-8") as f:
            current_xml = f.read()
        assert current_xml == original_xml
