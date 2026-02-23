#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

import json
import os
import subprocess
import xml.etree.ElementTree as ET
from unittest import mock

from netCDF4 import Dataset  # pylint: disable=no-name-in-module
import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR, USER_REQ_DELETE, OPEN_DS_KWARGS
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
TEST_FILL_VALUE = -123.4
NCATTED_CMD = "ncatted"
NCATTED_FLAG = "-a"


class TestLoadNewFillvalues:
    """Test the load_new_fillvalues function."""

    def test_load_valid_json(self, tmp_path):
        """Test loading a valid JSON file."""
        test_file = tmp_path / "test_fillvalues.json"
        test_data = {
            "/path/to/file1.nc": {"var1": -999.0, "var2": USER_REQ_DELETE},
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
            output_file,
            **OPEN_DS_KWARGS,
        )

        # Check temp variable - should have new fill value
        assert TEST_VAR_TEMP in ds_out
        assert ATTR in ds_out[TEST_VAR_TEMP].encoding
        assert ds_out[TEST_VAR_TEMP].encoding[ATTR] == TEST_FILL_VALUE
        # Data should be unchanged except fill value
        assert ds_out[TEST_VAR_TEMP].shape == (4,)

        # Check pressure variable - fill value should be deleted
        assert TEST_VAR_PRESSURE in ds_out
        assert ATTR not in ds_out[TEST_VAR_PRESSURE].encoding

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
            output_file,
            **OPEN_DS_KWARGS,
        )
        assert TEST_VAR_TEMP in ds_out
        assert ATTR in ds_out[TEST_VAR_TEMP].encoding
        assert ds_out[TEST_VAR_TEMP].encoding[ATTR] == TEST_FILL_VALUE
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

    @pytest.fixture(name="create_test_nc")
    def fixture_create_test_nc(self, tmp_path):
        """
        Factory fixture to create a test netCDF file with given or default parameters

        Returns:
            A function that creates the test netCDF file and returns its path.
        """

        def _create(
            *, netcdf_format: str = "NETCDF4_CLASSIC", first_value: np.float32 = 1.0
        ) -> str:
            test_nc = str(tmp_path / "test_input.nc")
            nan_fill = np.float32(np.nan)
            ds = xr.Dataset(
                {
                    TEST_VAR_TEMP: xr.DataArray(
                        np.array([first_value, 2.0, 3.0], dtype=np.float32),
                        dims=["time"],
                    ),
                }
            )
            ds.to_netcdf(test_nc, format=netcdf_format, encoding={TEST_VAR_TEMP: {ATTR: nan_fill}})
            ds.close()

            # Check fill value
            ds = xr.open_dataset(test_nc, mask_and_scale=True, **OPEN_DS_KWARGS)
            assert ATTR in ds[TEST_VAR_TEMP].encoding
            assert np.isnan(ds[TEST_VAR_TEMP].encoding[ATTR])
            ds.close()

            return test_nc

        return _create

    @pytest.fixture(name="mock_update_xml_file", autouse=True)
    def fixture_mock_update_xml_file(self):
        """Every test in this class will have update_xml_file() mocked; that's tested elsewhere"""
        with mock.patch("ctsm.no_nans_in_inputs.replace_fill_values.update_xml_file") as _fixture:
            yield _fixture

    # TODO: This doesn't actually touch the file system, so it's really more of a unit test
    @mock.patch("subprocess.run")
    @pytest.mark.parametrize("xml_file, exp_n_calls", [("dummy.xml", 1), (None, 0)])
    def test_update_xml_file(
        self, _mock_subprocess_run, xml_file, exp_n_calls, mock_update_xml_file
    ):
        """Test that execute_command does or doesn't call update_xml_file(), as appropriate"""

        # Build and execute the command
        execute_command(xml_file, "dummy_input.nc", "dummy_output.nc", "dummy command --option")

        # update_xml_file() should have been called exp_n_calls times
        assert mock_update_xml_file.call_count == exp_n_calls

    @pytest.mark.parametrize(
        "netcdf_format",
        [
            "NETCDF4",
            "NETCDF4_CLASSIC",
            "NETCDF3_CLASSIC",
        ],
    )
    def test_execute_preserves_format(self, tmp_path, create_test_nc, netcdf_format):
        """Test execute_command preserves NetCDF format from input to output."""
        # Create input file with specified format
        input_file = create_test_nc(netcdf_format=netcdf_format)

        # Build and execute the command without XML update (because we're not testing that here)
        output_file = str(tmp_path / "test_output.nc")
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE}
        cmd = build_ncatted_command(input_file, output_file, var_fillvalues)
        files_processed = execute_command(None, input_file, output_file, cmd)

        # Should have processed 1 file
        assert files_processed == 1

        # Output file should exist
        assert os.path.exists(output_file)

        # Verify the output file is valid NetCDF with correct fill value
        ds_out = xr.open_dataset(
            output_file,
            **OPEN_DS_KWARGS,
        )
        assert TEST_VAR_TEMP in ds_out
        assert ATTR in ds_out[TEST_VAR_TEMP].encoding
        assert ds_out[TEST_VAR_TEMP].encoding[ATTR] == TEST_FILL_VALUE
        ds_out.close()

        # Verify formats match
        with Dataset(input_file, "r") as nc_in:
            input_format = nc_in.data_model
        with Dataset(output_file, "r") as nc_out:
            output_format = nc_out.data_model
        assert output_format == input_format
        assert output_format == netcdf_format

    @pytest.mark.parametrize(
        "first_value, expect_filled",
        [
            (np.nan, True),
            (TEST_FILL_VALUE, True),
            (150, False),
            (7.24, False),
        ],
    )
    def test_execute_different_data(self, tmp_path, create_test_nc, first_value, expect_filled):
        """Test execute_command with different data values."""
        # Create input file with specified value first
        input_file = create_test_nc(first_value=first_value)

        # Build and execute the command without XML update (because we're not testing that here)
        output_file = str(tmp_path / "test_output.nc")
        var_fillvalues = {TEST_VAR_TEMP: TEST_FILL_VALUE}
        cmd = build_ncatted_command(input_file, output_file, var_fillvalues)
        execute_command(None, input_file, output_file, cmd)

        # Verify the output file is valid NetCDF with correct fill value. mask_and_scale True means
        # that the variable's DataArray's _FillValue attribute will be populated and any filled
        # values will be NaN.
        assert os.path.exists(output_file)
        ds_out = xr.open_dataset(output_file, mask_and_scale=True, **OPEN_DS_KWARGS)
        assert TEST_VAR_TEMP in ds_out
        assert ATTR in ds_out[TEST_VAR_TEMP].encoding
        assert ds_out[TEST_VAR_TEMP].encoding[ATTR] == TEST_FILL_VALUE

        # If we used a value we expect to be filled, then...
        if expect_filled:
            # It should be NaN after reading with mask_and_scale True
            assert np.isnan(ds_out[TEST_VAR_TEMP].values[0]), str(ds_out[TEST_VAR_TEMP].values)
            # It should be the fill value after reading with mask_and_scale False
            ds_out_no_ms = xr.open_dataset(output_file, mask_and_scale=False, **OPEN_DS_KWARGS)
            assert ds_out_no_ms[TEST_VAR_TEMP].values[0] == TEST_FILL_VALUE
            ds_out_no_ms.close()
        # Otherwise, we expect it to be unchanged
        else:
            assert ds_out[TEST_VAR_TEMP].values[0] == first_value

        ds_out.close()
