#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

import os
import subprocess
import xml.etree.ElementTree as ET
from unittest.mock import patch

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR, USER_REQ_DELETE, OPEN_DS_KWARGS
from ctsm.no_nans_in_inputs.json_io import NoNanFillValueProgress
from ctsm.no_nans_in_inputs.replace_fill_values import (
    get_output_filename,
    main,
)
from ctsm.no_nans_in_inputs.netcdf_utils import build_ncatted_command


# Test constants
TEST_VAR_TEMP = "temp"
TEST_VAR_PRESSURE = "pressure"
TEST_OUTPUT_FILE = "output.nc"
TEST_FILL_VALUE = -123.4
NCATTED_CMD = "ncatted"
NCATTED_FLAG = "-a"


class TestReplaceFullWorkflow:
    """Test the complete workflow of replacing fill values."""

    @pytest.fixture
    def test_setup(self, tmp_path, create_mock_xml_file):
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

        # Write the XML file
        xml_content = f"""<?xml version="1.0"?>
    <namelist_defaults>
        <paramfile>{input_file}</paramfile>
    </namelist_defaults>
    """
        xml_file = create_mock_xml_file(xml_content)

        # Create fillvalues JSON file
        fillvalues_file = tmp_path / "test_fillvalues.json"
        fillvalues_data = NoNanFillValueProgress(progress_file=str(fillvalues_file))
        input_file_dict = fillvalues_data[input_file]
        input_file_dict["found_in_files"] = {xml_file: {str(input_file)}}
        input_file_dict["new_fill_values"] = {
            TEST_VAR_TEMP: TEST_FILL_VALUE,
            TEST_VAR_PRESSURE: USER_REQ_DELETE,
        }
        fillvalues_data.save()

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
        fillvalues = NoNanFillValueProgress(
            progress_file=str(fillvalues_file), load_without_asking=True
        )
        var_fillvalues = fillvalues[input_file]["new_fill_values"]
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
        fillvalues_file = str(test_setup["fillvalues_file"])
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
        with patch("builtins.input", return_value="y"):  # continue after replacing
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
