#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

from ctsm.no_nans_in_inputs.replace_fill_values import get_output_filename


class TestGetOutputFilename:
    """Test the get_output_filename function."""

    def test_simple_nc_file(self):
        """Test with simple .nc file."""
        result = get_output_filename("/path/to/file.nc")
        assert result == "/path/to/file.no_nan_fill.nc"

    def test_double_extension(self):
        """Test with double extension like .tar.gz."""
        result = get_output_filename("/path/to/file.tar.gz")
        assert result == "/path/to/file.tar.no_nan_fill.gz"

    def test_no_extension(self):
        """Test with file without extension."""
        result = get_output_filename("/path/to/file")
        assert result == "/path/to/file.no_nan_fill"

    def test_no_directory(self):
        """Test with filename only (no directory)."""
        result = get_output_filename("file.nc")
        assert result == "file.no_nan_fill.nc"
