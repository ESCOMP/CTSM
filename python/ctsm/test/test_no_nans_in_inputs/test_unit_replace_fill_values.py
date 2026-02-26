#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

import numpy as np
import pytest

from ctsm.no_nans_in_inputs.replace_fill_values import (
    get_ncatted_type_code,
    get_output_filename,
)


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


class TestGetNcattedTypeCode:
    """Test the get_ncatted_type_code function."""

    def test_float64(self):
        """Test float64 dtype."""
        assert get_ncatted_type_code(np.dtype("float64")) == "d"

    def test_float32(self):
        """Test float32 dtype."""
        assert get_ncatted_type_code(np.dtype("float32")) == "f"

    def test_int64(self):
        """Test int64 dtype raises error."""
        with pytest.raises(ValueError, match="Integer dtype detected"):
            get_ncatted_type_code(np.dtype("int64"))

    def test_int32(self):
        """Test int32 dtype raises error."""
        with pytest.raises(ValueError, match="Integer dtype detected"):
            get_ncatted_type_code(np.dtype("int32"))

    def test_int16(self):
        """Test int16 dtype raises error."""
        with pytest.raises(ValueError, match="Integer dtype detected"):
            get_ncatted_type_code(np.dtype("int16"))

    def test_int8(self):
        """Test int8 dtype raises error."""
        with pytest.raises(ValueError, match="Integer dtype detected"):
            get_ncatted_type_code(np.dtype("int8"))

    def test_unknown_dtype(self):
        """Test that unknown dtype raises ValueError."""
        with pytest.raises(ValueError, match="Unknown dtype"):
            get_ncatted_type_code(np.dtype("complex128"))
