#!/usr/bin/env python3
"""
System tests for replace_fill_values.py script.

Tests the functionality of replacing NaN fill values in NetCDF files.
"""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ctsm.no_nans_in_inputs.replace_fill_values import (
    get_ncatted_type_code,
    get_output_filename,
    update_xml_file,
)


# Test constants
TEST_PARAM_CLM60 = "lnd/clm2/paramdata/ctsm60_params.c260108.nc"
TEST_PARAM_CLM50 = "lnd/clm2/paramdata/clm50_params.c250311.nc"
TEST_PARAM_CLM45 = "lnd/clm2/paramdata/clm45_params.c250311.nc"
TEST_PHYS_CLM60 = "clm6_0"
TEST_PHYS_CLM50 = "clm5_0"
TEST_PHYS_CLM45 = "clm4_5"


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


class TestUpdateXmlFile:
    """Test the update_xml_file function."""

    def test_update_xml_path(self, mock_xml_file_path):
        """Test updating a file path in XML."""
        # Create XML content in the mocked path
        xml_content = """<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>lnd/clm2/paramdata/test_params.nc</paramfile>
    <surfdata>lnd/clm2/surfdata/test_surf.nc</surfdata>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        # Update one of the paths
        old_path = "lnd/clm2/paramdata/test_params.nc"
        new_path = "lnd/clm2/paramdata/test_params.no_nan_fill.nc"

        update_xml_file(mock_xml_file_path, old_path, new_path)

        # Read and verify the updated XML
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()
        paramfile = root.find("paramfile")
        assert paramfile is not None
        assert paramfile.text == new_path

        # Other paths should be unchanged
        surfdata = root.find("surfdata")
        assert surfdata.text == "lnd/clm2/surfdata/test_surf.nc"

    def test_update_xml_path_not_found(self, mock_xml_file_path):
        """Test that updating non-existent path raises ValueError."""
        xml_content = """<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>lnd/clm2/paramdata/test_params.nc</paramfile>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        with pytest.raises(ValueError, match="not found"):
            update_xml_file(mock_xml_file_path, "nonexistent/path.nc", "new/path.nc")

    def test_update_xml_with_multiple_same_tag(self, mock_xml_file_path):
        """Test updating path when multiple elements have the same tag name."""
        # Simulate the real XML structure with multiple paramfile elements
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="{TEST_PHYS_CLM60}">{TEST_PARAM_CLM60}</paramfile>
    <paramfile phys="{TEST_PHYS_CLM50}">{TEST_PARAM_CLM50}</paramfile>
    <paramfile phys="{TEST_PHYS_CLM45}">{TEST_PARAM_CLM45}</paramfile>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        # Update one specific path
        old_path = TEST_PARAM_CLM60
        new_path = TEST_PARAM_CLM60.replace(".nc", ".no_nan_fill.nc")

        update_xml_file(mock_xml_file_path, old_path, new_path)

        # Read and verify the updated XML
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()

        # Find all paramfile elements
        paramfiles = root.findall("paramfile")
        assert len(paramfiles) == 3

        # Check that only the clm6_0 one was updated
        clm60_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM60][0]
        assert clm60_param.text == new_path

        # Check that others are unchanged
        clm50_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM50][0]
        assert clm50_param.text == TEST_PARAM_CLM50

        clm45_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM45][0]
        assert clm45_param.text == TEST_PARAM_CLM45

    def test_update_xml_replaces_within_element_text(self, mock_xml_file_path):
        """Test that all occurrences within a single element's text are replaced."""
        # Create XML with a path appearing multiple times in one element's text
        test_path = "lnd/clm2/test/file.nc"
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <multi_path>{test_path} {test_path}</multi_path>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        new_path = test_path.replace(".nc", ".no_nan_fill.nc")
        update_xml_file(mock_xml_file_path, test_path, new_path)

        # Verify both occurrences were replaced
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()
        multi_path = root.find("multi_path")
        assert multi_path is not None
        assert multi_path.text == f"{new_path} {new_path}"
        # Verify old path is completely gone
        assert test_path not in multi_path.text

    def test_update_xml_replaces_across_different_tags(self, mock_xml_file_path):
        """Test that same path in different element types are all replaced."""
        # Create XML with same path in multiple different element types
        test_path = "lnd/clm2/test/shared_file.nc"
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{test_path}</paramfile>
    <surfdata>{test_path}</surfdata>
    <initdata>{test_path}</initdata>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        new_path = test_path.replace(".nc", ".no_nan_fill.nc")
        update_xml_file(mock_xml_file_path, test_path, new_path)

        # Verify all three elements were updated
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()

        paramfile = root.find("paramfile")
        assert paramfile is not None
        assert paramfile.text == new_path

        surfdata = root.find("surfdata")
        assert surfdata is not None
        assert surfdata.text == new_path

        initdata = root.find("initdata")
        assert initdata is not None
        assert initdata.text == new_path

    def test_update_xml_replaces_across_same_tag_different_attrs(self, mock_xml_file_path):
        """Test that same path in elements with same tag but different attributes are all replaced."""
        # Simulate scenario where two paramfile elements with different attributes point to same file
        # (like if lines 625 and 626 in the real XML both had the same path)
        test_path = "lnd/clm2/paramdata/shared_params.nc"
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="{TEST_PHYS_CLM60}">{test_path}</paramfile>
    <paramfile phys="{TEST_PHYS_CLM50}">{test_path}</paramfile>
    <surfdata>lnd/clm2/surfdata/different_file.nc</surfdata>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        new_path = test_path.replace(".nc", ".no_nan_fill.nc")
        update_xml_file(mock_xml_file_path, test_path, new_path)

        # Verify both paramfile elements were updated
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()

        paramfiles = root.findall("paramfile")
        assert len(paramfiles) == 2

        # Both should have the new path
        clm60_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM60][0]
        assert clm60_param.text == new_path

        clm50_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM50][0]
        assert clm50_param.text == new_path

        # Surfdata should be unchanged
        surfdata = root.find("surfdata")
        assert surfdata.text == "lnd/clm2/surfdata/different_file.nc"
