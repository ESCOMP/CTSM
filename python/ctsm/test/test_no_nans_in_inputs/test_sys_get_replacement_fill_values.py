#!/usr/bin/env python3
"""
System tests for get_replacement_fill_values.py script.

Tests functions that require file I/O.
"""

import pytest

from ctsm.no_nans_in_inputs.get_replacement_fill_values import extract_file_paths_from_xml


# Test constants
TEST_PATH_PARAM = "lnd/clm2/paramdata/test_params.nc"
TEST_PATH_SURF = "lnd/clm2/surfdata/test_surf.nc"
TEST_PATH_INIT = "lnd/clm2/initdata/test_init.nc"
TEST_PATH_OTHER = "share/meshes/test_mesh.nc"


class TestExtractFilePathsFromXml:
    """Test the extract_file_paths_from_xml function."""

    def test_extracts_paths_from_default_content(self, create_mock_xml_file):
        """Test extracting paths from the default mock XML content."""
        result = extract_file_paths_from_xml(create_mock_xml_file)
        assert TEST_PATH_PARAM in result
        assert TEST_PATH_SURF in result

    def test_extracts_multiple_paths(self, mock_xml_file_path):
        """Test extracting multiple file paths from XML."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{TEST_PATH_PARAM}</paramfile>
    <surfdata>{TEST_PATH_SURF}</surfdata>
    <initdata>{TEST_PATH_INIT}</initdata>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == {TEST_PATH_PARAM, TEST_PATH_SURF, TEST_PATH_INIT}

    def test_ignores_non_lnd_clm2_paths(self, mock_xml_file_path):
        """Test that paths not containing lnd/clm2/ are ignored."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{TEST_PATH_PARAM}</paramfile>
    <meshfile>{TEST_PATH_OTHER}</meshfile>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == {TEST_PATH_PARAM}
        assert TEST_PATH_OTHER not in result

    def test_ignores_non_path_text(self, mock_xml_file_path):
        """Test that non-path text content is ignored."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{TEST_PATH_PARAM}</paramfile>
    <some_setting>42</some_setting>
    <another_setting>.true.</another_setting>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == {TEST_PATH_PARAM}

    def test_deduplicates_paths(self, mock_xml_file_path):
        """Test that duplicate paths are returned only once."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="clm6_0">{TEST_PATH_PARAM}</paramfile>
    <paramfile phys="clm5_0">{TEST_PATH_PARAM}</paramfile>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == {TEST_PATH_PARAM}

    def test_handles_whitespace_around_path(self, mock_xml_file_path):
        """Test that leading/trailing whitespace around paths is handled."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>  {TEST_PATH_PARAM}  </paramfile>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == {TEST_PATH_PARAM}

    def test_handles_multiline_path(self, mock_xml_file_path):
        """Test extracting a path that spans multiple lines (with leading newline)."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>
{TEST_PATH_PARAM}
</paramfile>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == {TEST_PATH_PARAM}

    def test_empty_xml(self, mock_xml_file_path):
        """Test with XML that has no file paths."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write("""<?xml version="1.0"?>
<namelist_defaults>
    <some_setting>42</some_setting>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == set()

    def test_file_not_found_exits(self):
        """Test that a missing XML file causes SystemExit."""
        with pytest.raises(SystemExit):
            extract_file_paths_from_xml("/nonexistent/file.xml")

    def test_handles_elements_with_attributes(self, mock_xml_file_path):
        """Test extracting paths from elements that have attributes."""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="clm6_0">{TEST_PATH_PARAM}</paramfile>
    <paramfile phys="clm5_0">{TEST_PATH_SURF}</paramfile>
</namelist_defaults>
""")
        result = extract_file_paths_from_xml(mock_xml_file_path)
        assert result == {TEST_PATH_PARAM, TEST_PATH_SURF}
