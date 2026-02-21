"""
Pytest configuration and fixtures for test_no_nans_in_inputs tests.
"""

import pytest

from ctsm.no_nans_in_inputs import constants
from ctsm.no_nans_in_inputs import replace_fill_values


@pytest.fixture(autouse=True, name="mock_xml_file_path")
def fixture_mock_xml_file_path(tmp_path, monkeypatch):
    """
    Auto-used fixture to mock XML_FILE constant with a temporary path.

    This prevents tests from accidentally modifying the real XML file.
    Does not create the file - use create_mock_xml_file fixture for that.
    """
    # Define the test XML file path (but don't create it yet)
    test_xml = tmp_path / "test_namelist_defaults.xml"

    # Monkeypatch the XML_FILE constant in both the constants module and replace_fill_values
    monkeypatch.setattr(constants, "XML_FILE", str(test_xml))

    # Also patch it in replace_fill_values since it imports XML_FILE
    monkeypatch.setattr(replace_fill_values, "XML_FILE", str(test_xml))

    return str(test_xml)


@pytest.fixture(name="create_mock_xml_file")
def fixture_create_mock_xml_file(mock_xml_file_path):
    """
    Factory fixture to create the mock XML file with given or default content.

    Use this fixture in tests that need an actual XML file to exist.
    Call with no arguments for default content, or pass custom XML content.

    Returns:
        A function that creates the XML file and returns its path.
    """
    default_content = """<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="clm6_0">lnd/clm2/paramdata/test_params.nc</paramfile>
    <surfdata>lnd/clm2/surfdata/test_surf.nc</surfdata>
</namelist_defaults>
"""

    def _create(content: str = default_content) -> str:
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(content)
        return mock_xml_file_path

    return _create
