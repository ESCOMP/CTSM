"""
Pytest configuration and fixtures for test_no_nans_in_inputs tests.
"""

import os
import sys

import pytest

# Add the python directory to sys.path so we can import ctsm modules
_CTSM_PYTHON = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir)
)
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
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
    Fixture to create the mock XML file with default content.

    Use this fixture in tests that need an actual XML file to exist.
    """
    xml_content = """<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="clm6_0">lnd/clm2/paramdata/test_params.nc</paramfile>
    <surfdata>lnd/clm2/surfdata/test_surf.nc</surfdata>
</namelist_defaults>
"""
    with open(mock_xml_file_path, "w", encoding="utf-8") as f:
        f.write(xml_content)

    return mock_xml_file_path
