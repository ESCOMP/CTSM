"""
Pytest configuration and fixtures for test_no_nans_in_inputs tests.
"""

import pytest

from ctsm.no_nans_in_inputs import constants
from ctsm.no_nans_in_inputs import get_replacement_fill_values
from ctsm.no_nans_in_inputs import replace_fill_values
from ctsm.no_nans_in_inputs import get_replacement_fill_values


@pytest.fixture(autouse=True, name="mock_progress_file")
def fixture_mock_progress_file(tmp_path, monkeypatch):
    """
    Auto-used fixture to mock NEW_FILLVALUES_FILE constant with a temporary path.

    This prevents tests from accidentally modifying the real file.
    Does not create the file; that happens in get_replacement_fill_values.py.
    """
    # Define the test NEW_FILLVALUES_FILE
    test_progress = tmp_path / "progress.json"

    # Monkeypatch the XML_FILE constant in both the constants module and replace_fill_values
    monkeypatch.setattr(constants, "NEW_FILLVALUES_FILE", str(test_progress))

    # Also patch it in things that import XML_FILE
    monkeypatch.setattr(get_replacement_fill_values, "NEW_FILLVALUES_FILE", str(test_progress))
    monkeypatch.setattr(replace_fill_values, "NEW_FILLVALUES_FILE", str(test_progress))

    return str(test_progress)


@pytest.fixture(autouse=True, name="mock_inputdata_prefix")
def fixture_mock_inputdata_prefix(tmp_path, monkeypatch):
    """
    Auto-used fixture to mock INPUTDATA_PREFIX constant with a temporary path.
    """
    # Monkeypatch
    monkeypatch.setattr(get_replacement_fill_values, "INPUTDATA_PREFIX", str(tmp_path))


@pytest.fixture(autouse=True, name="mock_progress_file", scope="function")
def fixture_mock_progress_file(tmp_path, monkeypatch):
    """
    Auto-used fixture to mock NEW_FILLVALUES_FILE constant with a temporary path.
    """
    # Define the test NEW_FILLVALUES_FILE file path (but don't create it yet)
    test_progress = tmp_path / "progress.json"

    # Monkeypatch the NEW_FILLVALUES_FILE constant where it's defined
    monkeypatch.setattr(constants, "NEW_FILLVALUES_FILE", str(test_progress))

    # Also where it's imported
    monkeypatch.setattr(get_replacement_fill_values, "NEW_FILLVALUES_FILE", str(test_progress))
    monkeypatch.setattr(replace_fill_values, "NEW_FILLVALUES_FILE", str(test_progress))

    return str(test_progress)


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

    # Also patch it in things that import XML_FILE
    monkeypatch.setattr(get_replacement_fill_values, "XML_FILE", str(test_xml))
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
