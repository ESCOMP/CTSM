"""
Pytest configuration and fixtures for test_no_nans_in_inputs tests.
"""

import os
from pathlib import Path
from dataclasses import dataclass

import pytest

from ctsm.no_nans_in_inputs import constants
from ctsm.no_nans_in_inputs import get_replacement_fill_values
from ctsm.no_nans_in_inputs import json_io
from ctsm.no_nans_in_inputs import namelist_utils
from ctsm.no_nans_in_inputs import replace_fill_values
from ctsm.no_nans_in_inputs import shared
from ctsm.no_nans_in_inputs import user_inputs


@pytest.fixture(autouse=True, name="mock_inputdata_prefix")
def fixture_mock_inputdata_prefix(tmp_path, monkeypatch):
    """
    Auto-used fixture to mock INPUTDATA_PREFIX constant with a temporary path.
    """
    # Monkeypatch
    monkeypatch.setattr(constants, "INPUTDATA_PREFIX", str(tmp_path))
    monkeypatch.setattr(namelist_utils, "INPUTDATA_PREFIX", str(tmp_path))
    monkeypatch.setattr(shared, "INPUTDATA_PREFIX", str(tmp_path))
    monkeypatch.setattr(replace_fill_values, "INPUTDATA_PREFIX", str(tmp_path))


@pytest.fixture(autouse=True, name="mock_dir_to_search_for_usernl_files")
def fixture_mock_dir_with_usernl_files(tmp_path, monkeypatch):
    """
    Auto-used fixture to mock DIR_TO_SEARCH_FOR_USER_NL_FILES constant with a temporary path.
    """
    # Monkeypatch
    monkeypatch.setattr(
        get_replacement_fill_values, "DIR_TO_SEARCH_FOR_USER_NL_FILES", str(tmp_path)
    )


@pytest.fixture(autouse=True, name="mock_cesm_top")
def fixture_mock_cesm_top(tmp_path, monkeypatch):
    """
    Auto-used fixture to mock CESM_TOP constant with a temporary path.
    """
    # Monkeypatch
    monkeypatch.setattr(constants, "CESM_TOP", str(tmp_path))
    monkeypatch.setattr(replace_fill_values, "CESM_TOP", str(tmp_path))


@pytest.fixture(autouse=True, name="mock_progress_file")
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
    monkeypatch.setattr(user_inputs, "NEW_FILLVALUES_FILE", str(test_progress))

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


@pytest.fixture(name="create_mock_user_nl_file")
def fixture_create_mock_user_nl_file(
    tmp_path, mock_inputdata_prefix
):  # pylint: disable=unused-argument
    """
    Factory fixture to create the mock user_nl_ file with given content.

    Use this fixture in tests that need an actual user_nl_ file to exist.

    Returns:
        A function that creates the file and returns its path.
    """

    @dataclass
    class NcPaths:
        """Class for holding our test paths"""

        rel_path: str = f"{constants.OUR_PATH}found.nc"
        not_found: str = f"{constants.OUR_PATH}not_found.nc"
        abs_path: str = str(tmp_path / Path(rel_path))
        # These two should resolve to abs_path because we've mocked INPUTDATA_PREFIX to tmp_path
        abs_path_dinlocroot: str = os.path.join("$DIN_LOC_ROOT", rel_path)
        abs_path_dinlocrootcurly: str = os.path.join("${DIN_LOC_ROOT}", rel_path)

    nc_paths = NcPaths()

    default_content = f"""
relative_path = '{nc_paths.rel_path}'
absolute_path='{nc_paths.abs_path}'
dlr = "{nc_paths.abs_path_dinlocroot}"
        dlr_in_curly_brackets    ='{nc_paths.abs_path_dinlocrootcurly}'   ! some comment
! but_not = {nc_paths.not_found}
"""

    def _create(content: str = default_content) -> str:
        mock_usernl_file_path = tmp_path / Path("user_nl_dummy")
        with open(mock_usernl_file_path, "w", encoding="utf-8") as f:
            f.write(content)
        return mock_usernl_file_path, nc_paths

    return _create


@pytest.fixture(name="example_progress")
def fixture_example_progress(tmp_path) -> json_io.NoNanFillValueProgress:
    """Some test data to represent a progress object"""
    progress_file = str(tmp_path / "progress.json")
    data_to_save = json_io.NoNanFillValueProgress(progress_file=progress_file)
    file1 = "/path/to/file1.nc"
    data_to_save[file1]["new_fill_values"] = {"var1": -999, "var2": constants.USER_REQ_DELETE}
    data_to_save[file1]["found_in_files"] = {"dummy.xml": {file1}}
    file2 = "/path/to/file2.nc"
    data_to_save[file2]["new_fill_values"] = {"var3": -999.0}
    data_to_save[file2]["found_in_files"] = {"user_nl_dummy": {file2}}
    print(f"{data_to_save=}")
    return data_to_save
