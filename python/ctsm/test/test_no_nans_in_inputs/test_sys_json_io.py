"""System tests of functions in json_io.py (things that do involve filesystem)"""

# pylint: disable=too-few-public-methods

import pytest

from ctsm.no_nans_in_inputs.constants import USER_REQ_DELETE
from ctsm.no_nans_in_inputs.json_io import (
    save_progress,
    load_progress,
    create_empty_progress_dict_onefile,
)


@pytest.fixture(name="example_data")
def fixture_example_data():
    """Some test data to represent a progress object"""
    data_to_save = {}
    file1 = "/path/to/file1.nc"
    data_to_save[file1] = create_empty_progress_dict_onefile()
    data_to_save[file1]["new_fill_values"] = {"var1": -999, "var2": USER_REQ_DELETE}
    data_to_save[file1]["found_in_files"] = {"dummy.xml": {file1}}
    file2 = "/path/to/file2.nc"
    data_to_save[file2] = create_empty_progress_dict_onefile()
    data_to_save[file2]["new_fill_values"] = {"var3": -999.0}
    data_to_save[file2]["found_in_files"] = {"user_nl_dummy": {file2}}
    return data_to_save


class TestLoadProgress:
    """Test the _load_progress() function"""

    def test_load_corrupted_json(self, tmp_path, capsys):
        """Test loading a corrupted JSON file returns an empty dict and warns."""
        progress_file = tmp_path / "progress.json"
        with open(progress_file, "w", encoding="utf-8") as f:
            f.write("this is not valid json")

        loaded_data = load_progress(str(progress_file))

        assert loaded_data == {}
        captured = capsys.readouterr()
        assert "Warning: Could not load progress file" in captured.err


class TestProgressFunctions:
    """Test the save_progress and _load_progress functions together"""

    def test_save_and_load_progress(self, tmp_path, example_data):
        """Test that data is saved and loaded correctly."""
        progress_file = tmp_path / "progress.json"

        # Save data
        save_progress(example_data, str(progress_file))

        # Load data
        loaded_data = load_progress(str(progress_file))

        # Check that loaded data is identical to saved data
        assert loaded_data == example_data
