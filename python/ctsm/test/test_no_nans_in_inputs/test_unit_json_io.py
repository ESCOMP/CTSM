"""Unit tests of functions in json_io.py (things that don't involve filesystem)"""

# pylint: disable=too-few-public-methods

from unittest.mock import patch
import pytest

from ctsm.no_nans_in_inputs.constants import USER_REQ_DELETE
from ctsm.no_nans_in_inputs.json_io import (
    _convert_fif_dict_sets,
    create_empty_progress_dict_onefile,
    save_progress,
    load_progress,
)

TEST_NC_ABSPATH = "/abs/path/abc123.nc"


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


class TestConvertFifDictSets:
    """Tests of _convert_fif_dict_sets()"""

    @pytest.mark.parametrize(
        "obj_in, expected",
        [
            ([TEST_NC_ABSPATH], {TEST_NC_ABSPATH}),
            ({TEST_NC_ABSPATH}, [TEST_NC_ABSPATH]),
        ],
    )
    def test_convert_fif_dict_sets_set2list(self, obj_in, expected):
        """Test converting set to list and vice versa"""

        progress_in = {TEST_NC_ABSPATH: create_empty_progress_dict_onefile()}
        file_containing = "user_nl_clm"
        progress_in[TEST_NC_ABSPATH]["found_in_files"][file_containing] = obj_in
        progress_out = _convert_fif_dict_sets(progress_in, dest_type=type(expected))
        assert progress_out[TEST_NC_ABSPATH]["found_in_files"][file_containing] == expected


class TestLoadProgress:
    """Test the _load_progress() function"""

    def test_load_nonexistent_file(self):
        """Test loading a nonexistent progress file returns an empty dict."""
        loaded_data = load_progress("/nonexistent/path/progress.json")
        assert loaded_data == {}


class TestSaveProgress:
    """Test the save_progress() function"""

    def test_save_io_error(self, example_data, capsys):
        """Test that an IOError during save is caught and a warning is printed."""
        with patch("builtins.open", side_effect=IOError("Permission denied")):
            save_progress(example_data, "/readonly/path/progress.json")
            captured = capsys.readouterr()
            assert "Warning: Could not save progress" in captured.err
