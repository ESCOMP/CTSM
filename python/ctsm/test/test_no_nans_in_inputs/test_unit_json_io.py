"""Unit tests of functions in json_io.py (things that don't involve filesystem)"""

# pylint: disable=too-few-public-methods

from unittest.mock import patch
import pytest

from ctsm.no_nans_in_inputs.constants import USER_REQ_DELETE
from ctsm.no_nans_in_inputs.json_io import (
    _convert_fif_dict_sets,
    create_empty_progress_dict_onefile,
    NoNanFillValueProgress,
)

TEST_NC_ABSPATH = "/abs/path/abc123.nc"


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
    """Test loading"""

    def test_load_nonexistent_file(self):
        """Test loading a nonexistent progress file returns an empty dict."""
        loaded_data = NoNanFillValueProgress(progress_file="/nonexistent/path/progress.json")
        assert not loaded_data


class TestSaveProgress:
    """Test the save() function"""

    def test_save_io_error(self, example_progress: NoNanFillValueProgress, capsys):
        """Test that an IOError during save is caught and a warning is printed."""
        with patch("builtins.open", side_effect=IOError("Permission denied")):

            example_progress.save()
            captured = capsys.readouterr()
            assert "Warning: Could not save progress" in captured.err


class TestDoneAndCleanup:
    """Test the done_with_file() and cleanup() functions"""

    def test_done_and_cleanup(self, example_progress: NoNanFillValueProgress):
        """Test the done_with_file() and cleanup() functions"""
        # Save the file with both keys
        keys_orig = list(example_progress.keys())
        k0 = keys_orig[0]
        k1 = keys_orig[1]
        example_progress.save()

        # Mark the 0th key for deletion
        example_progress.done_with_file(k0)
        assert not example_progress[k0]

        # Call cleanup
        example_progress.cleanup()

        # Make sure the 0th key was deleted from the object
        assert k0 not in example_progress
        assert k1 in example_progress

        # Make sure the 0th key was deleted from the file
        result = NoNanFillValueProgress(
            progress_file=example_progress.progress_file, load_without_asking=True
        )
        assert result == example_progress
