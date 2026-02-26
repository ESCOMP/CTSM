"""System tests of functions in json_io.py (things that do involve filesystem)"""

# pylint: disable=too-few-public-methods

import os
from unittest.mock import patch

from ctsm.no_nans_in_inputs.json_io import NoNanFillValueProgress


class TestSaveProgress:
    """Test the save() function"""

    def test_save(self, example_progress: NoNanFillValueProgress):
        example_progress.save()
        assert os.path.exists(example_progress.progress_file)
        with open(example_progress.progress_file, "r", encoding="utf8") as f:
            contents = f.read()
            assert contents != "{}"


class TestLoadProgress:
    """Test the load() function"""

    def test_load_corrupted_json(self, tmp_path, capsys):
        """Test loading a corrupted JSON file returns an empty dict and warns."""
        progress_file = tmp_path / "progress.json"
        with open(progress_file, "w", encoding="utf-8") as f:
            f.write("this is not valid json")

        loaded_data = NoNanFillValueProgress(progress_file=str(progress_file))

        assert loaded_data == NoNanFillValueProgress()
        captured = capsys.readouterr()
        assert "Warning: Could not load progress file" in captured.err


class TestProgressFunctions:
    """Test the save() and load() functions together"""

    def test_save_and_load_progress(self, example_progress: NoNanFillValueProgress):
        """Test that data is saved and loaded correctly."""

        # Save data
        example_progress.save()

        # Load data
        with patch("builtins.input", return_value="y"):
            loaded_data = NoNanFillValueProgress(progress_file=example_progress.progress_file)

        # Check that loaded data is identical to saved data
        assert loaded_data == example_progress
