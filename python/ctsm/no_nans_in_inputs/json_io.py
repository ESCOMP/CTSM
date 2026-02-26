"""
Module handling the JSON file we use for saving progress and passing info between scripts
"""

from pathlib import Path
import os
import sys
from copy import deepcopy
from typing import Type
import json
from collections import defaultdict

# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)


def create_empty_progress_dict_onefile():
    """Return a dictionary for one netCDF file"""
    return {"found_in_files": {}, "new_fill_values": {}, "vars_with_nan_fills": []}


class NoNanFillValueProgress(defaultdict):
    """Defaultdict-like for tracking progress in getting/replacing NaN fill values"""

    def __init__(
        self,
        default_factory=create_empty_progress_dict_onefile,
        progress_file: str = "fkqenreqorni",
        load_without_asking: bool = False,
    ):
        """
        Initialize our progress file: Either load an existing one or start a new one
        """
        super().__init__(default_factory)

        if isinstance(progress_file, Path):
            progress_file = str(progress_file)
        self.progress_file = progress_file

        if progress_file and os.path.exists(progress_file):
            try:
                with open(progress_file, "r", encoding="utf-8") as f:
                    progress = json.load(f)

                    # This is serialized as a list, but the code needs it as a set
                    progress = _convert_fif_dict_sets(progress, set)

                    print(f"\nLoaded progress from {progress_file}")
                    total_vars = _get_n_vars_in_progress(self)
                    print(f"Already processed {total_vars} variable(s) in {len(self)} file(s)")
                    if load_without_asking:
                        self.update(progress)
                    else:
                        response = (
                            input("Continue from where you left off? [Y/n]: ").strip().lower()
                        )
                        if response and response not in ("y", "yes"):
                            print("Starting fresh...")
                        else:
                            self.update(progress)
            except (IOError, OSError, json.JSONDecodeError) as e:
                print(f"Warning: Could not load progress file: {e}", file=sys.stderr)

    def __setitem__(self, key, value):
        """Ensure all keys are strings"""
        super().__setitem__(str(key), value)

    def update(self, *args, **kwargs):
        """Convert keys to str for update operations"""
        # Handle dict or iterable of key/value pairs
        if args:
            other = args[0]
            if hasattr(other, "items"):
                for k, v in other.items():
                    self[str(k)] = v
            else:  # iterable of (k, v)
                for k, v in other:
                    self[str(k)] = v
        for k, v in kwargs.items():
            self[str(k)] = v

    def save(self) -> None:
        """
        Save progress to a JSON file.
        """
        # Can't serialize sets. deepcopy() is needed so that caller's progress isn't affected
        # .copy() isn't sufficient since we have nested mutables.
        progress_out = _convert_fif_dict_sets(deepcopy(self), list)

        try:
            with open(self.progress_file, "w", encoding="utf-8") as f:
                json.dump(progress_out, f, indent=2)
            print(f"  [Progress saved to {self.progress_file}]")
        except (IOError, OSError) as e:
            print(f"  Warning: Could not save progress: {e}", file=sys.stderr)

    def done_with_file(self, netcdf_path: str) -> None:
        """After we're done with a netCDF file, mark for removal from progress object/file"""
        if netcdf_path not in self:
            raise KeyError(netcdf_path)
        self[netcdf_path] = None

    def cleanup(self) -> None:
        """Remove keys marked for deletion, then update progress file"""
        keys_to_remove = [k for k in self if not self[k]]
        for key in keys_to_remove:
            self.pop(key)
        self.save()


def _convert_fif_dict_sets(
    progress: NoNanFillValueProgress, dest_type: Type
) -> NoNanFillValueProgress:
    """
    The code needs the "found_in_files" dictionary to contain sets, but the JSON serializer can only
    handle lists. This function allows the conversion of items in that dictionary between lists and
    sets.

    Args:
        progress: Dictionary of found locations and collected fill values
        dest_type: Type to convert to
    """
    for abs_path in progress:
        fif_dict = progress[abs_path]["found_in_files"]
        for file_containing in fif_dict:
            fif_dict[file_containing] = dest_type(fif_dict[file_containing])
    return progress


def _get_n_vars_in_progress(progress: NoNanFillValueProgress) -> int:
    n_vars = 0
    for file in progress.keys():
        n_vars += len(progress[file]["new_fill_values"].keys())
    return n_vars
