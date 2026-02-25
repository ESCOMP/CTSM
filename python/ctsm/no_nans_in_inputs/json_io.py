"""
Module handling the JSON file we use for saving progress and passing info between scripts
"""

import os
import sys
from copy import deepcopy
from typing import Dict, Type
import json

# Add the python directory to sys.path for direct script execution
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _CTSM_PYTHON not in sys.path:
    sys.path.insert(1, _CTSM_PYTHON)

from ctsm.no_nans_in_inputs.constants import (  # pylint: disable=wrong-import-position
    NEW_FILLVALUES_FILE,
)


def _convert_fif_dict_sets(progress: Dict, dest_type: Type) -> Dict:
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


def create_empty_progress_dict_onefile():
    """Return a dictionary for one netCDF file"""
    return {"found_in_files": {}, "new_fill_values": {}, "vars_with_nan_fills": []}


def _get_n_vars_in_progress(progress: dict) -> int:
    n_vars = 0
    for file in progress.keys():
        n_vars += len(progress[file]["new_fill_values"].keys())
    return n_vars


def init_progress() -> dict:
    """
    Initialize our progress file: Either load an existing one or start a new one
    """
    progress = load_progress(NEW_FILLVALUES_FILE)
    if progress:
        print(f"\nLoaded progress from {NEW_FILLVALUES_FILE}")
        total_vars = _get_n_vars_in_progress(progress)
        print(f"Already processed {total_vars} variable(s) in {len(progress)} file(s)")
        response = input("Continue from where you left off? [Y/n]: ").strip().lower()
        if response and response not in ("y", "yes"):
            progress = {}
            print("Starting fresh...")
    return progress


def load_progress(progress_file: str) -> dict:
    """
    Load progress from a JSON file if it exists.

    Args:
        progress_file: Path to the progress file

    Returns:
        Previously saved progress, or empty dict if file doesn't exist
    """
    if not os.path.exists(progress_file):
        return {}

    try:
        with open(progress_file, "r", encoding="utf-8") as f:
            progress = json.load(f)

            # This is serialized as a list, but the code needs it as a set
            progress = _convert_fif_dict_sets(progress, set)

            return progress
    except (IOError, OSError, json.JSONDecodeError) as e:
        print(f"Warning: Could not load progress file: {e}", file=sys.stderr)
        return {}


def save_progress(progress: dict, progress_file: str) -> None:
    """
    Save progress to a JSON file.

    Args:
        progress: Dictionary of found locations and collected fill values
        progress_file: Path to the progress file
    """
    # Can't serialize sets. deepcopy() is needed so that caller's progress isn't affected (.copy()
    # isn't sufficient since we have nested mutables).
    progress_out = _convert_fif_dict_sets(deepcopy(progress), list)

    try:
        with open(progress_file, "w", encoding="utf-8") as f:
            json.dump(progress_out, f, indent=2)
        print(f"  [Progress saved to {progress_file}]")
    except (IOError, OSError) as e:
        print(f"  Warning: Could not save progress: {e}", file=sys.stderr)
