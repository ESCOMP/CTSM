"""
Shared constants for the no_nans_in_inputs module.
"""

# File paths
NEW_FILLVALUES_FILE = "new_fillvalues.json"  # File to save/load new fill values
XML_FILE = "bld/namelist_files/namelist_defaults_ctsm.xml"  # CTSM namelist defaults XML
OUR_PATH = "lnd/clm2/"  # String to be found in files we're responsible for

# Filename suffix after fixing NaN fills
NONANFILL_SUFFIX = "no_nan_fill"

# NetCDF attribute name
ATTR = "_FillValue"

# Special commands for user input
USER_REQ_QUIT = "quit"
USER_REQ_SKIP_VAR = "skip"
USER_REQ_SKIP_FILE = "skipfile"
USER_REQ_DELETE = "delete"

# Error strings corresponding to special user commands
ERR_STR_SKIP_VAR = "SKIP_VARIABLE"
ERR_STR_SKIP_FILE = "SKIP_FILE"

# Output formatting
SEP_LENGTH = 80  # Length of horizontal separators in stdout

# Keyword arguments we want to include in every xarray.open_dataset() call.
OPEN_DS_KWARGS = {"decode_timedelta": False, "decode_times": False}
