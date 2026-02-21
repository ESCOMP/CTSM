"""
Shared constants for the no_nans_in_inputs module.
"""

# File paths
NEW_FILLVALUES_FILE = "new_fillvalues.json"  # File to save/load new fill values
XML_FILE = "bld/namelist_files/namelist_defaults_ctsm.xml"  # CTSM namelist defaults XML

# NetCDF attribute name
ATTR = "_FillValue"

# Special commands for user input
USER_REQ_QUIT = "quit"
USER_REQ_SKIP_VAR = "skip"
USER_REQ_SKIP_FILE = "skipfile"
USER_REQ_DELETE = "delete"
