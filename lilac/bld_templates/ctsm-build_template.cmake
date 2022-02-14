# This, together with config_machines_template.xml, provides a machine port for building
# CTSM.
#
# If you are looking at the template file: Variable names prefixed with a dollar sign will
# be replaced with machine-specific values. A double dollar sign gets replaced with a
# single dollar sign, so something like $$MYVAR refers to the MYVAR cime variable.

if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " $GPTL_CPPDEFS")
endif()

set(NETCDF_PATH "$NETCDF_PATH")

# If PIO_FILESYSTEM_HINTS is provided, this will set a PIO_FILESYSTEM_HINTS variable; if
# not provided, this will just be a blank line.
$PIO_FILESYSTEM_HINTS

# If PNETCDF_PATH is provided, this will set a PNETCDF_PATH variable; if not provided,
# this will just be a blank line.
$PNETCDF_PATH

string(APPEND CFLAGS " $EXTRA_CFLAGS")
string(APPEND FFLAGS " $EXTRA_FFLAGS")
