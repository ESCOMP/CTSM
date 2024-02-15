# This script should be called in shell_commands with
#    . "${SRCROOT}"/cime_config/testdefs/testmods_dirs/clm/_scripts/get_fates_shell_commands_vars.sh
# where the leading period ensures it's run in the same shell.

CASE=`./xmlquery CASE --value`
FATESDIR="${SRCROOT}/src/fates"
FATESPARAMFILE="${SRCROOT}/src/fates/parameter_files/binaries/${CASE}-params.nc"

# No exit status because it should be called in the same shell.
