#!/usr/bin/env bash
set -e

# Check that clm6* compset aliases return CLM6* longnames

# Change to top level of clone
cd "$(git rev-parse --show-toplevel)"

# Check that query_config can run without error
cime/scripts/query_config --compsets 1>/dev/null

# Find bad compsets
OLD_IFS=$IFS
IFS='\n'
set +e
# Relies on case sensitivity here: Alias should have Clm6 and longname should have CLM6
bad_compsets="$(cime/scripts/query_config --compsets | sort | uniq | grep Clm6 | grep -v CLM6)"
set -e
if [[ "${bad_compsets}" != "" ]]; then
    echo "One or more compsets with Clm6 alias but not CLM6 longname:" >&2
    echo $bad_compsets  >&2
    exit 1
fi

exit 0