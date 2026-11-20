#!/usr/bin/env bash

# Fail on any non-zero exit code
set -e

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}/"

# Compare docs built with container vs. ctsm_pylib
./test_container_eq_ctsm_pylib.sh

# Check that Makefile method works
cd "${SCRIPT_DIR}/"
rm -rf _build
conda run --no-capture-output -n ctsm_pylib ./test_makefile_method.sh

exit 0
