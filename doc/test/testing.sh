#!/usr/bin/env bash

# Fail on any non-zero exit code
set -e

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}/"

# Compare docs built with container vs. ctsm_pylib
./test_container_eq_ctsm_pylib.sh

# Check that -r -v works (Docker)
# Also do a custom --conf-py-path and other stuff
cd "${SCRIPT_DIR}/"
rm -rf _build
./test_build_docs_-r-v.sh docker

# Check that Makefile method works
cd "${SCRIPT_DIR}/"
rm -rf _build
conda run --no-capture-output -n ctsm_pylib ./test_makefile_method.sh

# Check that -b works
cd "${SCRIPT_DIR}/"
rm -rf _build
./test_build_docs_-b.sh docker

# Check that doc-builder tests pass
# Don't run if on a GitHub runner; failing 🤷. Trust that doc-builder does this test.
if [[ "${GITHUB_ACTIONS}" == "" ]]; then
    cd "${SCRIPT_DIR}/"
    ./test_doc-builder_tests.sh
fi

exit 0
