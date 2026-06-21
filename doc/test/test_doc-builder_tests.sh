#!/usr/bin/env bash

# Fail on any non-zero exit code
set -e

cli_tool="${1:-conda}"

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}"

echo "~~~~~ Check that doc-builder tests pass"
cd ../doc-builder/test
set -x
${cli_tool} run --no-capture-output -n ctsm_pylib make test

exit 0
