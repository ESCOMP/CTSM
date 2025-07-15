#!/usr/bin/env bash

# Fail on any non-zero exit code
set -e

cli_tool="$1"

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}/.."

echo "~~~~~ Check that Makefile method works"
set -x
make SPHINXOPTS="-W --keep-going" BUILDDIR=${PWD}/_build html

exit 0
