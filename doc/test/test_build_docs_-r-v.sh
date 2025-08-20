#!/usr/bin/env bash

# Fail on any non-zero exit code
set -e

cli_tool="$1"

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}/.."

msg="~~~~~ Check that -r -v works"
cmd="./build_docs -r _build -v latest -c --conf-py-path doc-builder/test/conf.py --static-path ../_static --templates-path ../_templates"

. test/compose_test_cmd.sh
set -x
$cmd

exit 0
