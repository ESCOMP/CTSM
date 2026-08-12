#!/usr/bin/env bash

# Fail on any non-zero exit code
set -e

cli_tool="$1"

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}/.."

msg="~~~~~ Check that -b works"
cmd="./build_docs -b _build -c"

. test/compose_test_cmd.sh
set -x
$cmd

exit 0
