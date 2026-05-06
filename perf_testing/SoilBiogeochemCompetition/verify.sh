#!/bin/bash
# Build the standalone driver and confirm both --fast and --all
# checksums MATCH their committed baselines.
#
# Usage: ./verify.sh [extra_make_args...]
#   ./verify.sh
#   ./verify.sh EXTRA_FFLAGS="-acc=multicore"
#   ./verify.sh EXTRA_FFLAGS="-acc=gpu -gpu=cc80"

set -euo pipefail

cd "$(dirname "$0")"

# shellcheck disable=SC1091
. ../env.sh >/dev/null 2>&1

build_log=$(mktemp)
trap 'rm -f "$build_log"' EXIT

make clean >/dev/null
if ! make "$@" >"$build_log" 2>&1; then
    echo "BUILD FAILED:"
    cat "$build_log"
    exit 1
fi

run_mode() {
    local label="$1"
    shift
    echo "=== $label ==="
    ./driver "$@" 2>&1 | grep -E '^\s+(checksum|baseline|elapsed|per call)\s*='
}

run_mode "--fast" --fast
run_mode "--all"
