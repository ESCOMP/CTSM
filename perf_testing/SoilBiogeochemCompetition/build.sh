#!/bin/bash
# Source the project env file, then `make clean && make "$@"`.
# Use this instead of inlining `. ../env.sh && make ...` in shell commands —
# scripted entry points stay whitelistable across runs.
#
# Usage:
#   ./build.sh                                        # serial
#   ./build.sh EXTRA_FFLAGS="-mp"                     # OpenMP
#   ./build.sh EXTRA_FFLAGS="-acc=gpu -gpu=cc80"      # GPU
#   ./build.sh EXTRA_FFLAGS="-acc=gpu -gpu=cc80" INNER_TIMING=1
#
# Filter output (grep, tail, head, etc.) at the call site, not in here.

set -euo pipefail
cd "$(dirname "$0")"

# shellcheck disable=SC1091
. ../env.sh >/dev/null 2>&1

make clean
make "$@"
