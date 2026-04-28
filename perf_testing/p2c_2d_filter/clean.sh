#!/usr/bin/env bash
# Remove build artifacts and per-run output. Does NOT delete source or
# baseline_checksum.txt.
set -eu
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"
rm -f driver *.o *.mod last_run.txt
