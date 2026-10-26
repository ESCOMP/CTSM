#!/usr/bin/env bash
#
# Build the image and run the full validation (the build-ctsm-ci-image skill's
# five wrappers, smoke-test-pfunit.sh, and the strings acceptance test on
# cesm.exe) in one go.
#
# Run it inside an INTERACTIVE Casper allocation: podman storage is node-local,
# so validation has to happen in the allocation that built the image. A
# qsub'ed build-on-casper.sh cannot validate.
#   execcasper -A <account> -l select=1:ncpus=16:mem=256GB -l walltime=06:00:00
#   docker/ctsm-ci-derecho-gnu/build-and-validate.sh [extra 'podman build' args]
#
# The build still saves the image tarball to GLADE (build-on-casper.sh's
# default), so a validation failure does not cost the build. A build failure
# stops here; every validation step runs even if an earlier one failed, and
# the summary at the end lists each result. Everything is also written to
# build-and-validate.log.<timestamp> in the repo root.
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../.." && pwd)"
cd "${repo}"

log="${repo}/build-and-validate.log.$(date +%Y%m%d%H%M%S)"
exec > >(tee -a "${log}") 2>&1
echo "log: ${log}"

scratch="${SCRATCH_DIR:-${SCRATCH:?neither SCRATCH_DIR nor SCRATCH is set}/cases_devcontainer}"
test_name=SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs
marker="$(mktemp)"
trap 'rm -f "${marker}"' EXIT

declare -a names=() results=()
record() { names+=("$1"); results+=("$2"); }
step() {
    local name="$1"; shift
    echo
    echo "======== ${name}"
    if "$@"; then record "${name}" PASS; else record "${name}" "FAIL (exit $?)"; fi
}

echo "======== build"
"${here}/build-on-casper.sh" "$@" || { echo "BUILD FAILED; not validating"; exit 1; }

step "smoke-test"      "${here}/smoke-test.sh"
step "pFUnit smoke-test" "${here}/smoke-test-pfunit.sh"
step "unit tests"      "${here}/run-unit-tests-in-container.sh"
step "create_test"     "${here}/run-test-in-container.sh"
step "run_sys_tests"   "${here}/run-sys-tests-in-container.sh" -t "${test_name}"
step "case"            "${here}/run-case-in-container.sh" \
    --case "brazil_$(date +%Y%m%d_%H%M%S)" --compset IHistClm60BgcQianRsGs \
    --res 1x1_brazil --run-unsupported

# The acceptance test, against the cesm.exe from run-test-in-container.sh
# (output-root /scratch -> ${scratch}/<TEST>.<testid>/bld). Only executables
# written by this run count, hence the marker.
echo
echo "======== strings acceptance test"
exe="$(find "${scratch}" -maxdepth 3 -path "*/${test_name}.*/bld/cesm.exe" -newer "${marker}" -printf '%T@ %p\n' 2>/dev/null \
       | sort -rn | head -1 | cut -d' ' -f2-)"
if [ -z "${exe}" ]; then
    echo "no cesm.exe from this run found under ${scratch}"
    record "strings test" "FAIL (no cesm.exe)"
else
    echo "cesm.exe: ${exe}"
    n_mpis=$(strings -a "${exe}" | grep -c /ctsm/libraries/mpi-serial)
    n_pio=$(strings -a "${exe}" | grep -c /ctsm/libraries/parallelio)
    n_local=$(strings -a "${exe}" | grep -c /usr/local)
    echo "/ctsm/libraries/mpi-serial: ${n_mpis} (must be 0)"
    echo "/ctsm/libraries/parallelio: ${n_pio} (must be 0)"
    echo "/usr/local:                 ${n_local} (must be nonzero)"
    if [ "${n_mpis}" -eq 0 ] && [ "${n_pio}" -eq 0 ] && [ "${n_local}" -gt 0 ]; then
        record "strings test" PASS
    else
        record "strings test" FAIL
    fi
fi

echo
echo "======== summary"
failed=0
for i in "${!names[@]}"; do
    printf '%-18s %s\n' "${names[$i]}" "${results[$i]}"
    [ "${results[$i]}" = PASS ] || failed=1
done
echo "log: ${log}"
exit "${failed}"
