#!/usr/bin/env bash
#
# THROWAWAY -- delete together with Dockerfile.cprnc-probe.
#
# Builds the cprnc probe image and runs one system test in it, then says
# whether the case used the image's cprnc instead of building its own. Run it
# inside an interactive Casper allocation, because the probe image is not
# saved and dies with the allocation:
#   execcasper -A <account> -l select=1:ncpus=16:mem=256GB -l walltime=02:00:00
#   docker/ctsm-ci-derecho-gnu/probe-cprnc.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${here}/../.."

test_name=SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs
want_cprnc=/usr/local/bin/cprnc

echo "=== building the probe image"
DOCKERFILE=Dockerfile.cprnc-probe \
IMAGE_TAG=ctsm-ci-derecho-gnu:cprnc-probe \
CTSM_BUILD_NO_SAVE=1 \
    "${here}/build-on-casper.sh"

echo "=== running ${test_name} in it"
set +e
IMAGE_TAG=localhost/ctsm-ci-derecho-gnu:cprnc-probe \
    "${here}/run-sys-tests-in-container.sh" -t "${test_name}"
rc=$?
set -e
echo "run-sys-tests-in-container.sh exited ${rc}"

# Same scratch rule as container-common.sh.
scratch="${SCRATCH_DIR:-${SCRATCH:?neither SCRATCH_DIR nor SCRATCH is set}/cases_devcontainer}"
testroot="$(ls -td "${scratch}"/tests_*ct 2>/dev/null | head -1)"
case_dir="$(ls -d "${testroot}/${test_name}".* 2>/dev/null | head -1)"
if [ -z "${case_dir}" ]; then
    echo "RESULT: no case directory found under ${testroot:-${scratch}}; see the output above" >&2
    exit 1
fi
echo "=== case: ${case_dir}"

ok=1
cprnc_setting="$(grep -ho 'id="CCSM_CPRNC" value="[^"]*"' "${case_dir}"/env_*.xml | head -1)"
echo "CCSM_CPRNC in case: ${cprnc_setting:-<not found>}"
case "${cprnc_setting}" in
    *"\"${want_cprnc}\"") ;;
    *) echo "FAIL: case does not point CCSM_CPRNC at ${want_cprnc}"; ok=0 ;;
esac

if ls "${case_dir}"/bld/cprnc.bldlog.* >/dev/null 2>&1; then
    echo "FAIL: the case built cprnc itself: $(ls "${case_dir}"/bld/cprnc.bldlog.*)"
    ok=0
else
    echo "ok: no cprnc.bldlog in the case's bld/"
fi

echo "--- TestStatus"
cat "${case_dir}/TestStatus"
grep -q '^PASS .* SHAREDLIB_BUILD' "${case_dir}/TestStatus" \
    || { echo "FAIL: SHAREDLIB_BUILD did not pass"; ok=0; }

if [ "${ok}" = 1 ]; then
    echo "RESULT: probe passed (the test's own outcome is in TestStatus above)"
else
    echo "RESULT: probe FAILED"
    exit 1
fi
