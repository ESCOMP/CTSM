#!/usr/bin/env bash
#PBS -N ctsm-devcontainer-sys-tests
#PBS -q casper
#PBS -l select=1:ncpus=8:mem=96GB
#PBS -l walltime=12:00:00
#PBS -j oe
#
# Run CTSM's run_sys_tests inside the ctsm-ci-derecho-gnu container.
#
# A thin wrapper around ./run_sys_tests: every argument is passed straight
# through to it, and this script only adds the container plumbing (mounts, the
# $HOME/.cime macro copy) plus the flags that tell run_sys_tests it is in a
# container.
#
# Why this rather than run-test-in-container.sh, which already wraps
# create_test: run_sys_tests adds the bookkeeping -- a dated testroot,
# cs.status / cs.status.fails aggregation, a recorded SRCROOT_GIT_STATUS,
# baseline compare/generate, and retry. It is also what a CTSM developer
# actually types.
#
# The larger PBS reservation than the other wrappers is deliberate: a suite
# runs its tests serially in one container (BATCH_SYSTEM=none), and podman
# load alone peaks around 54-56 GB because the vfs driver duplicates layers.
#
# Usage:
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh [run_sys_tests args]
# Examples:
#   # one test by name
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
#       -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc
#   # derecho's mpi-serial suite, filtered to this image's one compiler
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial
#   # see what would run, without running it
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial --dry-run
#
# These defaults are injected only if you did not pass them yourself:
#   --machine-name container  picks up MACHINE_DEFAULTS["container"]: the
#                             no-batch launcher, /scratch as the testroot base,
#                             /scratch/baselines, and no account requirement
#   --wait                    run_sys_tests otherwise backgrounds create_test and
#                             returns, and podman would tear the container down
#                             while it was still running
#   --extra-create-test-args "--machine container --compiler gnu"
#                             run_sys_tests never passes --machine to
#                             create_test itself, only --xml-machine, so
#                             without this CIME tries to auto-detect the
#                             machine from the container hostname and fails.
#                             A value you pass is MERGED with this, not
#                             replaced.
# and, only when you asked for a suite with -s/--suite-name:
#   --xml-machine derecho     whose testlist to read. The container has no
#                             entries in cime_config/testdefs/testlist_clm.xml,
#                             and replicating derecho is the image's purpose.
#                             Override with $XML_MACHINE.
#   --suite-compiler gnu      the only compiler in the image
#
# Overridable via environment (see container-common.sh for the full list):
#   IMAGE_TAG  INPUTDATA  CASES_DIR  SCRATCH_DIR  MAKE_J  MOUNT_OPTS  DRY_RUN
#   XML_MACHINE  the machine whose testlist a suite is read from
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=container-common.sh
source "${here}/container-common.sh"

set -u

ctsm_host_setup

# run_sys_tests takes --extra-create-test-args as a SINGLE string, which it
# later splits. So a value the user passed has to be merged with ours rather
# than either one clobbering the other. Pull it out of the argument list here
# and re-add the combined value below.
injected_extra="--machine container --compiler gnu"
user_extra=""
args=()
skip_next=0
for arg in "$@"; do
    if [ "${skip_next}" -eq 1 ]; then
        user_extra="${arg}"
        skip_next=0
        continue
    fi
    case "${arg}" in
        --extra-create-test-args=*) user_extra="${arg#--extra-create-test-args=}" ;;
        --extra-create-test-args)   skip_next=1 ;;
        *)                          args+=("${arg}") ;;
    esac
done

cime_args=()
ctsm_has_arg --machine-name "${args[@]+"${args[@]}"}" \
    || cime_args+=(--machine-name container)
ctsm_has_arg --wait "${args[@]+"${args[@]}"}" \
    || cime_args+=(--wait)
cime_args+=(--extra-create-test-args "${injected_extra}${user_extra:+ ${user_extra}}")

# Suite-only injections. --xml-machine and --suite-compiler are meaningless
# without -s, and passing them anyway would just be confusing output.
if ctsm_has_arg -s "${args[@]+"${args[@]}"}" \
   || ctsm_has_arg --suite-name "${args[@]+"${args[@]}"}"; then
    ctsm_has_arg --xml-machine "${args[@]+"${args[@]}"}" \
        || cime_args+=(--xml-machine "${XML_MACHINE:-derecho}")
    ctsm_has_arg --suite-compiler "${args[@]+"${args[@]}"}" \
        || ctsm_has_arg --suite-compilers "${args[@]+"${args[@]}"}" \
        || cime_args+=(--suite-compiler gnu)
fi

cime_args+=("${args[@]+"${args[@]}"}")

# run_sys_tests names the testroot tests_<MMDD-HHMMSS><first 2 chars of the
# machine name>, so "co" here. Unlike the other wrappers there is no separate
# --test-root / --output-root: run_sys_tests passes --output-root <testroot> to
# create_test, so case dirs, bld/ and run/ all land inside the testroot.
echo "host-side directories:"
echo "  testroot   ${CTSM_SCRATCH_DIR}/tests_<MMDD-HHMMSS>co"
echo "             (cases, bld/ and run/ all land inside it)"
echo "  baselines  ${CTSM_SCRATCH_DIR}/baselines"
echo "note: run_sys_tests also links the testroot into ${CTSM_CASES_DIR},"
echo "      but that link points at the CONTAINER path and so is dangling"
echo "      when read from the host. Use the testroot path above."

log="${CTSM_SCRATCH_DIR}/run_sys_tests.$(date +%Y%m%d%H%M%S).log"
echo "log        ${log}"
echo "running    run_sys_tests ${cime_args[*]}"

set +e
ctsm_podman_run bash -s -- "${cime_args[@]}" <<'INSIDE' 2>&1 | tee "${log}"
set -eu

source /ctsm/docker/ctsm-ci-derecho-gnu/container-common.sh
ctsm_container_setup

echo "HOME=$HOME"

# run_sys_tests makes a symlink to the testroot in the working directory.
cd /cases

exec /ctsm/run_sys_tests "$@"
INSIDE
rc=${PIPESTATUS[0]}
set -e

# Same failure mode as run-test-in-container.sh: CIME cannot fetch missing
# inputdata through a read-only mount, and the first suspect is dangling
# symlinks rather than genuinely absent data. See that script for the full
# explanation of the wording matched here.
if [ "${rc}" -ne 0 ] \
   && grep -qE 'Could not find all inputdata on any server|wget failed|Read-only file system|Errno 30' "${log}"; then
    cat >&2 <<MSG

=====================================================================
ERROR: CIME could not find input data and failed trying to download it.

  host:      ${CTSM_INPUTDATA}
  container: /opt/cesmdata/inputdata (read-only)

FIRST SUSPECT: dangling symlinks, not genuinely missing data. Much of
the inputdata tree is absolute symlinks into sibling campaign
collections, which resolve only if that tree is mounted too. Check the
"symlinks" line in this script's output: if it is absent, set
INPUTDATA_SYMLINK_ROOT (default /glade/campaign) to the tree the
symlinks point into.

Full log: ${log}
=====================================================================
MSG
fi

exit "${rc}"
