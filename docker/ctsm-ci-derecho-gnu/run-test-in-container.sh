#!/usr/bin/env bash
#PBS -N ctsm-devcontainer-test
#PBS -q casper
#PBS -l select=1:ncpus=8:mem=32GB
#PBS -l walltime=02:00:00
#PBS -j oe
#
# Run CTSM system tests inside the ctsm-ci-derecho-gnu container.
#
# A thin wrapper around cime/scripts/create_test: every argument is passed
# straight through to it, and this script only adds the container plumbing
# (mounts, the $HOME/.cime macro copy).
#
# Unlike the older CI job, this does NOT pass --no-run: the tests build AND
# run, which is what the read-only inputdata mount is for.
#
# BATCHING: none is needed, and none should be requested. CIME infers no-batch
# mode from the machine config's BATCH_SYSTEM=none --
# cime/CIME/test_scheduler.py sets
#   self._no_batch = no_batch or not self._machobj.has_batch_system()
# which also makes create_test block until the tests actually finish and
# report a real PASS/FAIL (rather than "submitted"), so --wait is redundant.
# Do NOT add --queue: test_scheduler.py asserts that a queue is never combined
# with no-batch mode, so it is a hard error here.
#
# The run happens in the foreground, so this wants a compute allocation: qsub
# it, or run it inside an existing interactive session.
#
# Usage:
#   docker/ctsm-ci-derecho-gnu/run-test-in-container.sh [create_test args]
# Examples:
#   # the default single-point test
#   docker/ctsm-ci-derecho-gnu/run-test-in-container.sh
#   # a specific test
#   docker/ctsm-ci-derecho-gnu/run-test-in-container.sh \
#       SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc
#   # a whole suite, filtered to the one compiler this image has
#   docker/ctsm-ci-derecho-gnu/run-test-in-container.sh \
#       --xml-category prealpha --xml-machine derecho
#
# These defaults are injected only if you did not pass them yourself:
#   --machine container    CIME cannot detect the machine from a container hostname
#   --compiler gnu         the only compiler in the image. Fills in the compiler
#                          for tests named on the command line; a test name that
#                          encodes one (..._gnu) still wins.
#   --xml-compiler gnu     the filter for --xml-* suite queries, which is a
#                          DIFFERENT knob from --compiler. Ignored unless the
#                          query branch is taken (create_test.py takes
#                          `if args.testargs:` whenever test names are given),
#                          so injecting both is safe.
#   --test-root /cases     case dirs land on the persistent host mount
#   --output-root /scratch bld/ and run/ land on the persistent host mount
#
# With no arguments at all, runs DEFAULT_TEST (below).
#
# Overridable via environment (see container-common.sh for the full list):
#   IMAGE_TAG  INPUTDATA  CASES_DIR  SCRATCH_DIR  MAKE_J  MOUNT_OPTS  DRY_RUN
#   DEFAULT_TEST  the test used when no arguments are given
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=container-common.sh
source "${here}/container-common.sh"

set -u

# mpi-serial, matching how CTSM's own single-point tests are written in
# cime_config/testdefs/testlist_clm.xml (they all use _Mmpi-serial). CIME
# builds mpi-serial from CTSM's libraries/mpi-serial submodule, forces
# NTASKS=1 -- which a 1x1 grid needs, since it cannot decompose over the
# container machine's four default tasks -- and selects PIO_TYPENAME=netcdf.
#
# Mirrors testlist_clm.xml's SMS_Ld5_Mmpi-serial.1x1_brazil.IHistClm60Bgc,
# minus its testmods and shortened to one day.
default_test="${DEFAULT_TEST:-SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc}"

ctsm_host_setup

cime_args=()
ctsm_has_arg --machine "$@"      || cime_args+=(--machine container)
ctsm_has_arg --compiler "$@"     || cime_args+=(--compiler gnu)
ctsm_has_arg --xml-compiler "$@" || cime_args+=(--xml-compiler gnu)
ctsm_has_arg --test-root "$@"    || cime_args+=(--test-root /cases)
ctsm_has_arg --output-root "$@"  || cime_args+=(--output-root /scratch)

if [ "$#" -eq 0 ]; then
    echo "no arguments given; using default test ${default_test}"
    cime_args+=("${default_test}")
else
    cime_args+=("$@")
fi

ctsm_print_host_dirs "${CTSM_CASES_DIR}/<TEST>.<testid>" "<TEST>.<testid>"

log="${CTSM_SCRATCH_DIR}/create_test.$(date +%Y%m%d%H%M%S).log"
echo "log        ${log}"
echo "running    create_test ${cime_args[*]}"

set +e
ctsm_podman_run bash -s -- "${cime_args[@]}" <<'INSIDE' 2>&1 | tee "${log}"
set -eu

source /ctsm/docker/ctsm-ci-derecho-gnu/container-common.sh
ctsm_container_setup

echo "HOME=$HOME"

exec /ctsm/cime/scripts/create_test "$@"
INSIDE
rc=${PIPESTATUS[0]}
set -e

# create_test drives its own phases, so unlike run-case-in-container.sh there
# is no place to slip a read-only check_input_data in first. Recognize the
# failure after the fact instead. Observed 2026-08-31, the wording that
# actually appears is CIME's download failing per file --
#   wget failed with output:  and errput <path>: No such file or directory
#   ERROR: Could not find all inputdata on any server
# The os.makedirs() under DIN_LOC_ROOT that a :ro mount refuses (Errno 30) is
# also matched, since which one surfaces depends on whether the parent
# directory already exists. Falls through to the raw error otherwise.
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
symlinks point into. To confirm, pick a file the log names and run
"readlink -f" on its host path -- if it resolves outside
${CTSM_INPUTDATA}, that is the problem.

Otherwise the data really is missing. From a normal host shell (not
this container), in the failed case under
  ${CTSM_CASES_DIR}
run "./check_input_data" (no --download) for the list, then
"./check_input_data --download" to fetch it, and re-run this script.

Full log: ${log}
=====================================================================
MSG
fi

exit "${rc}"
