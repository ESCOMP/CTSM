#!/usr/bin/env bash
#PBS -N ctsm-devcontainer-case
#PBS -q casper
#PBS -l select=1:ncpus=8:mem=32GB
#PBS -l walltime=02:00:00
#PBS -j oe
#
# Create, build and run a CTSM case inside the ctsm-ci-derecho-gnu container.
#
# A thin wrapper around cime/scripts/create_newcase: every argument is passed
# straight through to it, and this script only adds the container plumbing
# (mounts, the $HOME/.cime macro copy) plus the case.setup / preview_namelists
# / check_input_data / case.build / case.submit chain afterwards.
#
# Because the machine config declares BATCH_SYSTEM=none, case.submit runs the
# model in the FOREGROUND rather than queueing it -- see
# cime/CIME/case/case_submit.py, which forces no_batch when the batch system
# type is "none". So this script wants a compute allocation: qsub it, or run
# it inside an existing interactive session. A single-point Ld1 case is small
# enough for the latter.
#
# Usage:
#   docker/ctsm-ci-derecho-gnu/run-case-in-container.sh --case NAME [create_newcase args]
# Example:
#   docker/ctsm-ci-derecho-gnu/run-case-in-container.sh \
#       --case brazil_test \
#       --compset IHistClm60Bgc --res 1x1_brazil --mpilib mpi-serial \
#       --run-unsupported
#
# --case may be a bare name, which lands in /cases (i.e. the host's
# $CASES_DIR), or an absolute container path, which is used as given.
#
# These defaults are injected only if you did not pass them yourself:
#   --machine container   CIME cannot detect the machine from a container hostname
#   --compiler gnu        the only compiler in the image
#   --output-root /scratch  bld/ and run/ land on the persistent host mount
#
# Overridable via environment (see container-common.sh for the full list):
#   IMAGE_TAG  INPUTDATA  CASES_DIR  SCRATCH_DIR  MAKE_J  MOUNT_OPTS  DRY_RUN
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=container-common.sh
source "${here}/container-common.sh"

set -u

ctsm_host_setup

# Read --case (without rewriting it) so we know where to cd for the
# setup/build/submit chain. The container's working directory is /cases, so a
# bare name given to create_newcase already lands there; we only need to
# resolve it to the same absolute path create_newcase will use.
args=("$@")
case_val=""
idx=0
while [ ${idx} -lt ${#args[@]} ]; do
    case "${args[idx]}" in
        --case=*) case_val="${args[idx]#--case=}" ;;
        --case)   case_val="${args[idx+1]:-}" ;;
    esac
    idx=$((idx + 1))
done

if [ -z "${case_val}" ]; then
    echo "ERROR: --case is required (create_newcase needs it, and this script" >&2
    echo "       needs it to know where to run case.build / case.submit)." >&2
    exit 2
fi

case "${case_val}" in
    /*) casedir="${case_val}" ;;
    *)  casedir="/cases/${case_val}" ;;
esac

# The same directory as seen from the host, for error messages. Only the part
# of the tree under the /cases mount has a host equivalent; an absolute --case
# pointing elsewhere in the container does not.
case "${casedir}" in
    /cases/*) host_casedir="${CTSM_CASES_DIR}/${casedir#/cases/}" ;;
    *)        host_casedir="(container path ${casedir}, outside the /cases mount)" ;;
esac

cime_args=()
ctsm_has_arg --machine "$@"     || cime_args+=(--machine container)
ctsm_has_arg --compiler "$@"    || cime_args+=(--compiler gnu)
ctsm_has_arg --output-root "$@" || cime_args+=(--output-root /scratch)
cime_args+=("$@")

echo "case       ${casedir}"
ctsm_print_host_dirs "${host_casedir}" "$(basename "${casedir}")"
echo "running    create_newcase ${cime_args[*]}"

set +e
ctsm_podman_run bash -s -- "${casedir}" "${cime_args[@]}" <<'INSIDE'
set -eu
casedir="$1"; shift

source /ctsm/docker/ctsm-ci-derecho-gnu/container-common.sh
ctsm_container_setup

echo "HOME=$HOME"

cd /cases
/ctsm/cime/scripts/create_newcase "$@"

cd "${casedir}"

# The container machine config sets DOUT_S_ROOT=$ENV{HOME}/archive/$CASE,
# which inside the container is /root/archive -- ephemeral. Redirect it onto
# the persistent mount. Single-quoted so $CASE reaches the XML unexpanded and
# CIME resolves it per-case, exactly as the machine config does.
./xmlchange DOUT_S_ROOT='/scratch/archive/$CASE'

./case.setup

# preview_namelists runs buildnml, which writes the Buildconf/*.input_data_list
# files that check_input_data reads. Without it there is nothing to check.
./preview_namelists

# check_input_data WITHOUT --download only reports; it never writes, so it is
# safe against the read-only inputdata mount and gives a clean list of what is
# missing. Exit 42 is a private signal to the host half of this script, which
# knows the host-side paths worth naming in the error message.
if ! ./check_input_data; then
    exit 42
fi

./case.build
./case.submit
INSIDE
rc=$?
set -e

if [ "${rc}" -eq 42 ]; then
    cat >&2 <<MSG

=====================================================================
ERROR: input data missing from the read-only inputdata mount.

  host:      ${CTSM_INPUTDATA}
  container: /opt/cesmdata/inputdata (read-only)

check_input_data listed the missing files above. The mount is read-only
on purpose, so CIME cannot fetch them; left to itself it would have
failed mid-run with a wall of "wget failed" errors.

FIRST SUSPECT: dangling symlinks rather than genuinely missing data.
Much of the inputdata tree is absolute symlinks into sibling campaign
collections, which resolve only if that tree is mounted too. Check the
"symlinks" line in this script's output; if it is absent, set
INPUTDATA_SYMLINK_ROOT (default /glade/campaign). Confirm by running
"readlink -f" on a named file's host path.

To fix, from a normal host shell (not this container), in the case at
  ${host_casedir}
run "./check_input_data --download", or stage the files by hand, then
re-run this script.
=====================================================================
MSG
    exit 1
fi

exit "${rc}"
