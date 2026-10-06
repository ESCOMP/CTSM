#!/usr/bin/env bash
#PBS -N ctsm-ci-derecho-gnu-build
#PBS -q casper
#PBS -l select=1:ncpus=16:mem=256GB
#PBS -l walltime=03:00:00
#PBS -j oe
#
# Build the ctsm-ci-derecho-gnu image with podman on an NCAR HPC node (Casper).
#
# This wrapper exists because of a host-side requirement that CANNOT live
# in the Dockerfile: rootless podman/buildah must place the build
# container's rootfs on a node-local filesystem. When TMPDIR points at a
# parallel filesystem (e.g. glade scratch) the build fails while creating
# the rootfs, before any Dockerfile instruction runs, e.g.:
#   creating directory ".../buildahNNN/mnt/rootfs": permission denied
# TMPDIR is read by buildah in the host shell, so it must be exported here
# rather than set via ENV inside the image. See:
#   https://ncar-hpc-docs.readthedocs.io/en/latest/environment-and-software/user-environment/containers/working_with_containers/
#
# NCAR does not provision subuid/subgid ranges (rootless podman runs in
# single-UID mapping); that is fine for this image because the build never
# switches users or chowns to other UIDs.
#
# Usage, interactively on a compute node (e.g. after execcasper):
#   docker/ctsm-ci-derecho-gnu/build-on-casper.sh [extra 'podman build' args]
#   docker/ctsm-ci-derecho-gnu/build-on-casper.sh --build-arg MAKE_JOBS=8
# or as a batch job, submitted FROM THE REPO ROOT so PBS_O_WORKDIR locates the
# build context (PBS runs a copy of this script from its own spool directory):
#   qsub -A <account> docker/ctsm-ci-derecho-gnu/build-on-casper.sh
# Overridable via environment:
#   CTSM_BUILD_CONTEXT build context dir     (default: wherever this file is,
#                                             else found via PBS_O_WORKDIR)
#   CTSM_BUILD_TMPDIR  node-local scratch dir (default /var/tmp/$USER)
#   IMAGE_TAG          image tag to build     (default ctsm-ci-derecho-gnu:dev)
#   DOCKERFILE         Dockerfile to use      (default Dockerfile)
set -eo pipefail

if [ "$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    batch=1
    logdest=/dev/null
else
    batch=0
    logdest="build-on-casper.log.$(date +%Y%m%d%H%M%S%N)"
fi

set -u

module load podman

user="${USER:-$(id -un)}"
dockerfile="${DOCKERFILE:-Dockerfile}"

# Where the build context lives. Normally that is this script's own directory,
# but under `qsub` PBS executes a COPY of the script from its spool area
# (/var/spool/pbs/mom_priv/jobs), so BASH_SOURCE points there and nothing of
# this repo is alongside it. Recover the real directory from PBS_O_WORKDIR,
# which is wherever qsub was invoked -- accepting either the repo root or this
# directory. CTSM_BUILD_CONTEXT overrides for anything stranger.
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -f "${here}/${dockerfile}" ]; then
    for cand in "${CTSM_BUILD_CONTEXT:-}" \
                "${PBS_O_WORKDIR:-}/docker/ctsm-ci-derecho-gnu" \
                "${PBS_O_WORKDIR:-}"; do
        if [ -n "${cand}" ] && [ -f "${cand}/${dockerfile}" ]; then
            here="$(cd "${cand}" && pwd)"
            break
        fi
    done
fi
if [ ! -f "${here}/${dockerfile}" ]; then
    echo "ERROR: cannot find ${dockerfile}. Looked next to this script" >&2
    echo "       (${here}) and under PBS_O_WORKDIR=${PBS_O_WORKDIR:-<unset>}." >&2
    echo "       Submit with qsub from the repo root, or set" >&2
    echo "       CTSM_BUILD_CONTEXT to this directory." >&2
    exit 1
fi

# Force a node-local TMPDIR (ignore any inherited parallel-FS value).
TMPDIR="${CTSM_BUILD_TMPDIR:-/var/tmp/${user}}"
export TMPDIR
mkdir -p "${TMPDIR}"

image="${IMAGE_TAG:-ctsm-ci-derecho-gnu:dev}"

echo "Building ${image} from ${dockerfile} (TMPDIR=${TMPDIR})"
podman build \
    -f "${here}/${dockerfile}" \
    -t "${image}" \
    "$@" \
    "${here}" 2>&1 | tee -p "${logdest}"
