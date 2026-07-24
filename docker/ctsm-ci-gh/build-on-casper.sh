#!/usr/bin/env bash
#PBS -N ctsm-ci-gh-build
#PBS -q casper
#PBS -l select=1:ncpus=16:mem=256GB
#PBS -l walltime=03:00:00
#PBS -j oe
#
# Build the ctsm-ci-gh image with podman on an NCAR HPC node (Casper).
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
# Usage:
#   docker/ctsm-ci-gh/build-on-casper.sh [extra 'podman build' args]
# Examples:
#   docker/ctsm-ci-gh/build-on-casper.sh
#   docker/ctsm-ci-gh/build-on-casper.sh --build-arg MAKE_JOBS=8
# Overridable via environment:
#   CTSM_BUILD_TMPDIR  node-local scratch dir (default /var/tmp/$USER)
#   IMAGE_TAG          image tag to build     (default ctsm-ci-gh:dev)
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

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
user="${USER:-$(id -un)}"

# Force a node-local TMPDIR (ignore any inherited parallel-FS value).
TMPDIR="${CTSM_BUILD_TMPDIR:-/var/tmp/${user}}"
export TMPDIR
mkdir -p "${TMPDIR}"

image="${IMAGE_TAG:-ctsm-ci-gh:dev}"
dockerfile="${DOCKERFILE:-Dockerfile}"

echo "Building ${image} from ${dockerfile} (TMPDIR=${TMPDIR})"
podman build \
    -f "${here}/${dockerfile}" \
    -t "${image}" \
    "$@" \
    "${here}" 2>&1 | tee -p "${logdest}"
