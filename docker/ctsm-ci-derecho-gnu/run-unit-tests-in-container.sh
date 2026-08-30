#!/usr/bin/env bash
# Step 2 of the pFUnit probe: run CTSM's own Fortran unit tests inside the
# container, which is what actually answers "is the pFUnit in the image
# usable by CIME?". Step 1 (smoke-test-pfunit.sh) only proves the mechanism.
#
# Usage:
#   docker/ctsm-ci-derecho-gnu/run-unit-tests-in-container.sh [extra run_tests.py args]
# Overridable via environment:
#   IMAGE_TAG   image to run    (default localhost/ctsm-ci-derecho-gnu-pfunit:probe)
#   BUILD_DIR   build dir INSIDE the container (default /tmp/unit_tests.temp).
#               Defaults to container-local scratch rather than the mounted
#               repo: the build writes many small files and glade is a
#               parallel FS. Point it at /ctsm/src/unit_tests.temp if you want
#               the incremental rebuild that src/README.unit_testing describes.
#   MOUNT_OPTS  extra flags on the repo mount, e.g. ":z" if SELinux relabeling
#               turns out to be needed (default none)
#   MAKE_J      build parallelism (default 16)
set -eo pipefail

module load podman 2>/dev/null || true

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../.." && pwd)"

image="${IMAGE_TAG:-localhost/ctsm-ci-derecho-gnu-pfunit:probe}"
build_dir="${BUILD_DIR:-/tmp/unit_tests.temp}"

echo "Running CTSM unit tests in ${image}"
echo "  repo      ${repo} -> /ctsm"
echo "  build-dir ${build_dir} (inside the container)"

podman run --rm -i \
    -v "${repo}:/ctsm${MOUNT_OPTS:-}" \
    -w /ctsm/src \
    -e "MAKE_J=${MAKE_J:-16}" \
    "${image}" \
    bash -s -- "${build_dir}" "$@" <<'INSIDE'
set -eu
build_dir="$1"; shift

# GitHub Actions container jobs override HOME, so the /root/.cime copy baked
# into the image is not necessarily where CIME will look. Place it wherever
# HOME actually points. This is the same two lines the CI job will need.
echo "HOME=$HOME"
if [ ! -f "$HOME/.cime/gnu_container.cmake" ]; then
    echo "installing gnu_container.cmake into $HOME/.cime"
    mkdir -p "$HOME/.cime"
    cp /opt/ctsm-container/cime-macros/gnu_container.cmake "$HOME/.cime/"
fi

# --machine container: CIME cannot detect the machine from the hostname in
#   here, and the container machine is what the image is built to satisfy.
# --enable-genf90: CTSM has 12 .F90.in templates whose generated .F90 files
#   are NOT committed. Without this, genf90_utils.cmake sets the output
#   variables "as if the generation had occurred" but generates nothing, and
#   the build fails on missing sources. The flag also makes run_tests.py pass
#   -DCMAKE_PROGRAM_PATH=<cime>/CIME/non_py/externals/genf90 so find_program
#   locates genf90.pl, which is otherwise not on PATH.
# --test-spec-dir defaults to "." (i.e. /ctsm/src), so it is left implicit,
#   matching src/README.unit_testing.
exec ../cime/scripts/fortran_unit_testing/run_tests.py \
    --build-dir "$build_dir" \
    --machine container \
    --compiler gnu \
    --enable-genf90 \
    --make-j "${MAKE_J:-16}" \
    "$@"
INSIDE
