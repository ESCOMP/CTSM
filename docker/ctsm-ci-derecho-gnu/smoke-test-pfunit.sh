#!/usr/bin/env bash
# Smoke-test the pFUnit half of the ctsm-ci-derecho-gnu image. Verifies that
# pFUnit installed where the CIME macro drop-in says it did, that a real pFUnit
# test
# preprocesses/compiles/links/runs, and -- the part that actually answers the
# design question -- that PFUNIT_PATH survives the variable filter in
# ccs_config's Macros.cmake, which is how run_tests.py's find_pfunit() sees it.
#
# This needs no CTSM checkout: it exercises the mechanism, not the CTSM tests.
# Run it alongside smoke-test.sh, which covers the rest of the stack; then
# run-unit-tests-in-container.sh for CTSM's own unit tests.
#
# Usage: docker/ctsm-ci-derecho-gnu/smoke-test-pfunit.sh
#   IMAGE_TAG overrides the image (default localhost/ctsm-ci-derecho-gnu:dev)
set -eo pipefail

module load podman 2>/dev/null || true

image="${IMAGE_TAG:-localhost/ctsm-ci-derecho-gnu:dev}"
echo "Smoke-testing ${image}"

# Runs INSIDE a fresh container so we exercise exactly the baked-in
# environment (non-login shell, image ENV only) that GitHub Actions sees.
podman run --rm -i "${image}" bash -s <<'INSIDE'
set -u
fail() { echo "PFUNIT SMOKE FAIL: $*" >&2; exit 1; }

macro=/opt/ctsm-container/cime-macros/gnu_container.cmake

echo "### macro drop-in is present"
[ -f "$macro" ] || fail "missing $macro"
echo "HOME=$HOME  (GitHub Actions overrides this; see the macro file)"
[ -f "$HOME/.cime/gnu_container.cmake" ] \
    || fail "missing \$HOME/.cime/gnu_container.cmake"

echo "### pFUnit install layout"
pf="$(sed -n 's/^ *set(PFUNIT_PATH "\([^"]*\)").*/\1/p' "$macro" | tail -1)"
echo "PFUNIT_PATH (macro default) = $pf"
# PFUNIT_PATH must name the SHARED PREFIX holding all four sibling packages,
# not the PFUNIT-4.8 subdirectory. PFUNITConfig.cmake finds its dependencies
# only via GFTL_ROOT / GFTL_SHARED_ROOT / FARGPARSE_ROOT, which CMake honors
# only under CMP0074 -- and CTSM's src/CMakeLists.txt declares
# cmake_minimum_required(VERSION 2.8), leaving that policy unset. Checking
# only PFUNIT here would pass while the real CTSM build failed on GFTL.
for pkg in PFUNIT GFTL GFTL_SHARED FARGPARSE; do
    cfg="$(ls -d "$pf"/${pkg}-*/cmake/${pkg}Config.cmake 2>/dev/null | head -1)"
    [ -n "$cfg" ] || fail "no ${pkg}Config.cmake under $pf/${pkg}-*/cmake/"
    echo "  $pkg -> $cfg"
done

pflib="$(ls -d "$pf"/PFUNIT-*/lib 2>/dev/null | head -1)"
[ -n "$pflib" ] || fail "no PFUNIT-*/lib under $pf"
funitproc="$(ls -d "$pf"/PFUNIT-*/bin/funitproc 2>/dev/null | head -1)"
[ -n "$funitproc" ] && [ -x "$funitproc" ] || fail "no funitproc under $pf/PFUNIT-*/bin"
ls "$pflib" | head

echo "### pFUnit is a serial build (no MPI symbols linked in)"
if ls "$pflib" | grep -qi 'mpi'; then
    fail "found an MPI-flavored pFUnit library; expected SKIP_MPI=YES"
fi

echo "### CIME's Macros.cmake filter still lets PFUNIT_PATH through"
# ccs_config/machines/cmake_macros/Macros.cmake snapshots the CMake variable
# set before including the macro files and emits only what appeared after.
# That is why -DPFUNIT_PATH=... and an env var alone do NOT work, and it is
# the exact mechanism find_pfunit() scrapes. Replicate it here.
mkdir -p /tmp/macrocheck && cd /tmp/macrocheck
cat > CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.5)
project(cime LANGUAGES NONE)
get_cmake_property(VARS_BEFORE_BUILD_INTERNAL_IGNORE VARIABLES)
include(/opt/ctsm-container/cime-macros/gnu_container.cmake)
get_cmake_property(VARS_AFTER VARIABLES)
foreach (VAR_AFTER IN LISTS VARS_AFTER)
  if (NOT VAR_AFTER IN_LIST VARS_BEFORE_BUILD_INTERNAL_IGNORE)
    message("CIME_SET_MAKEFILE_VAR ${VAR_AFTER} := ${${VAR_AFTER}}")
  endif()
endforeach()
EOF
out="$(cmake -DCONVERT_TO_MAKE=ON . 2>&1)" || { echo "$out"; fail "macro cmake"; }
echo "$out" | grep "CIME_SET_MAKEFILE_VAR PFUNIT_PATH :=" \
    || { echo "$out" | grep CIME_SET_MAKEFILE_VAR || true
         fail "PFUNIT_PATH was not emitted; find_pfunit() would abort"; }

echo "### and the PFUNIT_PATH env override reaches it too"
rm -rf /tmp/macrocheck2 && cp -r /tmp/macrocheck /tmp/macrocheck2
cd /tmp/macrocheck2 && rm -rf CMakeCache.txt CMakeFiles
out2="$(PFUNIT_PATH=/an/override/path cmake -DCONVERT_TO_MAKE=ON . 2>&1)" \
    || { echo "$out2"; fail "macro cmake (override)"; }
echo "$out2" | grep -q "CIME_SET_MAKEFILE_VAR PFUNIT_PATH := /an/override/path" \
    || fail "\$PFUNIT_PATH env override did not reach PFUNIT_PATH"
echo "override OK"

echo "### a real pFUnit test preprocesses, compiles, links and runs"
mkdir -p /tmp/pfprobe && cd /tmp/pfprobe
cat > test_trivial.pf <<'EOF'
module test_trivial
  use funit
  implicit none
contains

  @test
  subroutine test_arithmetic()
    @assertEqual(4, 2 + 2)
  end subroutine test_arithmetic

end module test_trivial
EOF
# cmake_minimum_required(VERSION 2.8) deliberately MATCHES CTSM's
# src/CMakeLists.txt. Declaring 3.12 here would set policy CMP0074 to NEW,
# make CMake honor pFUnit's GFTL_ROOT hints, and let this probe pass against
# an install layout that CTSM's own build cannot use -- which is exactly the
# false pass this check exists to prevent.
cat > CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 2.8)
project(pfunit_probe LANGUAGES Fortran)
find_package(PFUNIT REQUIRED)
enable_testing()
add_pfunit_ctest(trivial TEST_SOURCES test_trivial.pf)
EOF
# -DCMAKE_PREFIX_PATH=<pfunit_path> and FC=gfortran mirror what run_tests.py
# does with the value it scrapes (it forces mpilib=mpi-serial, so the serial
# compiler is the right one here).
FC=gfortran cmake -DCMAKE_PREFIX_PATH="$pf" . || fail "find_package(PFUNIT)"
make                                          || fail "build pFUnit test"
ctest --output-on-failure                     || fail "run pFUnit test"

echo "ALL PFUNIT SMOKE TESTS PASSED"
INSIDE
