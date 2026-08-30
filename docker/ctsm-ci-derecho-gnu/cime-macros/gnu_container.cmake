# CIME cmake macro drop-in for the CTSM CI container: tell CIME's Fortran
# unit-test driver where pFUnit is.
#
# WHY THIS FILE EXISTS
# cime/scripts/fortran_unit_testing/run_tests.py locates pFUnit in
# find_pfunit(), which scrapes PFUNIT_PATH out of the CMake macros. It reads
# neither the environment nor --cmake-args:
# ccs_config/machines/cmake_macros/Macros.cmake snapshots the CMake variable
# set BEFORE including the macro files and emits only variables that appeared
# afterwards, so anything passed as -DPFUNIT_PATH=... is filtered right back
# out. PFUNIT_PATH has to be set by a macro *file*.
#
# HOW IT GETS PICKED UP
# CIME's copy_local_macros_to_dir() (cime/CIME/utils.py) copies
# $HOME/.cime/*.cmake into the case's cmake_macros/ directory. Macros.cmake
# includes ${COMPILER}_${MACH}.cmake last, so gnu_container.cmake takes
# precedence over everything ccs_config ships -- and ccs_config has no
# gnu_container.cmake, so nothing is being shadowed.
#
# This keeps the container self-contained: no ccs_config or cime edit needed.
# The eventual upstream fix is a one-line set(PFUNIT_PATH ...) in
# ccs_config/machines/container/container.cmake.
#
# NOTE ON $HOME: GitHub Actions container jobs override HOME (to /github/home),
# so the copy baked into /root/.cime is not found in CI. The workflow step has
# to place this file itself:
#   mkdir -p "$HOME/.cime"
#   cp /opt/ctsm-container/cime-macros/gnu_container.cmake "$HOME/.cime/"
#
# WHY THE SHARED PREFIX AND NOT THE PFUNIT-4.8 SUBDIRECTORY
# run_tests.py turns this value into -DCMAKE_PREFIX_PATH for the test build.
# pFUnit installs four sibling packages under one prefix -- PFUNIT-4.8,
# GFTL-1.11, GFTL_SHARED-1.7, FARGPARSE-1.6 -- and PFUNITConfig.cmake locates
# its own three dependencies ONLY by set()-ing GFTL_ROOT / GFTL_SHARED_ROOT /
# FARGPARSE_ROOT before find_dependency(). Those <pkg>_ROOT variables are
# honored only under policy CMP0074, and CTSM's src/CMakeLists.txt opens with
# cmake_minimum_required(VERSION 2.8), which leaves CMP0074 unset -- so CMake
# ignores them and reports "PFUNIT could not be found because dependency GFTL
# could not be found", after which add_pfunit_ctest is never defined.
#
# Naming the shared prefix sidesteps the policy entirely: find_package locates
# all four as siblings. Do not "tighten" this to the PFUNIT-4.8 subdirectory.
#
# The env-var branch lets an alternate pFUnit build be tested without
# rebuilding the image. An environment variable is not a CMake variable, so
# set()-ing it here is precisely what makes it visible to find_pfunit().
if (DEFINED ENV{PFUNIT_PATH})
  set(PFUNIT_PATH "$ENV{PFUNIT_PATH}")
else()
  set(PFUNIT_PATH "/usr/local/pfunit-4.8.0")
endif()
