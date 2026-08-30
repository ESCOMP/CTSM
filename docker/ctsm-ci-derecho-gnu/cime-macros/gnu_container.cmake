# CIME cmake macro drop-in for the CTSM CI container. Carries the two
# container-specific settings CIME needs and ccs_config does not provide:
# where pFUnit is, and the gfortran flags ccs_config's own version guard
# cannot apply during a unit-test build (see the two sections below).
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

# ---------------------------------------------------------------------------
# gfortran >= 10 argument-mismatch flags
# ---------------------------------------------------------------------------
# CTSM's mpi-serial code calls mpi_bcast through an implicit interface with
# different actual argument types -- e.g. src/biogeochem/ch4varcon.F90 passes
# a LOGICAL at line 153 and an INTEGER at line 157 -- which gfortran 10+
# rejects outright:
#
#   Error: Type mismatch between actual argument at (1) and actual argument
#          at (2) (INTEGER(4)/LOGICAL(4)).
#
# ccs_config's cmake_macros/gnu.cmake adds -fallow-argument-mismatch for
# exactly this, but guards it on CMAKE_Fortran_COMPILER_VERSION >= 10. That
# guard cannot fire here: CTSM's src/CMakeLists.txt includes
# CIME_initial_setup (which pulls in Macros.cmake, and so this file) at line
# 4, but does not call project() until line 10 -- so when the guard is
# evaluated CMake has not yet probed the compiler and the version variable is
# empty. gnu.cmake's own "Fortran compiler version is" message prints blank in
# the build log, which is the visible symptom.
#
# The container pins gfortran to 12.2.0 in the Dockerfile, so the version test
# has a known answer and the flags can be set unconditionally. This file is
# included last (${COMPILER}_${MACH}.cmake), after gnu.cmake, so appending to
# FFLAGS here reaches CIME_utils.cmake's
#     set(CMAKE_Fortran_FLAGS "${CPPDEFS} ${FFLAGS}")
#
# Not a container-only problem: any machine running CTSM's unit tests with
# gnu hits it. It goes unnoticed upstream because ccs_config defines
# PFUNIT_PATH only for the intel builds on derecho, casper and izumi, so the
# gnu unit-test path is effectively untravelled.
string(APPEND FFLAGS " -fallow-argument-mismatch -fallow-invalid-boz")
