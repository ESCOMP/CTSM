# This CMakeLists.txt file is currently used just for building unit tests.

cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
include(CIME_initial_setup)

project(mksurfdat_tests Fortran)

include(CIME_utils)

# Build library containing stuff needed for the unit tests
list(APPEND mksurfdat_sources
  shr_kind_mod.F90
  shr_log_mod.F90
  nanMod.F90
  shr_string_mod.F90
  fileutils.F90
  shr_timer_mod.F90
  shr_file_mod.F90
  mkgridmapMod.F90
  mkindexmapMod.F90
  mkpftConstantsMod.F90
  mkpctPftTypeMod.F90
  mkpftMod.F90
  mkdomainMod.F90
  mkgridmapMod.F90
  mkutilsMod.F90
  mkpftUtilsMod.F90
  mksoilUtilsMod.F90
  mkvarctl.F90
  mkvarpar.F90
  shr_const_mod.F90
  shr_sys_mod.F90
  unit_test_stubs/abort.F90
  unit_test_stubs/mkncdio.F90)
add_library(mksurfdat ${mksurfdat_sources})

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add the test directory
add_subdirectory(test)
