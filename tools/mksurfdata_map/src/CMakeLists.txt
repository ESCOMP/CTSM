# This CMakeLists.txt file is currently used just for building unit tests.

cmake_minimum_required(VERSION 2.8)
project(mksurfdat_tests Fortran)

list(APPEND CMAKE_MODULE_PATH ${CESM_CMAKE_MODULE_DIRECTORY})
include(CESM_utils)

# Build library containing stuff needed for the unit tests
list(APPEND mksurfdat_sources shr_kind_mod.F90 mkpftConstantsMod.F90 mkpctPftTypeMod.F90 mkutilsMod.F90 mkpftUtilsMod.F90)
add_library(mksurfdat ${mksurfdat_sources})

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add the test directory
add_subdirectory(test)