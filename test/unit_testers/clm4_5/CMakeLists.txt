cmake_minimum_required(VERSION 2.8)
project(clm45_tests Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CESM_CMAKE_MODULE_DIRECTORY})
include(CESM_utils)

set(CLM_ROOT "../../..")
set(CESM_ROOT "${CLM_ROOT}/../../..")

# Add source directories from other models (csm_share, etc.). This should be
# done first, so that in case of name collisions, the CLM versions take
# precedence (when there are two files with the same name, the one added later
# wins).
add_subdirectory(${CESM_ROOT}/models/csm_share/shr csm_share)

# Add CLM source directories
add_subdirectory(${CLM_ROOT}/src/clm4_5/biogeophys clm_biogeophys)
add_subdirectory(${CLM_ROOT}/src/clm4_5/dyn_subgrid clm_dyn_subgrid)
add_subdirectory(${CLM_ROOT}/src/clm4_5/main clm_main)

# Add general unit test directories (mocked out files, etc.)
add_subdirectory(mock)

# Add each main unit test directory
add_subdirectory(dyn_subgrid)
add_subdirectory(biogeophys)
