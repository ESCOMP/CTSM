cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
include(CIME_initial_setup)

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../tools/mksurfdata_esmf/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../share/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../components/cmeps/cmake")

project(clm_tests Fortran C)

include(CIME_utils)

set(CLM_ROOT "..")

# find needed external packages
# NetCDF is required -- because PIO and NetCDF are required by the standard default ESMF libraries
find_package(NetCDF 4.7.4 REQUIRED Fortran)

# The following - for finding PIO - is copied from the share CMakeLists.txt
if (DEFINED ENV{PIO_ROOT})
  message("PIO_ROOT is $ENV{PIO_ROOT}")
else()
  if (DEFINED PIO)
    set(PIO_PATH ${PIO})
  else()
    set(PIO_PATH $ENV{PIO})
  endif()
  find_package(PIO REQUIRED COMPONENT C Fortran PATH ${PIO_PATH})
endif()
# The following line differs a bit from what's in the share CMakeLists.txt
include_directories(${PIO_Fortran_INCLUDE_DIR})
link_libraries(${PIO_Fortran_LIBRARIES})

# The following - for finding ESMF - is copied from the share CMakeLists.txt
if (DEFINED ENV{ESMF_ROOT})
  list(APPEND CMAKE_MODULE_PATH $ENV{ESMF_ROOT}/cmake)
endif()
find_package(ESMF REQUIRED)
# This adds include directories needed for ESMF
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_F90COMPILEPATHS} ")
# The following line differs a bit from what's in the share CMakeLists.txt
link_directories(${ESMF_LIBSDIR})

# Add source directories from other share code (csm_share, etc.). This should be
# done first, so that in case of name collisions, the CLM versions take
# precedence (when there are two files with the same name, the one added later
# wins).
add_subdirectory(${CLM_ROOT}/share/src csm_share)
add_subdirectory(${CLM_ROOT}/share/unit_test_stubs/util csm_share_stubs)

# Add files needed from CMEPS
list ( APPEND drv_sources_needed
  ${CLM_ROOT}/components/cmeps/cesm/nuopc_cap_share/glc_elevclass_mod.F90
  ${CLM_ROOT}/components/cmeps/cesm/nuopc_cap_share/shr_dust_emis_mod.F90
  )

# Add CLM source directories
add_subdirectory(${CLM_ROOT}/src/utils clm_utils)
add_subdirectory(${CLM_ROOT}/src/biogeochem clm_biogeochem)
add_subdirectory(${CLM_ROOT}/src/soilbiogeochem clm_soilbiogeochem)
add_subdirectory(${CLM_ROOT}/src/biogeophys clm_biogeophys)
add_subdirectory(${CLM_ROOT}/src/dyn_subgrid clm_dyn_subgrid)
add_subdirectory(${CLM_ROOT}/src/main clm_main)
add_subdirectory(${CLM_ROOT}/src/init_interp clm_init_interp)
add_subdirectory(${CLM_ROOT}/src/self_tests clm_self_tests)

# Add general unit test directories (stubbed out files, etc.)
add_subdirectory(unit_test_stubs)
add_subdirectory(unit_test_shr)

# Remove shr_mpi_mod from share_sources.
# This is needed because we want to use the mock shr_mpi_mod in place of the real one
#
# TODO: this should be moved into a general-purpose function in Sourcelist_utils.
# Then this block of code could be replaced with a single call, like:
# remove_source_file(${share_sources} "shr_mpi_mod.F90")
foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_mpi_mod.F90" match_found ${sourcefile})
  if(match_found)
    list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()

# Build libraries containing stuff needed for the unit tests.
# Eventually, these add_library calls should probably be distributed into the correct location, rather than being in this top-level CMakeLists.txt file.
add_library(csm_share ${share_sources} ${drv_sources_needed})
declare_generated_dependencies(csm_share "${share_genf90_sources}")
add_library(clm ${clm_sources})
declare_generated_dependencies(clm "${clm_genf90_sources}")
add_dependencies(clm csm_share esmf)

# We need to look for header files here, in order to pick up shr_assert.h
include_directories(${CLM_ROOT}/share/include)

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories (${NETCDF}/include)

# Directories and libraries to include in the link step
link_directories(${CMAKE_CURRENT_BINARY_DIR})
link_libraries( netcdf )

# Add the test directories
# Note: it's possible that these could be added by each source directory that
# has tests in it. However, it appears that the order needs to be done
# carefully: for example, include_directories and link_directories needs to be
# done before adding the tests themselves.
add_subdirectory(${CLM_ROOT}/src/unit_test_shr/test clm_unit_test_shr_test)
add_subdirectory(${CLM_ROOT}/src/utils/test clm_utils_test)
add_subdirectory(${CLM_ROOT}/src/biogeophys/test clm_biogeophys_test)
add_subdirectory(${CLM_ROOT}/src/biogeochem/test clm_biogeochem_test)
add_subdirectory(${CLM_ROOT}/src/soilbiogeochem/test clm_soilbiogeochem_test)
add_subdirectory(${CLM_ROOT}/src/dyn_subgrid/test clm_dyn_subgrid_test)
add_subdirectory(${CLM_ROOT}/src/main/test clm_main_test)
add_subdirectory(${CLM_ROOT}/src/init_interp/test clm_init_interp_test)
add_subdirectory(${CLM_ROOT}/src/self_tests/test clm_self_tests_test)

# Add driver unit test directories
# (these should be moved to the appropriate submodule)
add_subdirectory(${CLM_ROOT}/src/drv_test drv_test)

