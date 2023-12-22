cmake_minimum_required(VERSION 2.8.12.1)

##include("/glade/work/negins/UFSCOMP/cime/tools/Macros.cmake")

set     (CIME_ROOT "/glade/work/negins/UFSCOMP/cime")
message ("CIME_ROOT:                   ${CIME_ROOT}")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/CMakeModules")
message ("CMAKE_MODULE_PATH:           ${CMAKE_MODULE_PATH}")


set     (CIME_CMAKE_MODULE_DIRECTORY "/glade/work/negins/UFSCOMP/cime/src/CMake/")
message ("CIME_CMAKE_MODULE_DIRECTORY: ${CIME_CMAKE_MODULE_DIRECTORY}")


list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
message ("CMAKE_MODULE_PATH:           ${CMAKE_MODULE_PATH}")

set     (MACRO_ROOT  "/glade/work/negins/UFSCOMP/cime/tools/")
include(${MACRO_ROOT}/Macros.cmake)

list(APPEND CMAKE_MODULE_PATH ${MACRO_ROOT})
message ("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")


set     (CLM_ROOT "/glade/work/negins/UFSCOMP/components/clm")

message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

include_directories (${CMAKE_SOURCE_DIR}/cmake/CMakeModules/)
include (${CMAKE_SOURCE_DIR}/cmake/CMakeModules/genf90_utils.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/CMakeModules/Sourcelist_utils.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/CMakeModules/pFUnit_utils.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/CMakeModules/FindpFUnit.cmake)


#include (Macros.cmake)
#include(CIME_initial_setup)

message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

### -------------------------------------------------------------

# project name
project(LILAC Fortran C)
enable_language(Fortran)


# This definition is needed to avoid having ESMF depend on mpi
add_definitions(-DHIDE_MPI)


message("----------------------------------------------------")
message ("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
message ("CMAKE_MODULE_PATH:        ${CMAKE_MODULE_PATH}")
message("----------------------------------------------------")



message("----------------------------------------------------")
# Add source directories from other share code (csm_share, etc.). This should be
# done first, so that in case of name collisions, the CLM versions take
# precedence (when there are two files with the same name, the one added later
# wins).
add_subdirectory(${CIME_ROOT}/../share/src csm_share)
add_subdirectory(${CIME_ROOT}/../share/unit_test_stubs/util csm_share_stubs)
add_subdirectory(${CIME_ROOT}/../share/esmf_wrf_timemgr esmf_wrf_timemgr)
add_subdirectory(${CIME_ROOT}/../components/cpl7/driver drv_share)
message("----------------------------------------------------")

# Extract just the files we need from drv_share
set (drv_sources_needed_base
      glc_elevclass_mod.F90
        )
extract_sources("${drv_sources_needed_base}" "${drv_sources}" drv_sources_needed)

message("~~~~~~~~~~~~~~~~~~~~~~CLM_ROOT~~~~~~~~~~~~~~~~~~~~~~")
# Add CLM source directories (these add their own test directories)
add_subdirectory(${CLM_ROOT}/src/utils clm_utils)
add_subdirectory(${CLM_ROOT}/src/biogeochem clm_biogeochem)
add_subdirectory(${CLM_ROOT}/src/soilbiogeochem clm_soilbiogeochem)
add_subdirectory(${CLM_ROOT}/src/biogeophys clm_biogeophys)
add_subdirectory(${CLM_ROOT}/src/dyn_subgrid clm_dyn_subgrid)
add_subdirectory(${CLM_ROOT}/src/main clm_main)
add_subdirectory(${CLM_ROOT}/src/init_interp clm_init_interp)
add_subdirectory(${CLM_ROOT}/src/fates/main fates_main)

# Add general unit test directories (stubbed out files, etc.)
add_subdirectory(unit_test_stubs)
add_subdirectory(unit_test_shr)


# Remove shr_mpi_mod from share_sources.
# This is needed because we want to use the mock shr_mpi_mod in place of the real one
#
# TODO: this should be moved into a general-purpose function in Sourcelist_utils.
# Then this block of code could be replaced with a single call, like:
# remove_source_file(${share_sources} "shr_mpi_mod.F90")}

foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_mpi_mod.F90" match_found ${sourcefile})
  if(match_found)
        list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()


# We rely on pio for cmake utilities like findnetcdf.cmake, so that we don't
# need to duplicate this cmake code
message ("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
list    (APPEND CMAKE_MODULE_PATH "${CIME_ROOT}/../libraries/parallelio/cmake")
message ("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")


add_subdirectory (${CIME_ROOT}/..libraries/parallelio/tests)

message("----------------------------------------------------")
option(ENABLE_PFUNIT "Enable pfUnit testing Framework" ON)
if (ENABLE_PFUNIT)
    find_package(pfUnit)
    include(pfUnit_utils)
    include_directories("${PFUNIT_INCLUDE_DIRS}")
endif (ENABLE_PFUNIT)
message("----------------------------------------------------")


find_package(MPI REQUIRED)
# TODO: This should be found from the find_package call but its not working
#set(CMAKE_Fortran_COMPILER "/usr/lib64/mpich/bin/mpif90")
find_package(ESMF REQUIRED)


message("------------include (CIME_utils)--------------------")
include(CIME_utils)
message("----------------------------------------------------")

find_package(NetCDF COMPONENTS C Fortran)
include_directories(${NetCDF_C_INCLUDE_DIRS} ${NetCDF_Fortran_INCLUDE_DIRS})
message("NetCDF_C_INCLUDE_DIRS: ${NetCDF_C_INCLUDE_DIRS}")
message("----------------------------------------------------")

##=======##
#set(CESM_ROOT "/glade/work/negins/UFSCOMP/")
#set(CSM_SHR "/glade/work/negins/UFSCOMP/components/clm/src/unit_test_stubs/csm_share/")

#add_subdirectory(${CESM_ROOT}/models/csm_share/shr csm_share)
#add_subdirectory(${CSM_SHR} )

message("----------------------------------------------------")



# -lclm libclm.a
SET(NAMES libclm.a)

#find_library(LIB_TO_INCLUDE 
#    libclm.a
#    PATHS /glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/lib/)
#find_library(LIB_TO_INCLUDE  /glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/lib/)

#message(STATUS "include_directories for ${NAMES}: ${LIB_TO_INCLUDE}")
#include_directories(${LIB_TO_INCLUDE})
#link_directories(${LIB_TO_INCLUDE})
#message(STATUS "include_directories for ${NAMES}: ${LIB_TO_INCLUDE}")
#find_library(LIB_TO_INCLUDE /glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/lib/)
#message(STATUS "include_directories: ${LIB_TO_INCLUDE}")
#target_link_libraries (${LIB_TO_INCLUDE})


# Local CMake modules

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_COMPILER_LINE}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_LINK_LINE} -g -cpp")


message(STATUS "==============================================================")
message(STATUS "Fortran    Compiler : ${CMAKE_Fortran_COMPILER}")
message(STATUS "cmake Fortran Flags : ${CMAKE_Fortran_FLAGS}")
message(STATUS "==============================================================")
message(STATUS "==============================================================")


#add_executable("lilac.exe" ../lilac/*.F90)

#
# Compile.
#

file(GLOB_RECURSE SOURCES  lilac/*.F90)
#add_subdirectory(lilac)
#add_executable(${PROJECT_NAME}.exe ../lilac/demo_driver.F90
#                ../lilac/lilac_mod.F90 ../lilac/atmos_cap.F90  ../lilac/lilac_utils.F90 
#                ../lilac/lnd_cap.F90 ../lilac/cpl_mod.F90)

add_executable (${PROJECT_NAME}.exe ${SOURCES})
target_link_libraries(${PROJECT_NAME}.exe ${LIB_TO_INCLUDE})

#add_subdirectory(lilac)
#add_subdirectory(tests)
