cmake_minimum_required(VERSION 2.8.12.1)



### ********** ###
set     (CIME_ROOT "/glade/work/negins/UFSCOMP/cime")
message ("CIME_ROOT:                   ${CIME_ROOT}")
set     (CIME_CMAKE_MODULE_DIRECTORY "/glade/work/negins/UFSCOMP/cime/src/CMake/")
message ("CIME_CMAKE_MODULE_DIRECTORY: ${CIME_CMAKE_MODULE_DIRECTORY}")

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
message ("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
include(CIME_initial_setup)

message("----------------------------------------------------")

project(LILAC Fortran C)
enable_language(Fortran)

# We rely on pio for cmake utilities like findnetcdf.cmake, so that we don't
# need to duplicate this cmake code
message ("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
list(APPEND CMAKE_MODULE_PATH "${CIME_ROOT}/src/externals/pio2/cmake")
message ("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

message("----------------------------------------------------")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/CMakeModules")
message ("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
message ("CMAKE_MODULE_PATH:        ${CMAKE_MODULE_PATH}")
message("----------------------------------------------------")

find_package(MPI REQUIRED)
# TODO: This should be found from the find_package call but its not working
#set(CMAKE_Fortran_COMPILER "/usr/lib64/mpich/bin/mpif90")
find_package(ESMF REQUIRED)


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

#==============###
set(SHARE_ROOT "${CIME_ROOT}/src/share")
#add_subdirectory(${CIME_ROOT}/src/share/util csm_share)
#add_subdirectory(${CIME_ROOT}/src/share/unit_test_stubs/util csm_share_stubs)
#add_subdirectory(${CIME_ROOT}/src/share/esmf_wrf_timemgr esmf_wrf_timemgr)
#add_subdirectory(${CIME_ROOT}/src/drivers/mct/shr drv_share)
#message("----------------------------------------------------")


add_subdirectory(${SHARE_ROOT}/util csm_share)
add_subdirectory(${SHARE_ROOT}/unit_test_stubs/util csm_share_stubs)
include_directories(${SHARE_ROOT}/include)

# esmf_wrf_timemgr not built here because it depends on csm_share.
#add_subdirectory(${SHARE_ROOT}/esmf_wrf_timemgr esmf_wrf_timemgr)
#include_directories(${SHARE_ROOT}/esmf_wrf_timemgr)




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
# TODO: This should not be necessary but certain header files are missing from the build
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I /usr/include/ -I/usr/src/lilac/external/esmf/src/Infrastructure/Util/include -I/usr/src/lilac/external/esmf/build_config/Linux.gfortran.default -I /usr/src/lilac/external/esmf/src/include")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/usr/include/ -I/usr/src/esmf/src/Infrastructure/Util/include/ -I/usr/src/esmf/build_config/Linux.gfortran.default -I /usr/src/esmf/src/include")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/include -L/glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/lib -lclm -L/glade/scratch/negins/test_clean/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/c1a1l1/lib -lcsm_share -L/glade/scratch/negins/test_clean/bld/intel/mpt/nodebug/nothreads/nuopc/lib -lpiof -lpioc -lgptl -lmct -lmpeu -mkl=cluster -L/glade/u/apps/ch/opt/pnetcdf/1.11.0/mpt/2.19/intel/19.0.2/lib -lpnetcdf -L/glade/u/home/dunlap/ESMF-INSTALL/intel19/8.0.0bs32/lib/libO/Linux.intel.64.mpt.default -L/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib -L/glade/u/home/dunlap/YAML-INSTALL/lib -Wl,-rpath,/glade/u/home/dunlap/ESMF-INSTALL/intel19/8.0.0bs32/lib/libO/Linux.intel.64.mpt.default -Wl,-rpath,/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib -Wl,-rpath,/glade/u/home/dunlap/YAML-INSTALL/lib -lesmf -cxxlib -lrt -ldl -lnetcdff -lnetcdf -lyaml-cpp -cxxlib ")


message(STATUS "==============================================================")
message(STATUS "Fortran    Compiler : ${CMAKE_Fortran_COMPILER}")
message(STATUS "cmake Fortran Flags : ${CMAKE_Fortran_FLAGS}")
message(STATUS "==============================================================")
message(STATUS "==============================================================")


#add_executable("lilac.exe" ../lilac/*.F90)

#
# Compile.
#
file(GLOB_RECURSE sources  lilac/*.F90)

#add_executable(${PROJECT_NAME}.exe ../lilac/demo_driver.F90
#                ../lilac/lilac_mod.F90 ../lilac/atmos_cap.F90  ../lilac/lilac_utils.F90 
#                ../lilac/lnd_cap.F90 ../lilac/cpl_mod.F90)

add_executable(${PROJECT_NAME}.exe ${sources})
target_link_libraries(${PROJECT_NAME}.exe ${LIB_TO_INCLUDE})
#emo_driver.o atmos_cap.o lilac_mod.o lilac_utils.o cpl_mod.o lnd_cap.o
#add_subdirectory(lilac)
#add_subdirectory(tests)
