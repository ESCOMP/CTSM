cmake_minimum_required(VERSION 2.8.12.1)

project(LILAC Fortran)
enable_language(Fortran)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/CMakeModules")

find_package(MPI REQUIRED)
# TODO: This should be found from the find_package call but its not working
#set(CMAKE_Fortran_COMPILER "/usr/lib64/mpich/bin/mpif90")
find_package(ESMF REQUIRED)

# -lclm libclm.a
SET(NAMES libclm.a)

find_library(LIB_TO_INCLUDE 
    libclm.a
    PATHS /glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/lib/)
#find_library(LIB_TO_INCLUDE  /glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/lib/)

message(STATUS "include_directories for ${NAMES}: ${LIB_TO_INCLUDE}")
#include_directories(${LIB_TO_INCLUDE})
#link_directories(${LIB_TO_INCLUDE})
message(STATUS "include_directories for ${NAMES}: ${LIB_TO_INCLUDE}")
#find_library(LIB_TO_INCLUDE /glade/scratch/negins/baghale6/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/lib/)
message(STATUS "include_directories: ${LIB_TO_INCLUDE}")
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
