cmake_minimum_required(VERSION 2.8.12.1)

project(LILAC Fortran)
enable_language(Fortran)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/CMakeModules")

find_package(MPI REQUIRED)
# TODO: This should be found from the find_package call but its not working
set(CMAKE_Fortran_COMPILER "/usr/lib64/mpich/bin/mpif90")
find_package(ESMF REQUIRED)

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
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I /usr/include/ -I/usr/src/lilac/external/esmf/src/Infrastructure/Util/include -I/usr/src/lilac/external/esmf/build_config/Linux.gfortran.default -I /usr/src/lilac/external/esmf/src/include")
message("CMAKE_Fortran_FLAGS:" ${CMAKE_Fortran_FLAGS})

add_subdirectory(lilac)
add_subdirectory(tests)
