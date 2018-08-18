cmake_minimum_required(VERSION 2.8.12.1)

project(LILAC Fortran)
enable_language(Fortran)

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

if(DEFINED PFUNIT_INSTALL)
    message(STATUS "Manual setup of variable PFUNIT_INSTALL: ${PFUNIT_INSTALL}")
    set(PFUNIT_DIR ${PFUNIT_INSTALL})
else()
    include(ExternalProject)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pfunit
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        )
    ExternalProject_Add(pfunit
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/pfunit
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/pfunit-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_BINARY_DIR}/external/pfunit/mod)
    set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pfunit)
endif()


if(DEFINED ESMF_INSTALL)
    message(STATUS "Manual setup of variable ESMF_INSTALL: ${ESMF_INSTALL}")
    set(ESMF_DIR ${ESMF_INSTALL})

else()
    message(STATUS "We will build ESMF")
    include(ExternalProject)

    # set(ENV{ESMFMKFILE} "${PROJECT_BINARY_DIR}/external/esmf/esmf.mk")

    set(ExternalProjectCMakeArgs
        ESMF_DIR=${ESMF_DIR}
        ESMF_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        ESMF_INSTALL_BINDIR=${PROJECT_BINARY_DIR}/bin/
        ESMF_INSTALL_DOCDIR=${PROJECT_BINARY_DIR}/doc/
        ESMF_INSTALL_HEADERDIR=${PROJECT_BINARY_DIR}/include/
        ESMF_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/
        ESMF_INSTALL_MODDIR=${PROJECT_BINARY_DIR}/mod/
        )

    ExternalProject_Add(esmf
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/esmf
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        BUILD_IN_SOURCE TRUE
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make
        INSTALL_COMMAND pwd && make install
        )

    # link to ESMF
    # find_package(ESMF)
    # include_directories(${PROJECT_BINARY_DIR}/include/)
    # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\ ${ESMF_COMPILER_LINE}")
    # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_LINK_LINE} -g -cpp")
    # TODO figure out how to get these out once ESMF is compiled
    set(ESMF_COMPILER_LINE, "-O   -m64 -mcmodel=small -ffree-line-length-none  -fopenmp -I/Users/jhamman/Dropbox/src/lilac/external/esmf/mod/modO/Darwin.gfortran.64.mpiuni.default -I/Users/jhamman/Dropbox/src/lilac/external/esmf/src/include -ffree-form")
    set(ESMF_LINK_LINE, "-m64 -mcmodel=small  -fopenmp -L/Users/jhamman/Dropbox/src/lilac/external/esmf/lib/libO/Darwin.gfortran.64.mpiuni.default  -L/opt/local/lib/gcc5/gcc/x86_64-apple-darwin17/5.5.0/../../../ -lesmf  -lstdc++")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_COMPILER_LINE}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_LINK_LINE} -g -cpp")



endif()

add_subdirectory(lilac)
add_subdirectory(tests)
