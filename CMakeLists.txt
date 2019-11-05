cmake_minimum_required (VERSION 3.5.2)
project (ctsm C Fortran)

# Local CMake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/CMakeModules)
find_package(ESMF REQUIRED)

# cime build
set(CIME_ROOT "/glade/work/jhamman/lilac/cime")  # TDOO: generalize this
file(GLOB CIME_SOURCES ${CIME_ROOT}/src/share/streams/*F90
                       ${CIME_ROOT}/src/share/util/*F90
                       ${CIME_ROOT}/src/share/esmf_wrf_timemgr/*F90
                       ${CIME_ROOT}/src/share/streams/*.inc
                       ${CIME_ROOT}/src/share/util/*.inc
                       ${CIME_ROOT}/src/share/esmf_wrf_timemgr/*.inc
                       ${CIME_ROOT}/src/share/streams/*.h
                       ${CIME_ROOT}/src/share/util/*.h
                       ${CIME_ROOT}/src/share/esmf_wrf_timemgr/*.h)
add_library(cime SHARED ${CIME_SOURCES})
target_link_libraries(cime ESMF)

# ctsm build
file(GLOB_RECURSE CTSM_SOURCES src/**/*.F90 src/**/*.h src/**/*.inc)
add_library(ctsm SHARED ${CTSM_SOURCES})
target_link_libraries(ctsm cime ESMF)
