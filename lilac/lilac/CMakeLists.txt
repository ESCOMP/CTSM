# Compile LILAC library
file(GLOB_RECURSE LILAC_SOURCES *.f90 *.h)
add_library(lilac ${LILAC_SOURCES})
target_include_directories(lilac PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lilac)
