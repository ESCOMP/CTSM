file(GLOB_RECURSE SOURCES *.f90 *.h)
add_executable("test_hello_world" ${SOURCES} )
target_link_libraries(lilac)
