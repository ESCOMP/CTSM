# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

set(genf90_files
  Assertions.F90.in
  )

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_BINARY_DIR} clm_genf90_sources)

sourcelist_to_parent(clm_genf90_sources)

list(APPEND clm_sources "${clm_genf90_sources}")

sourcelist_to_parent(clm_sources)
