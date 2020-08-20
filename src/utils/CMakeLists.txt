# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

set(genf90_files
  array_utils.F90.in
  )

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_SOURCE_DIR} clm_genf90_sources)

sourcelist_to_parent(clm_genf90_sources)

list(APPEND clm_sources "${clm_genf90_sources}")

list(APPEND clm_sources
  ctsm_AnnualFluxDribbler.F90
  ctsm_Domain.F90
  ctsm_NlUtils.F90
  ctsm_TimeManager.F90
  ctsm_QuadraticSolver.F90
  ctsm_FileUtils.F90
  ctsm_IssueFixedMetadataHandler.F90
  ctsm_Numerics.F90
  )

sourcelist_to_parent(clm_sources)
