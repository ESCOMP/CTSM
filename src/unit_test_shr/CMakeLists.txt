set(genf90_files
  unittestArrayMod.F90.in
  )

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_SOURCE_DIR} clm_genf90_sources)

sourcelist_to_parent(clm_genf90_sources)

list(APPEND clm_sources "${clm_genf90_sources}")

list(APPEND clm_sources
  ctsm_UnitTestFilterBuilder.F90
  ctsm_UnitTestGlcMec.F90
  ctsm_UnitTestSimpleSubgridSetups.F90
  ctsm_UnitTestSubgrid.F90
  ctsm_UnitTestTimeManager.F90
  ctsm_UnitTestUtils.F90
  ctsm_UnitTestctsm_WaterTypeFactory.F90
  )

sourcelist_to_parent(clm_sources)
