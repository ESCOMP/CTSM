set(genf90_files
  unittestArrayMod.F90.in
  )

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_BINARY_DIR} clm_genf90_sources)

sourcelist_to_parent(clm_genf90_sources)

list(APPEND clm_sources "${clm_genf90_sources}")

list(APPEND clm_sources
  unittestDustEmisInputs.F90
  unittestFilterBuilderMod.F90
  unittestGlcMec.F90
  unittestInitializeAndFinalize.F90
  unittestSimpleSubgridSetupsMod.F90
  unittestSubgridMod.F90
  unittestTimeManagerMod.F90
  unittestUtils.F90
  unittestWaterTypeFactory.F90
  )

sourcelist_to_parent(clm_sources)
