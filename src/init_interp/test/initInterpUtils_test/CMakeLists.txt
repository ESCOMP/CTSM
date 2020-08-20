set (pfunit_sources
  test_glc_elevclasses_are_same.pf)

create_pFUnit_test(ctsm_InitInterpUtils test_ctsm_InitInterpUtils_exe
  "${pfunit_sources}" "")

target_link_libraries(test_ctsm_InitInterpUtils_exe clm csm_share)