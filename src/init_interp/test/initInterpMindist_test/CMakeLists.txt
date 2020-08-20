set (pfunit_sources
  test_set_mindist.pf
  test_set_single_match.pf
  ctsm_InitInterpMindistTestUtils.pf)

create_pFUnit_test(ctsm_InitInterpMindist test_ctsm_InitInterpMindist_exe
  "${pfunit_sources}" "")

target_link_libraries(test_ctsm_InitInterpMindist_exe clm csm_share)