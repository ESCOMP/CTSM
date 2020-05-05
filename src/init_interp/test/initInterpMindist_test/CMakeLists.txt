set (pfunit_sources
  test_set_mindist.pf
  test_set_single_match.pf
  initInterpMindistTestUtils.pf)

create_pFUnit_test(initInterpMindist test_initInterpMindist_exe
  "${pfunit_sources}" "")

target_link_libraries(test_initInterpMindist_exe clm csm_share)