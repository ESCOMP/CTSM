set (pfunit_sources
  test_set_mindist.pf
  test_set_single_match.pf
  initInterpMindistTestUtils.pf)

add_pfunit_ctest(initInterpMindist
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share)
