set (pfunit_sources
  test_glc_elevclasses_are_same.pf)

add_pfunit_ctest(initInterpUtils
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share)
