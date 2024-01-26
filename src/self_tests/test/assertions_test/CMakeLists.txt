set (pfunit_sources
  test_assertions.pf)

add_pfunit_ctest(assertions
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
