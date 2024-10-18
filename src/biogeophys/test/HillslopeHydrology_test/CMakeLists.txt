set (pfunit_sources
  test_hillslopehydrologyUtils.pf)

add_pfunit_ctest(HillslopeHydrologyUtils
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
