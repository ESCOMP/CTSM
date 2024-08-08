set (pfunit_sources
  test_CNPhenology.pf
  test_CNPhenology_gregorian.pf)

add_pfunit_ctest(CNPhenology
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
