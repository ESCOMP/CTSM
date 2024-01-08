set (pfunit_sources
  test_ComputeSeedAmounts.pf)

add_pfunit_ctest(CNVegComputeSeed
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
