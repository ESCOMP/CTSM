set(pfunit_sources
  test_accumul.pf)

add_pfunit_ctest(accumul
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
