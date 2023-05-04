set (pfunit_sources
  test_water_type.pf)

add_pfunit_ctest(water_type
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
