set (pfunit_sources
  test_water_tracer_container.pf)

add_pfunit_ctest(water_tracer_container
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
