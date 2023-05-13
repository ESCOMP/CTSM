add_pfunit_ctest(clm_time_manager
  TEST_SOURCES "test_clm_time_manager.pf"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
