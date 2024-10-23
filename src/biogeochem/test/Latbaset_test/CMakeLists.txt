set (pfunit_sources
  test_Latbaset.pf)

add_pfunit_ctest(CropTypeLatbaset
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
