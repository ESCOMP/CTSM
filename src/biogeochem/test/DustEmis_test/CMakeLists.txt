set (pfunit_sources
  test_DustEmis.pf)

add_pfunit_ctest(DustEmis
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
