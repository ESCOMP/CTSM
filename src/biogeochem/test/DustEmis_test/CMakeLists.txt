set (pfunit_sources
  test_DustEmisBase.pf
  test_DustEmisZender2003.pf)

add_pfunit_ctest(DustEmis
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
