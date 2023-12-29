set(pfunit_sources
  test_dyn_cons_biogeophys.pf)

add_pfunit_ctest(dynConsBiogeophys
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf_wrf_timemgr)
