set (pfunit_sources
  test_CIsoAtmTimeSeries.pf)

add_pfunit_ctest(CIsoAtmTimeSeries
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share
  # I don't think we need esmf here
  #LINK_LIBRARIES clm csm_share esmf
  #EXTRA_FINALIZE unittest_finalize_esmf
  #EXTRA_USE unittestInitializeAndFinalize)
