add_pfunit_ctest(annual_flux_dribbler
  TEST_SOURCES "test_annual_flux_dribbler.pf"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
