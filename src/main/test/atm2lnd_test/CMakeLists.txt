set(pfunit_sources
  test_downscale_forcings.pf
  test_partition_precip.pf)

add_pfunit_ctest(atm2lnd
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
