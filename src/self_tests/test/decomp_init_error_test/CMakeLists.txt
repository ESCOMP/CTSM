set (pfunit_sources
  test_decomp_init_errors.pf)

add_pfunit_ctest(test_decomp_init_error
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf fates
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
