add_pfunit_ctest(initVertical
  TEST_SOURCES "test_initVertical.pf"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
