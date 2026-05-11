add_pfunit_ctest(dynTimeInfo
  TEST_SOURCES "test_dynTimeInfo.pf"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
