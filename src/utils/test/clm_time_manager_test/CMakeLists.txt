add_pfunit_ctest(clm_time_manager
  TEST_SOURCES "test_clm_time_manager.pf"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
