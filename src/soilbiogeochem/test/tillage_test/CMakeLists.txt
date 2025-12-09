add_pfunit_ctest(tillage
  TEST_SOURCES "test_tillage.pf"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
