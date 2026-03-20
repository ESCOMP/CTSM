add_pfunit_ctest(dynInitColumns
  TEST_SOURCES "test_init_columns.pf"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
