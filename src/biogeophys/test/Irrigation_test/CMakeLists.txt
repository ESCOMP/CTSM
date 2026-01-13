set (pfunit_sources
  test_irrigation.pf)

add_pfunit_ctest(irrigation
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
