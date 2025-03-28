set (pfunit_sources
  test_CNPhenology.pf)

add_pfunit_ctest(CNPhenology
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
