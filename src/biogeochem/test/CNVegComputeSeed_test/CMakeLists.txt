set (pfunit_sources
  test_ComputeSeedAmounts.pf)

add_pfunit_ctest(CNVegComputeSeed
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
