set (pfunit_sources
  test_DustEmisZender2003.pf
  test_DustEmisLeung2023.pf
)

add_pfunit_ctest(DustEmis
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
