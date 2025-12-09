set (pfunit_sources
  test_find_k_max_indices.pf
  test_convert_to_logical.pf
  )

add_pfunit_ctest(array_utils
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
