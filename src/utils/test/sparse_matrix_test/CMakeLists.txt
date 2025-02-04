set (pfunit_sources
  test_sparse_matrix.pf)

add_pfunit_ctest(sparse_matrix
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share)
