set (pfunit_sources
  test_matrix.pf)

add_pfunit_ctest(matrix
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share)
