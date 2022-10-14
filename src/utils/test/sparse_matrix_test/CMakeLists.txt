set (pfunit_sources
  test_sparse_matrix.pf)

set (extra_sources
  )

create_pFUnit_test(sparse_matrix test_sparse_matrix_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_sparse_matrix_exe clm csm_share)
