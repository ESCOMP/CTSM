set (pfunit_sources
  test_matrix.pf)

set (extra_sources
  )

create_pFUnit_test(matrix test_matrix_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_matrix_exe clm csm_share)
