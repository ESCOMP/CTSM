set (pfunit_sources
  test_find_k_max_indices.pf)

set (extra_sources
  )

create_pFUnit_test(array_utils test_array_utils_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_array_utils_exe clm csm_share esmf_wrf_timemgr)
