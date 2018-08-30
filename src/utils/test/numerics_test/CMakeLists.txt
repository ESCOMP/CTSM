set (pfunit_sources
  test_truncate_small_values.pf)

set (extra_sources
  )

create_pFUnit_test(numerics test_numerics_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_numerics_exe clm csm_share esmf_wrf_timemgr)
