set (pfunit_sources
  test_assertions.pf)

set (extra_sources
  )

create_pFUnit_test(assertions test_assertions_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_assertions_exe clm csm_share esmf_wrf_timemgr)
