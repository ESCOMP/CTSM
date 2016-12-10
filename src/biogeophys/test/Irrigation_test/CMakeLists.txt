set (pfunit_sources
  test_irrigation.pf)

# extra sources used for this test, which are not .pf files
# (currently none)
set (extra_sources
  )

create_pFUnit_test(irrigation test_irrigation_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_irrigation_exe clm csm_share esmf_wrf_timemgr)