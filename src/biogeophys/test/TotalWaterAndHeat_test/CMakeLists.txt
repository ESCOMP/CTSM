set (pfunit_sources
  test_total_water_and_heat.pf)

# extra sources used for this test, which are not .pf files
# (currently none)
set (extra_sources
  )

create_pFUnit_test(total_water_and_heat test_total_water_and_heat_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_total_water_and_heat_exe clm csm_share esmf_wrf_timemgr)