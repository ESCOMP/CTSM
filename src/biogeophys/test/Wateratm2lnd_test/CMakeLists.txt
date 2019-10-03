set (pfunit_sources
  test_set_tracers.pf)

# extra sources used for this test, which are not .pf files
# (currently none)
set (extra_sources
  )

create_pFUnit_test(water_atm2lnd test_water_atm2lnd_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_water_atm2lnd_exe clm csm_share esmf_wrf_timemgr)