set(pfunit_sources
  test_downscale_forcings.pf
  test_partition_precip.pf)

create_pFUnit_test(atm2lnd test_atm2lnd_exe
  "${pfunit_sources}" "")

target_link_libraries(test_atm2lnd_exe clm csm_share esmf_wrf_timemgr)