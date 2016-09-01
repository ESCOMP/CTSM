create_pFUnit_test(annual_flux_dribbler test_annual_flux_dribbler_exe
  "test_annual_flux_dribbler.pf" "")

target_link_libraries(test_annual_flux_dribbler_exe clm csm_share esmf_wrf_timemgr)
