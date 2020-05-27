create_pFUnit_test(initVertical test_initVertical_exe
  "test_initVertical.pf" "")

target_link_libraries(test_initVertical_exe clm csm_share esmf_wrf_timemgr)