create_pFUnit_test(clm_time_manager test_clm_time_manager_exe
  "test_clm_time_manager.pf" "")

target_link_libraries(test_clm_time_manager_exe clm csm_share esmf_wrf_timemgr)
