create_pFUnit_test(ctsm_TimeManager test_ctsm_TimeManager_exe
  "test_ctsm_TimeManager.pf" "")

target_link_libraries(test_ctsm_TimeManager_exe clm csm_share esmf_wrf_timemgr)
