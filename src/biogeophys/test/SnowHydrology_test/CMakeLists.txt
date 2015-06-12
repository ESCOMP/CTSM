create_pFUnit_test(SnowHydrology test_SnowHydrology_exe
  "test_SnowHydrology.pf" "")

target_link_libraries(test_SnowHydrology_exe clm csm_share esmf_wrf_timemgr)