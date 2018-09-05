create_pFUnit_test(balance test_balance_exe
  "test_Balance.pf" "")

target_link_libraries(test_balance_exe clm csm_share esmf_wrf_timemgr)
