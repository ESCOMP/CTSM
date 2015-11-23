create_pFUnit_test(acspinup test_acspinup_exe
  "test_acspinup.pf" "")

target_link_libraries(test_acspinup_exe clm csm_share)
