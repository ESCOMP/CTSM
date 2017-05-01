create_pFUnit_test(surfrdUtils test_surfrdUtils_exe
  "test_surfrdUtils.pf" "")

target_link_libraries(test_surfrdUtils_exe clm csm_share)
