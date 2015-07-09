create_pFUnit_test(ncdio_utils test_ncdio_utils_exe
  "test_ncdio_utils.pf" "")

target_link_libraries(test_ncdio_utils_exe clm csm_share)
