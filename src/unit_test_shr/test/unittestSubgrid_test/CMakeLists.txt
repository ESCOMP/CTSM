create_pFUnit_test(unittestSubgrid test_unittestSubgrid_exe
  "test_unittestSubgrid.pf" "")

target_link_libraries(test_unittestSubgrid_exe clm csm_share)