create_pFUnit_test(initInterpMindist test_initInterpMindist_exe
  "test_init_interp_mindist.pf" "")

target_link_libraries(test_initInterpMindist_exe clm csm_share)