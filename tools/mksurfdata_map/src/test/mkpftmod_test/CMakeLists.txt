set (pfunit_sources
  test_pftrun.pf
  test_pft_oride.pf
  test_pftInit.pf)

create_pFUnit_test(mkpftMod test_mkpft_exe
  "${pfunit_sources}" "")

target_link_libraries(test_mkpft_exe mksurfdat)
