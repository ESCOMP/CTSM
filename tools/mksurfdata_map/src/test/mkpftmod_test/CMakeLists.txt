set (pfunit_sources
  test_pft_override.pf)

create_pFUnit_test(mkpftMod test_mkpft_exe
  "${pfunit_sources}" "")

target_link_libraries(test_mkpft_exe mksurfdat)
