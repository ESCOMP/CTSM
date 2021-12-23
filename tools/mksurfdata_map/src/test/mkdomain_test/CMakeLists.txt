set (pfunit_sources
  test_mkdomain.pf)

create_pFUnit_test(mkdomainMod test_mkdomain_exe
  "${pfunit_sources}" "")

target_link_libraries(test_mkdomain_exe mksurfdat)
