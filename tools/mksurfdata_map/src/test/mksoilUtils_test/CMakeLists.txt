set (pfunit_sources
  test_dominant_soil_color.pf)

create_pFUnit_test(mksoilUtils test_mksoilUtils_exe
  "${pfunit_sources}" "")

target_link_libraries(test_mksoilUtils_exe mksurfdat)