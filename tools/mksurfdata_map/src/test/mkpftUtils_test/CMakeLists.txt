set (pfunit_sources
  test_adjust_total_veg_area.pf
  test_convert_from_p2g.pf)

create_pFUnit_test(mkpftUtils test_mkpftUtils_exe
  "${pfunit_sources}" "")

target_link_libraries(test_mkpftUtils_exe mksurfdat)