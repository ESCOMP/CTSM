set (pfunit_sources
  test_glc_elevclasses_are_same.pf)

create_pFUnit_test(initInterpUtils test_initInterpUtils_exe
  "${pfunit_sources}" "")

target_link_libraries(test_initInterpUtils_exe clm csm_share)