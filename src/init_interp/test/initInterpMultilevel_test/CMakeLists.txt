set (pfunit_sources
  initInterpMultilevelMock.pf
  test_init_interp_multilevel_interp.pf
  test_init_interp_multilevel_snow.pf
  test_init_interp_multilevel_split.pf)

# extra sources used for this test, which are not .pf files
set (extra_sources
  multilevel_interp_factory.F90)

create_pFUnit_test(initInterpMultilevel test_initInterpMultilevel_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_initInterpMultilevel_exe clm csm_share)