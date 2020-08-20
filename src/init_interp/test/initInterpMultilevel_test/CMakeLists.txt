set (pfunit_sources
  ctsm_InitInterpMultilevelMock.pf
  test_init_interp_multilevel_interp.pf
  test_init_interp_multilevel_snow.pf
  test_init_interp_multilevel_split.pf)

# extra sources used for this test, which are not .pf files
set (extra_sources
  ctsm_MultilevelInterpFactory.F90)

create_pFUnit_test(ctsm_InitInterpMultilevel test_ctsm_InitInterpMultilevel_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_ctsm_InitInterpMultilevel_exe clm csm_share)