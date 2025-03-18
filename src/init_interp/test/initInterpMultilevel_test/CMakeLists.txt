set (pfunit_sources
  initInterpMultilevelMock.pf
  test_init_interp_multilevel_interp.pf
  test_init_interp_multilevel_snow.pf
  test_init_interp_multilevel_split.pf)

# extra sources used for this test, which are not .pf files
set (extra_sources
  multilevel_interp_factory.F90)

add_pfunit_ctest(initInterpMultilevel
  TEST_SOURCES "${pfunit_sources}"
  OTHER_SOURCES "${extra_sources}"
  LINK_LIBRARIES clm csm_share)
