set (pfunit_sources
  test_dynVarTimeUninterp.pf
  test_dynVarTimeInterp.pf)

# extra sources used for this test, which are not .pf files
set (extra_sources
  test_dynVarShared.F90)

add_pfunit_ctest(dynVar
  TEST_SOURCES "${pfunit_sources}"
  OTHER_SOURCES "${extra_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
