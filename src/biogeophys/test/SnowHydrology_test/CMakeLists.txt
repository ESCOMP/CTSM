set (pfunit_sources
  test_SnowHydrology_initSnowLayers.pf
  test_SnowHydrology_newSnowBulkDensity.pf
  test_SnowHydrology_SnowCappingExcess.pf)

add_pfunit_ctest(SnowHydrology
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share esmf
  EXTRA_FINALIZE unittest_finalize_esmf
  EXTRA_USE unittestInitializeAndFinalize)
