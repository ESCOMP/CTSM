set(pfunit_sources
  test_update_landunit_weights_one_gcell.pf
  test_update_landunit_weights.pf)

add_pfunit_ctest(dynLandunitArea
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share)
