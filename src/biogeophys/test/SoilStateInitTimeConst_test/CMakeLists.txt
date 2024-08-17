set (pfunit_sources
  test_dust_soil_clay_functions.pf)

add_pfunit_ctest(SoilStateInit
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share)
