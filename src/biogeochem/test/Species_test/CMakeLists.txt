set (pfunit_sources
  test_SpeciesNonIsotope.pf
  test_SpeciesIsotope.pf)

add_pfunit_ctest(Species
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES clm csm_share)
